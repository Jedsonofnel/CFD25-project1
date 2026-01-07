//
// CFD Project 1
// Author: Jed Nelson <jed@nelson.ac>
//

package main

import (
	"encoding/csv"
	"flag"
	"io"
	"log"
	"math"
	"os"
	"strconv"
)

//
// Constants
//

const (
	// physical parameters
	RHO    float64 = 0.5 // density
	GAMMA  float64 = 0.5 // diffusivity
	LENGTH float64 = 1

	// boundaries
	PHI_0 float64 = 100 // "leftmost"
	PHI_L float64 = 20  // "rightmost"

	// sentinel values
	LEFT_BOUND_MARKER  int32 = -1
	RIGHT_BOUND_MARKER int32 = -2
)

//
// Main entrypoint and core functionality (ie look here first)
//

func main() {
	// Parse command-line flags
	numCells := flag.Int("n", 10, "number of cells")
	velocity := flag.Float64("u", 5.0, "velocity")
	flag.Parse()

	mesh := New1DMesh(LENGTH, *numCells)
	analytic := CalculateAnalytical(mesh, *velocity)
	numericCDS := CalculateNumeric(mesh, *velocity, DivOperatorCDS)
	numericUDS := CalculateNumeric(mesh, *velocity, DivOperatorUDS)
	numericPLDS := CalculateNumericPLDS(mesh, *velocity)

	WriteToCSV(os.Stdout, mesh, analytic, numericCDS, numericUDS, numericPLDS)
}

func CalculateAnalytical(mesh *Mesh, velocity float64) []float64 {
	field := make([]float64, mesh.numCells)

	// peclet number
	Pe := RHO * velocity * LENGTH / GAMMA

	for i := range field {
		x := mesh.centroids[i]

		if math.Abs(Pe) < 1e-10 { // pure diffusion therefore linear profile
			field[i] = PHI_0 + (x/LENGTH)*(PHI_L-PHI_0)
		} else {
			coeff := (math.Exp(Pe*x/LENGTH) - 1) / (math.Exp(Pe) - 1)
			field[i] = PHI_0 + coeff*(PHI_L-PHI_0)
		}
	}

	return field
}

func CalculateNumeric(mesh *Mesh, velocity float64, divOp DivOp) []float64 {
	field := make([]float64, mesh.numCells)
	system := NewFVSystem(mesh)

	// 1) apply laplacian operator
	LaplacianOperator(mesh, system.matrix, GAMMA)

	// 2) apply divergence operator (using function pointer parameter)
	divOp(mesh, system.matrix, RHO, velocity)

	// 3) apply boundary conditions (they contribute to diagonal and RHS)
	DirichletBC(system, PHI_0, LEFT_BOUND_MARKER)
	DirichletBC(system, PHI_L, RIGHT_BOUND_MARKER)

	// 4) solve the matrix
	system.SolveTDMA(field)
	// system.SolveCG(field, 1e-6, 50) // (alternative - more general than TDMA)

	return field
}

func CalculateNumericPLDS(mesh *Mesh, velocity float64) []float64 {
	field := make([]float64, mesh.numCells)
	system := NewFVSystem(mesh)

	// PLDSOperator combines laplacian and divergence
	PLDSOperator(mesh, system.matrix, RHO, velocity, GAMMA)
	DirichletBC(system, PHI_0, LEFT_BOUND_MARKER)
	DirichletBC(system, PHI_L, RIGHT_BOUND_MARKER)

	system.SolveTDMA(field)

	return field
}

//
// Meshing and geometry
//

// Each connection corresponds to a face connecting two centroids.  These are
// demarcated as owner/neighbour for distinguishing flux signs.
type Connection struct {
	owner     int32
	neighbour int32
}

type Mesh struct {
	// useful metadata
	numCells int
	numConns int

	// indexed by cell
	centroids []float64
	cellVols  []float64 // in 1D these are really lengths

	// indexed by face/connection (same thing really)
	connections []Connection
	connDists   []float64
	faceAreas   []float64 // in 1D this is moot but may as well be consistent with literature
	faceNormals []float64 // 1D (either 1 or -1 but kept as a float for consistency)
}

// Creates a 1D mesh WITHOUT half cells at boundaries
func New1DMesh(length float64, numCells int) *Mesh {
	// cellwise data
	dx := length / float64(numCells)
	centroids := make([]float64, numCells)
	cellVols := make([]float64, numCells)
	for i := range numCells {
		centroids[i] = (0.5 + float64(i)) * dx
		cellVols[i] = dx
	}

	// connectionwise data
	numInternalConnections := numCells - 1
	totalConnections := numInternalConnections + 2 // ie include boundary connections
	connections := make([]Connection, totalConnections)
	connDists := make([]float64, totalConnections)
	faceAreas := make([]float64, totalConnections)
	faceNormals := make([]float64, totalConnections)

	// internal connections (these are all left->right)
	for i := range numInternalConnections {
		connections[i+1] = Connection{owner: int32(i), neighbour: int32(i + 1)}
		connDists[i+1] = dx
		faceAreas[i+1] = 1
		faceNormals[i+1] = 1.0
	}

	// leftmost boundary
	connections[0] = Connection{owner: 0, neighbour: LEFT_BOUND_MARKER}
	connDists[0] = dx / 2 // half distance because of boundary
	faceAreas[0] = 1
	faceNormals[0] = -1 // negative normal (right->left)

	// rightmost boundary
	rmIdx := totalConnections - 1
	connections[rmIdx] = Connection{owner: int32(numCells - 1), neighbour: RIGHT_BOUND_MARKER}
	connDists[rmIdx] = dx / 2
	faceAreas[rmIdx] = 1
	faceNormals[rmIdx] = 1

	return &Mesh{
		numCells:    numCells,
		numConns:    totalConnections,
		centroids:   centroids,
		cellVols:    cellVols,
		connections: connections,
		connDists:   connDists,
		faceAreas:   faceAreas,
		faceNormals: faceNormals,
	}
}

//
// Output helpers
//

// Writes a CSV to the output writer supplied
func WriteToCSV(
	out io.Writer,
	mesh *Mesh,
	analytic []float64,
	numericCDS []float64,
	numericUDS []float64,
	numericPLDS []float64,
) {
	headers := []string{"x", "phi_analytical", "phi_CDS", "phi_UDS", "phi_PLDS"}

	w := csv.NewWriter(out)
	err := w.Write(headers)
	if err != nil {
		log.Fatalln("error writing record to csv:", err)
	}

	// manually add boundary condition row for plotting
	phi0String := strconv.FormatFloat(PHI_0, 'g', -1, 64)
	err = w.Write([]string{"0.0", phi0String, phi0String, phi0String, phi0String})
	if err != nil {
		log.Fatalln("error writing record to csv:", err)
	}

	for i := range mesh.numCells {
		record := []string{
			strconv.FormatFloat(mesh.centroids[i], 'g', -1, 64),
			strconv.FormatFloat(analytic[i], 'g', -1, 64),
			strconv.FormatFloat(numericCDS[i], 'g', -1, 64),
			strconv.FormatFloat(numericUDS[i], 'g', -1, 64),
			strconv.FormatFloat(numericPLDS[i], 'g', -1, 64),
		}

		err := w.Write(record)
		if err != nil {
			log.Fatalln("error writing record to csv:", err)
		}
	}

	// right boundary condition
	lengthString := strconv.FormatFloat(LENGTH, 'g', -1, 64)
	phiLString := strconv.FormatFloat(PHI_L, 'g', -1, 64)
	err = w.Write([]string{lengthString, phiLString, phiLString, phiLString, phiLString})
	if err != nil {
		log.Fatalln("error writing record to csv:", err)
	}

	w.Flush()

	if err := w.Error(); err != nil {
		log.Fatal(err)
	}
}
