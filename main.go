//
// CFD Project 1
// Author: Jed Nelson <jed@nelson.ac>
//

package main

import (
	"encoding/csv"
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
	BOUNDARY_INDEX int = -1
)

//
// Main entrypoint
//

// TODO get flags to work for velocity and number of cells (ostensibly dx)

func main() {
	numCells := 10
	var velocity float64 = 5
	mesh := New1DMesh(LENGTH, numCells)
	analytic := CalculateAnalytical(mesh, velocity)

	WriteToCSV(os.Stdout, mesh, analytic)
}

//
// Data structures
//

// Based on the following
// https://openfoamwiki.net/index.php/OpenFOAM_guide/Matrices_in_OpenFOAM
type LDUMatrix struct {
	diag  []float64 // indexed by cell
	lower []float64 // indexed by face
	upper []float64 // indexed by face
}

// Each connection corresponds to a face connecting two centroids.  These are
// demarcated as owner/neighbour for distinguishing flux signs.
type Connection struct {
	owner     int
	neighbour int
}

type Mesh struct {
	numCells int

	// indexed by cell
	centroids []float64
	cellVols  []float64 // in 1D these are really lengths

	// indexed by face/connection (same thing really)
	connections []Connection
	connDists   []float64
	faceAreas   []float64 // in 1D this is moot but may as well be consistent with literature
}

//
// "Meshing"
//

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

	// internal connections
	for i := range numInternalConnections {
		connections[i+1] = Connection{owner: i, neighbour: i + 1}
		connDists[i+1] = dx
		faceAreas[i+1] = 1
	}

	// leftmost boundary
	connections[0] = Connection{owner: 0, neighbour: BOUNDARY_INDEX}
	connDists[0] = dx / 2 // half distance because of boundary
	faceAreas[0] = 1

	// rightmost boundary
	rmIdx := totalConnections - 1
	connections[rmIdx] = Connection{owner: 0, neighbour: BOUNDARY_INDEX}
	connDists[rmIdx] = dx / 2
	faceAreas[rmIdx] = 1

	return &Mesh{
		numCells:    numCells,
		centroids:   centroids,
		cellVols:    cellVols,
		connections: connections,
		connDists:   connDists,
		faceAreas:   faceAreas,
	}
}

//
// Analytical solution
//

func CalculateAnalytical(mesh *Mesh, velocity float64) []float64 {
	field := make([]float64, mesh.numCells)

	for i := range field {
		x := mesh.centroids[i]
		coeff := (math.Exp(RHO*velocity*x/GAMMA) - 1) /
			(math.Exp(RHO*velocity*LENGTH/GAMMA) - 1)
		field[i] = PHI_0 + coeff*(PHI_L-PHI_0)
	}

	return field
}

//
// Output helpers
//

// Writes a CSV to the writer supplied
func WriteToCSV(
	out io.Writer,
	mesh *Mesh,
	analytic []float64,
) {
	headers := []string{"x", "phi_analytical"}

	w := csv.NewWriter(out)
	err := w.Write(headers)
	if err != nil {
		log.Fatalln("error writing record to csv:", err)
	}

	for i := range mesh.numCells {
		record := []string{
			strconv.FormatFloat(mesh.centroids[i], 'g', -1, 64),
			strconv.FormatFloat(analytic[i], 'g', -1, 64),
		}

		err := w.Write(record)
		if err != nil {
			log.Fatalln("error writing record to csv:", err)
		}
	}

	w.Flush()

	if err := w.Error(); err != nil {
		log.Fatal(err)
	}
}
