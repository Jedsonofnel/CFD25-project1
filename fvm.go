package main

import (
	"math"
)

//
// Data structures
//

// Based on the following
// https://openfoamwiki.net/index.php/OpenFOAM_guide/Matrices_in_OpenFOAM
type LDUMatrix struct {
	diag []float64 // indexed by cell

	// both indexed by face/connection (therefore some are internal faces and some are boundary)
	lower []float64
	upper []float64

	// copied slice header from mesh - required to index lower/upper correctly
	conns []Connection
}

func NewLDUMatrix(mesh *Mesh) *LDUMatrix {
	return &LDUMatrix{
		diag:  make([]float64, mesh.numCells),
		lower: make([]float64, mesh.numConns),
		upper: make([]float64, mesh.numConns),
		conns: mesh.connections,
	}
}

type FVSystem struct {
	matrix *LDUMatrix
	rhs    []float64
}

func NewFVSystem(mesh *Mesh) *FVSystem {
	matrix := NewLDUMatrix(mesh)
	return &FVSystem{
		matrix: matrix,
		rhs:    make([]float64, mesh.numCells),
	}
}

//
// Operators. These are simple data transformations of the system (no state)
//

// LaplacianOperator transforms the matrix according to the FVM discretisation
// of gamma d2x/dy2
func LaplacianOperator(mesh *Mesh, matrix *LDUMatrix, gamma float64) {
	for i, conn := range mesh.connections {
		faceArea := mesh.faceAreas[i]
		distance := mesh.connDists[i]
		fluxCoeff := gamma * faceArea / distance

		// crucially we are not SETTING the values just mutating them so this
		// operator can be part of a more complex transformation pipeline
		// (multiple operators and source terms etc)
		matrix.lower[i] -= fluxCoeff
		matrix.upper[i] -= fluxCoeff

		// only add to diagonals for internal faces (let boundaries handle
		// their own diags)
		if conn.neighbour >= 0 {
			matrix.diag[conn.owner] += fluxCoeff
			matrix.diag[conn.neighbour] += fluxCoeff
		}
	}
}

// DivOp is a function pointer type with the signature corresponding to a
// divergence operator tranformation function
type DivOp func(mesh *Mesh, matrix *LDUMatrix, rho, U float64)

func DivOperatorCDS(mesh *Mesh, matrix *LDUMatrix, rho, U float64) {
	for i, conn := range mesh.connections {
		faceArea := mesh.faceAreas[i]
		faceNormal := mesh.faceNormals[i]
		F := rho * U * faceArea * faceNormal

		matrix.lower[i] += F / 2
		matrix.upper[i] -= F / 2

		if conn.neighbour >= 0 {
			matrix.diag[conn.owner] += F / 2
			matrix.diag[conn.neighbour] -= F / 2
		}
	}
}

func DivOperatorUDS(mesh *Mesh, matrix *LDUMatrix, rho, U float64) {
	for i, conn := range mesh.connections {
		faceArea := mesh.faceAreas[i]
		faceNormal := mesh.faceNormals[i]
		F := rho * U * faceArea * faceNormal

		matrix.lower[i] -= max(-F, 0)
		matrix.upper[i] -= max(F, 0)

		if conn.neighbour >= 0 {
			matrix.diag[conn.owner] += max(F, 0)
			matrix.diag[conn.neighbour] += max(-F, 0)
		}
	}
}

// PLDSOperator implements the power law scheme which fundamentally blends
// diffusion with advection according to the Peclet number and so you cannot
// separate the diffusion and advection operators (it's its own sort of thing)
func PLDSOperator(mesh *Mesh, matrix *LDUMatrix, rho, U, gamma float64) {
	for i, conn := range mesh.connections {
		faceArea := mesh.faceAreas[i]
		faceNormal := mesh.faceNormals[i]
		distance := mesh.connDists[i]

		D := gamma * faceArea / distance
		F := rho * U * faceArea * faceNormal

		Pe := F / D

		// Power law function: A(|Pe|) = max(0, (1 - 0.1|Pe|)^5)
		absPe := math.Abs(Pe)
		A := max(0, math.Pow(1-0.1*absPe, 5))

		// Patankar eq5.34 page 91
		matrix.lower[i] -= D*A + max(-F, 0)
		matrix.upper[i] -= D*A + max(F, 0)

		// Only add to diagonals for internal faces
		if conn.neighbour >= 0 {
			matrix.diag[conn.owner] += D*A + max(F, 0)
			matrix.diag[conn.neighbour] += D*A + max(-F, 0)
		}
	}
}

// DirichletBC applies the diag and RHS contributions for the FVM
// discretisation to all connections where neighbour matches the marker
func DirichletBC(system *FVSystem, value float64, marker int32) {
	matrix := system.matrix
	for i, conn := range matrix.conns {
		if conn.neighbour == marker {
			bcCoeff := -matrix.lower[i]
			diagCoeff := -matrix.upper[i]

			matrix.diag[conn.owner] += diagCoeff
			system.rhs[conn.owner] += value * bcCoeff // phi_bc * a_n
		}
	}
}

//
// Solving
//

// SolveTDMA solves Ax = b with Thomas' algorithm and mutates the x vector given as an argument
func (sys *FVSystem) SolveTDMA(x []float64) {
	n := len(sys.matrix.diag)

	// extract tridiagonal structure from connection-based LDU
	sub := make([]float64, n)
	super := make([]float64, n)

	for i, conn := range sys.matrix.conns {
		if conn.neighbour >= 0 { // internal connection only
			owner := int(conn.owner)
			neighbour := int(conn.neighbour)

			// lower[i] = coefficient of phi[neighbour] in owner's equation
			// If neighbour = owner+1, this is the super-diagonal
			if neighbour == owner+1 {
				super[owner] = sys.matrix.lower[i]
			}

			// upper[i] = coefficient of phi[owner] in neighbour's equation
			// If owner = neighbour-1, this is the sub-diagonal
			if owner == neighbour-1 {
				sub[neighbour] = sys.matrix.upper[i]
			}
		}
	}

	// Scratch space for forward elimination
	cp := make([]float64, n) // modified super-diagonal
	dp := make([]float64, n) // modified RHS

	// Forward elimination
	cp[0] = super[0] / sys.matrix.diag[0]
	dp[0] = sys.rhs[0] / sys.matrix.diag[0]

	for i := 1; i < n-1; i++ {
		denom := sys.matrix.diag[i] - sub[i]*cp[i-1]
		cp[i] = super[i] / denom
		dp[i] = (sys.rhs[i] - sub[i]*dp[i-1]) / denom
	}

	// Last row (no super-diagonal term)
	denom := sys.matrix.diag[n-1] - sub[n-1]*cp[n-2]
	dp[n-1] = (sys.rhs[n-1] - sub[n-1]*dp[n-2]) / denom

	// Back substitution
	x[n-1] = dp[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}
}

// Wipe zeros all matrix coefficients
func (m *LDUMatrix) Wipe() {
	for i := range m.diag {
		m.diag[i] = 0
	}
	for i := range m.lower {
		m.lower[i] = 0
		m.upper[i] = 0
	}
}
