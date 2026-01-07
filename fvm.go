package main

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
	// TODO
}

func DivOperatorUDS(mesh *Mesh, matrix *LDUMatrix, rho, U float64) {
	// TODO
}

func DivOperatorPLDS(mesh *Mesh, matrix *LDUMatrix, rho, U float64) {
	// TODO
}

//
// BCs. These are treated like source terms that decorate the system (rather
// than stateful)
//

// DirichletBC applies the diag and RHS contributions for the FVM
// discretisation to all connections where neighbour matches the marker
func DirichletBC(system *FVSystem, value float64, marker int32) {
	matrix := system.matrix
	for i, conn := range matrix.conns {
		if conn.neighbour == marker {
			fluxCoeff := -matrix.lower[i] // a_n
			matrix.diag[conn.owner] += fluxCoeff
			system.rhs[conn.owner] += value * fluxCoeff // phi_bc * a_n
		}
	}
}

// NeumannFluxBC applies the flux source to the rhs inside system, is not used
// for this project but wanted to demonstrate how this "simple data
// transformation" architecture could be extended.  It's neat because it's not
// OOP.
func NeumannFluxBC(system *FVSystem, mesh *Mesh, marker int32, flux float64) {
	for i, conn := range system.matrix.conns {
		if conn.neighbour == marker {
			faceArea := mesh.faceAreas[i]
			system.rhs[conn.owner] += flux * faceArea
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

// SolveCG solves a system without the tridiagonal assumption - uses conjugate
// gradient iterative approach without a preconditioner.  Done for a laugh.
// See: https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
func (sys *FVSystem) SolveCG(x []float64, tolerance float64, maxIters int) {
	n := len(sys.matrix.diag)
	r := make([]float64, n)
	d := make([]float64, n)
	Ad := make([]float64, n)

	A := sys.matrix
	b := sys.rhs

	Ax := A.MatVec(x, r)
	for i, val := range Ax {
		r[i] = b[i] - val
		d[i] = r[i]
	}

	var rDotr float64 = 0
	for _, val := range r {
		rDotr += val * val
	}

	recomputeAxInterval := 50

	threshold := tolerance * tolerance * rDotr
	for iter := 0; iter < maxIters && rDotr > threshold; iter++ {
		Ad = A.MatVec(d, Ad)

		var dDotAd float64 = 0
		for i, val := range d {
			dDotAd += val * Ad[i]
		}

		alpha := rDotr / dDotAd

		for i, val := range x {
			x[i] = val + alpha*d[i]
		}

		if iter%recomputeAxInterval == 0 {
			Ax = A.MatVec(x, r)
			for i, val := range Ax {
				r[i] = b[i] - val
			}
		} else {
			for i, val := range r {
				r[i] = val - alpha*Ad[i]
			}
		}

		rDotrOld := rDotr
		rDotr = 0
		for _, val := range r {
			rDotr += val * val
		}

		beta := rDotr / rDotrOld

		for i, val := range d {
			d[i] = r[i] + beta*val
		}
	}
}

//
// Matrix helpers
//

// MatVec mutates y such that y = Ax, returns the slice header y
func (m *LDUMatrix) MatVec(x, y []float64) []float64 {
	// zero output
	for i := range y {
		y[i] = 0
	}

	// diagonal contribution
	for i := range m.diag {
		y[i] += m.diag[i] * x[i]
	}

	// off-diagonal contributions from internal faces only
	for i, conn := range m.conns {
		if conn.neighbour >= 0 { // internal only
			y[conn.owner] += m.lower[i] * x[conn.neighbour]
			y[conn.neighbour] += m.upper[i] * x[conn.owner]
		}
	}

	return y
}
