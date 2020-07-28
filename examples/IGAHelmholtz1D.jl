#script for solving a 1D BVP of the form
# -a0(x)*u''(x)+a1(x)*u(x)=f(x) for x ∈ (0,1)
# Here a0(x)=1, a1(x)=0, f(x)=4*π^2*sin(2*π*x)
# with Neumann-Robin boundary conditions u'(0) = ik, u'(1)-i*k*u(1)=0
include("../src/Solver1D.jl")

k=40
a0(x) = -1
a1(x) = k^2
f(x) = 0
exact_sol(x) = exp(im*k*x)

ptLeft = 0.
ptRight = 1.
uPtLeft = 0.
uPtRight = 1.
nrb = nrbline(ptLeft, ptRight)

bound_left = Boundary1D("Neumann", ptLeft, uPtLeft, -im*k)
α = -im*k
bound_right = Boundary1D("Robin", ptRight, uPtRight, α, 0)

numElem = 40
degP = 3
new_knots = collect(LinRange(ptLeft, ptRight, numElem+1)[2:end-1])
nrb = nrbdegelev(nrb, [degP-1])
nrb = nrbkntins(nrb, [new_knots])

prob_sin = Problem1D(f, a0, a1, [bound_left, bound_right], nrb)
IGAmesh = genMesh(prob_sin.domain)
plotBasisParam(IGAmesh)

gauss_rule = genGaussLegendre(degP+1)
stiff = assemble_stiff(IGAmesh, a0, gauss_rule)
mass = assemble_mass(IGAmesh, a1, gauss_rule)
rhs = assemble_rhs(IGAmesh, f, gauss_rule)
lhs, rhs = applyBCnurbs(prob_sin.boundary_cond, stiff, mass, rhs, nrb)
sol0 = lhs\rhs
plotSol(IGAmesh, real(sol0), "Helmholtz1D")
plotSolError(IGAmesh, real(sol0), exact_sol, "Helmholtz1D")
