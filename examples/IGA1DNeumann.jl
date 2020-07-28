#script for solving a 1D BVP of the form
# -a0(x)*u''(x)+a1(x)*u(x)=f(x) for x ∈ (0,1)
# Here a0(x)=1, a1(x)=0, f(x)=4*π^2*sin(2*π*x)
# with Dirichlet boundary conditions u(0) = 0,
# and Neumann boundary condition u'(1)=0
include("../src/Solver1D.jl")

a0(x) = 1
a1(x) = 0
f(x) = (4*π^2)*sin(2*π*x)
exact_sol(x) = sin(2*π*x)

ptLeft = 0.
ptRight = 1.
uPtLeft = 0.
uPtRight = 1.
nrb = nrbline(ptLeft, ptRight)

solPtLeft = exact_sol(ptLeft)
#solPtRight = exact_sol(ptRight)

bound_left = Boundary1D("Dirichlet", ptLeft, uPtLeft, solPtLeft)
bound_right = Boundary1D("Neumann", ptRight, uPtRight, 2*π)

numElem = 10
degP = 2
new_knots = collect(LinRange(ptLeft, ptRight, numElem+1)[2:end-1])
nrb = nrbdegelev(nrb, [degP-1])
nrb = nrbkntins(nrb, [new_knots])

#refine some more knots at the beginning
ptMid = (ptLeft+ptRight)/2
newer_knots = collect(LinRange(ptLeft, ptMid, numElem+1))[2:end-1]
newer_knots = setdifftol(newer_knots, new_knots)
nrb = nrbkntins(nrb, [newer_knots])

prob_sin = Problem1D(f, a0, a1, [bound_left, bound_right], nrb)
IGAmesh = genMesh(prob_sin.domain)
plotBasisParam(IGAmesh)

gauss_rule = genGaussLegendre(degP+1)
stiff = assemble_stiff(IGAmesh, a0, gauss_rule)
mass = assemble_mass(IGAmesh, a1, gauss_rule)
rhs = assemble_rhs(IGAmesh, f, gauss_rule)
lhs, rhs = applyBCnurbs(prob_sin.boundary_cond, stiff, mass, rhs, nrb)
sol0 = lhs\rhs
plotSol(IGAmesh, real(sol0), "Poisson1D")
plotSolError(IGAmesh, real(sol0), exact_sol, "Poisson1D")
