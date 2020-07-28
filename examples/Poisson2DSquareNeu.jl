#script for solving a 2D BVP of the form
# -a0(x)*\Delta u(x,y)+a1(x)*u(x,y)=f(x,y) for (x,y) ∈ Ω:= (0,1)^2
# Here a0(x,y)=1, a1(x,y)=0, f(x,y)=-8*π^2*sin(2*π*x)*cos(2*π*y)
# with Dirichlet boundary conditions u(x,y)=0 for x=0 and
# Neumann boundary conditions du/dx = 2*π*cos(2*π*y)
include("../src/Solver2D.jl")
include("../src/SplinePlotting.jl")

a0(x,y) = 1
a1(x,y) = 0
f(x,y) = 8*π^2*sin(2*π*x)*cos(2*π*y)
exact_sol(x,y) = sin(2*π*x)*cos(2*π*y)
deriv_exact_sol(x,y) = [2*π*cos(2*π*x)*cos(2*π*y), -2*π*sin(2*π*x)*sin(2*π*y)]

#Define the boundary conditions
u_bound_dir(x,y) = 0 # for x = 0
u_bound_neu(x,y) = 2*π*cos(2*π*y)
bound_left = Boundary2D("Dirichlet", "Left", u_bound_dir)
bound_right = Boundary2D("Neumann", "Right", u_bound_neu)
bound_all = [bound_left, bound_right]

#Define the domain geometry
cornerLowerLeft = [0., 0.]
lengthx = 1.
lengthy = 1.
degP = [3, 3]
numSubdiv = 20
nrb = nrbsquare(cornerLowerLeft, lengthx, lengthy, degP, numSubdiv)
#Plots.surface(reuse = false)
#plt = nrbctrlplot(nrb)
#plt = nrbkntplot(nrb)
#display(plt)
IEN, elemVertex = makeIEN(nrb.knots, numSubdiv^2, degP)
IGAmesh = genMesh(nrb)
#plotBasisParam(IGAmesh)

gauss_rule = [genGaussLegendre(degP[1]+1), genGaussLegendre(degP[2]+1)]
stiff = assemble_stiff2D(IGAmesh, a0, gauss_rule)
rhs = assemble_rhs2D(IGAmesh, f, gauss_rule)
bcdof_all, elem_all = classifyBoundary2D(IGAmesh)
lhs, rhs = applyBCnurbs(IGAmesh, bound_all, stiff, rhs, bcdof_all, elem_all, gauss_rule)

@time sol0 = lhs\rhs
plotSol(IGAmesh, sol0, "Poisson2DSquare")
plotSolError(IGAmesh, sol0, exact_sol, "Poisson2DSquareError")
relL2Err, relH1Err = compErrorNorm(IGAmesh, sol0, exact_sol, deriv_exact_sol, a0, gauss_rule)
println("Relative L2-norm error is $relL2Err")
println("Relative energy-norm error is $relH1Err")
print("Done!")
