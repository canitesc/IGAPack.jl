#script for solving the acoustic duct problem
# -a0(x)*\Delta u(x,y)+a1(x)*u(x,y)=f(x,y) for (x,y) ∈ Ω:= (0,2)×(0,1)
# Here a0(x,y)=-1, a1(x,y)=k^2, f(x,y)=0
# with Neumann boundary conditions ∂u/∂x(x,y)=cos(mπy) for x = 0 and
# Robin boundary conditions ∂u/∂x(x,y) + iku = 0 for x=2

include("../src/Solver2D.jl")
include("../src/SplinePlotting.jl")

a0(x,y) = 1
a1(x,y) = -k^2
f(x,y) = 0
k = 40
m = 2
kx = sqrt(k^2 - (m*π)^2)
LHS = [im*kx -im*kx; (k-kx)*exp(-2*im*kx) (k+kx)*exp(2*im*kx)]
RHS = [1, 0]
A = LHS\RHS

exact_sol(x,y) = cos(m*π*y)*(A[1]*exp(-im*kx*x)+A[2]*exp(im*kx*x))
deriv_exact_sol(x,y) = [cos(m*π*y)*(A[1]*(-im)*kx*exp(-im*kx*x)+A[2]*im*kx*exp(im*kx*x)),
                        -sin(m*π*y)*(m*π)*(A[1]*exp(-im*kx*x)+A[2]*exp(im*kx*x))]

#Define the boundary conditions
u_bound_left(x,y) = cos(m*π*y)
u_bound_right(x,y) = 0.
bound_left = Boundary2D("Neumann", "Left", u_bound_left)
bound_right = Boundary2D("Robin", "Right", -im*k, u_bound_right)
bound_all = [bound_left, bound_right]

#Define the domain geometry
cornerLowerLeft = [0., 0.]
lengthx = 2.
lengthy = 1.
degP = [3, 3]
numSubdiv = 80
nrb = nrbsquare(cornerLowerLeft, lengthx, lengthy, degP, numSubdiv)
IGAmesh = genMesh(nrb)
#plotBasisParam(IGAmesh)

gauss_rule = [genGaussLegendre(degP[1]+1), genGaussLegendre(degP[2]+1)]
stiff = assemble_stiff2D(IGAmesh, a0, gauss_rule)
mass = assemble_mass2D(IGAmesh, a1, gauss_rule)
rhs = assemble_rhs2D(IGAmesh, f, gauss_rule)
bcdof_all, elem_all = classifyBoundary2D(IGAmesh)
lhs = mass + stiff
lhs, rhs = applyBCnurbs(IGAmesh, bound_all, lhs, rhs, bcdof_all, elem_all, gauss_rule)

println("Solving linear system")
@time sol0 = lhs\rhs
println("Plotting...")
@time plotSol(IGAmesh, real(sol0), "AcousticDuctReal")
@time plotSolError(IGAmesh, real(sol0), exact_sol, "AcousticDuctRealError")
println("Computing error")
@time relL2Err, relH1Err = compErrorNorm(IGAmesh, real(sol0), exact_sol, deriv_exact_sol, a0, gauss_rule)
println("Relative L2-norm error is $relL2Err")
println("Relative H1-seminorm error is $relH1Err")
print("Done!")
