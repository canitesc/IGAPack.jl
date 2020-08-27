#script for the 2D Linear elasticity problem
include("../src/Solver2D.jl")
include("../src/SplinePlotting.jl")

#Material properties
Emod = 1e5
nu = 0.3
Cmat = Emod/(1-nu^2)*[1  nu  0; nu  1  0; 0  0  (1-nu)/2]
material = Material2D(Emod, nu, Cmat)

#symmetry Dirichlet B.C., u_y = 0 for y=0 and u_x=0 for x=0
u_bound_dir_symy0(x,y) = [undef, 0]
u_bound_dir_symx0(x,y) = [0, undef]
# pressure Neumann B.C. τ(x,y) = [x, y] on the circular boundary
bound_press = 10
u_bound_neu(x,y) = bound_press*[x, y]

bound_down = Boundary2D("Dirichlet", "Down", u_bound_dir_symy0)
bound_up = Boundary2D("Dirichlet", "Up", u_bound_dir_symx0)
bound_left = Boundary2D("Neumann", "Left", u_bound_neu)
bound_all = [bound_down, bound_up, bound_left]

#define the domain geometry
rad_int = 1
rad_ext = 4
center = [0., 0.]
sang = 0.
eang = π/2

left_crv =  nrbcirc(rad_int, center, sang, eang)
right_crv = nrbcirc(rad_ext, center, sang, eang)
down_crv = nrbline([rad_int, 0.], [rad_ext, 0.])
up_crv = nrbline([0., rad_int], [0., rad_ext])
nrb = nrbcoons(down_crv, up_crv, left_crv, right_crv)

#degree elevate and refine
degP = [2, 2]
num_subdiv = [20, 20]
new_knots_u = LinRange(0, 1, num_subdiv[1])[2:end-1]
new_knots_v = LinRange(0, 1, num_subdiv[2])[2:end-1]
nrb = nrbdegelev(nrb, degP.+1-nrb.order)
nrb = nrbkntins(nrb, [collect(new_knots_u), collect(new_knots_v)])
# Plots.surface(reuse = false)
# plt = nrbctrlplot(nrb)
# display(plt)

#define the exact solution
function exact_disp(x,y)
    cart2pol(x, y) = hypot(x, y), atan(y, x)
    r, th = cart2pol(x, y)
    u_r = rad_int^2*bound_press*r/(Emod*(rad_ext^2-rad_int^2))*(1-nu+(rad_ext/r)^2*(1+nu))
    return [u_r*cos(th), u_r*sin(th)]
end

function exact_stress(x,y)
    cart2pol(x, y) = hypot(x, y), atan(y, x)
    r, th = cart2pol(x, y)
    sigma_rr = rad_int^2*bound_press/(rad_ext^2-rad_int^2)*(1-rad_ext^2/r^2) 
    sigma_tt = rad_int^2*bound_press/(rad_ext^2-rad_int^2)*(1+rad_ext^2/r^2)
    sigma_rt = 0

    A = [cos(th)^2 sin(th)^2 2*sin(th)*cos(th); sin(th)^2 cos(th)^2 -2*sin(th)*cos(th);
        -sin(th)*cos(th) sin(th)*cos(th) cos(th)^2-sin(th)^2]
    stress = A\[sigma_rr;sigma_tt;sigma_rt]
    return stress
end

#generate the NURBS mesh
IGAmesh = genMesh(nrb)

gauss_rule = [genGaussLegendre(degP[1]+1), genGaussLegendre(degP[2]+1)]
stiff = assemble_stiff_elast2D(IGAmesh, Cmat, gauss_rule)
#Assume no body forces  
rhs = zeros(2*IGAmesh.numBasis)
bcdof_all, elem_all = classifyBoundary2D(IGAmesh)
lhs, rhs = applyBCnurbsElast(IGAmesh, bound_all, stiff, rhs, bcdof_all, elem_all, gauss_rule)
println("Solving linear system")
@time sol0 = lhs\rhs
@time plotSolElast(IGAmesh, Cmat, sol0, "QuarterAnnulus")
@time plotSolErrorElast(IGAmesh, Cmat, sol0, exact_disp, exact_stress, "QuarterAnnulusErr")
@time relL2Err, relH1Err = compErrorNormElast(IGAmesh, Cmat, sol0, exact_disp, exact_stress, gauss_rule)
println("Relative L2-norm error is $relL2Err")
println("Relative energy-norm error is $relH1Err")

