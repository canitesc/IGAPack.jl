#script for the 2D Linear elasticity problem
#plate with a hole under tension at infinity benchmark
include("../src/Solver2D.jl")
include("../src/SplinePlotting.jl")

#Material properties
Emod = 1e5
nu = 0.3
Cmat = Emod/(1-nu^2)*[1  nu  0; nu  1  0; 0  0  (1-nu)/2]
material = Material2D(Emod, nu, Cmat)

#define the domain geometry
R = 1.
L = 4.
center = [0., 0.]
sang = π/2
eang = π

bottom_crv =  nrbcirc(R, center, sang, eang)
bottom_crv = nrbreverse(bottom_crv)
right_seg = nrbline([0, R], [0, L])
left_seg = nrbline([-R, 0], [-L, 0])
up_seg = nrbline([-L, 0], [-L, L], [0, L])
nrb = nrbcoons(bottom_crv, up_seg, left_seg, right_seg)

#degree elevate and refine
degP = [3, 3]
num_subdiv = [10, 20]
nrb = nrbdegelev(nrb, degP.+1-nrb.order)
_, _, new_knots = kntrefine(nrb.knots, num_subdiv, degP, degP.-1)
nrb = nrbkntins(nrb, new_knots)
# Plots.surface(reuse = false)
# plt = nrbctrlplot(nrb)
# display(plt)

#symmetry Dirichlet B.C., u_y = 0 for y=0 and u_x=0 for x=0
u_bound_dir_symy0(x,y) = [undef, 0]
u_bound_dir_symx0(x,y) = [0, undef]
# traction at infinity
tx = 10


#define the exact solution
function exact_disp(x,y)
    cart2pol(x, y) = hypot(x, y), atan(y, x)
    r, th = cart2pol(x, y)
    ux = (1+nu)/Emod*tx*(1/(1+nu)*r*cos(th)+2*R^2/((1+nu)*r)*cos(th)+R^2/(2*r)*cos(3*th)-R^4/(2*r^3)*cos(3*th))
    uy = (1+nu)/Emod*tx*(-nu/(1+nu)*r*sin(th)-(1-nu)*R^2/((1+nu)*r)*sin(th)+
            R^2/(2*r)*sin(3*th)-R^4/(2*r^3)*sin(3*th))
    return [ux, uy]
end

function exact_stress(x,y)
    cart2pol(x, y) = hypot(x, y), atan(y, x)
    r, th = cart2pol(x, y)
    sigma_rr = tx/2*(1-R^2/r^2)+tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*th) 
    sigma_tt = tx/2*(1+R^2/r^2)-tx/2*(1+3*R^4/r^4)*cos(2*th)
    sigma_rt = -tx/2*(1+2*R^2/r^2-3*R^4/r^4)*sin(2*th)

    A = [cos(th)^2 sin(th)^2 2*sin(th)*cos(th); sin(th)^2 cos(th)^2 -2*sin(th)*cos(th);
        -sin(th)*cos(th) sin(th)*cos(th) cos(th)^2-sin(th)^2]
    stress = A\[sigma_rr;sigma_tt;sigma_rt]
    return stress
end
"""
Define the Neumann boundary condition derived from the exact stresses
"""
function u_bound_neu(x,y)
    tol_eq=1e-10
    ex_stress = exact_stress(x,y)
    if abs(x+L)<tol_eq 
        #on the boundary x=-L the outer normal is (-1, 0)
        nx = -1
        ny = 0
    elseif abs(y-L)<tol_eq
        #on the boundary y=L, the outer normal is (0, 1)
        nx = 0
        ny = 1
    else
        error("Point not on the Neumann boundary")
    end
    tx = nx*ex_stress[1]+ny*ex_stress[3]
    ty = nx*ex_stress[3]+ny*ex_stress[2]
    return tx, ty
end

bound_left = Boundary2D("Dirichlet", "Left", u_bound_dir_symy0)
bound_right = Boundary2D("Dirichlet", "Right", u_bound_dir_symx0)
bound_up = Boundary2D("Neumann", "Up", u_bound_neu)
bound_all = [bound_left, bound_right, bound_up]

#generate the NURBS mesh
IGAmesh = genMesh(nrb)

gauss_rule = [genGaussLegendre(degP[1]+1), genGaussLegendre(degP[2]+1)]
stiff = assemble_stiff_elast2D(IGAmesh, Cmat, gauss_rule)
#Assume no body forces  
rhs = zeros(2*IGAmesh.numBasis)
bcdof_all, elem_all = classifyBoundary2D(IGAmesh)
lhs, rhs = applyBCnurbsElast(IGAmesh, bound_all, stiff, rhs, bcdof_all, elem_all, gauss_rule)
println("Solving linear system with $(length(rhs)) DOFs")
@time sol0 = lhs\rhs
@time plotSolElast(IGAmesh, Cmat, sol0, "PlateHoleC0")
@time plotSolErrorElast(IGAmesh, Cmat, sol0, exact_disp, exact_stress, "PlateHoleC0Err")
@time relL2Err, relH1Err = compErrorNormElast(IGAmesh, Cmat, sol0, exact_disp, exact_stress, gauss_rule)
println("Relative L2-norm error is $relL2Err")
println("Relative energy-norm error is $relH1Err")

