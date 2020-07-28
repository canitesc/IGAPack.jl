
using SparseArrays
using FastGaussQuadrature

include("SplineModel.jl")
include("MeshUtils.jl")

struct Boundary1D
    type::String #can be "Dirichlet, Neumann or Robin"
    x_val::Float64 #boundary point in the physical space
    u_val::Float64 #boundary point in the parameter space
    α::Number #parameter value for Robin boundary condition (αu+u')
    op_val::Number #boundary operator value
end

function Boundary1D(type::String, x_val::Float64, u_val::Float64, op_val::Number)
    return Boundary1D(type, x_val, u_val, 0., op_val)
end

struct Problem1D
    rhs_f::Function
    a0::Function
    a1::Function
    boundary_cond::Array{Boundary1D}
    domain::NURBS
end

struct GaussQuad
    nodes::Array{Float64,1}
    weights::Array{Float64,1}
end

"""
Generates the Gauss-Legendre points of order n
"""
function genGaussLegendre(n)
    nodes, weights = gausslegendre(n)
    gauss_rule = GaussQuad(nodes, weights)
    return gauss_rule
end

"""
Aseembles the stiffness matrix K_ij = ∫ a0(x)ϕ_i'(x)ϕ_j'(x) dΩ
"""
function assemble_stiff(mesh::Mesh, a0::Function, gauss_rule)
    B, dB = bernsteinBasis(gauss_rule.nodes, degP)
    domainLength = 0
    stiff = spzeros(mesh.numBasis, mesh.numBasis)
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 2]
        Jac_ref_par = (uMax-uMin)/2

        #compute the (B-)spline basis functions and derivatives with Bezier extraction
        N_mat = B * mesh.C[iElem]'
        dN_mat = dB * mesh.C[iElem]'/Jac_ref_par

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1, curNodes]
        wgts = mesh.weights[curNodes]
        localStiff = zeros(numNodes, numNodes)
        for iGauss = 1:length(gauss_rule.nodes)
            #compute the rational basis
            RR = N_mat[iGauss,:].* wgts
            dR = dN_mat[iGauss,:].* wgts
            w_sum = sum(RR)
            dw_xi = sum(dR)
            dR = dR/w_sum - RR*dw_xi/w_sum^2

            #compute the derivatives w.r.t to the physical space
            dxdxi = dR' * cpts
            dR = dxdxi\dR
            Jac_par_phys = det(dxdxi)

            RR /= w_sum
            phys_pt = RR'*cpts
            localStiff += Jac_par_phys * Jac_ref_par * a0(phys_pt) * (dR*dR') * gauss_rule.weights[iGauss]
            domainLength += Jac_par_phys * Jac_ref_par * gauss_rule.weights[iGauss]

        end
        stiff[curNodes, curNodes] += localStiff

    end
    @show domainLength
    return stiff
end

"""
Aseembles the mass matrix M_ij = ∫ a1(x)ϕ_i(x)ϕ_j(x) dΩ
"""
function assemble_mass(mesh::Mesh, a1::Function, gauss_rule)
    B, dB = bernsteinBasis(gauss_rule.nodes, degP)
    domainLength = 0
    mass = spzeros(Complex{Float64}, mesh.numBasis, mesh.numBasis)
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 2]
        Jac_ref_par = (uMax-uMin)/2

        #compute the (B-)spline basis functions and derivatives with Bezier extraction
        N_mat = B * mesh.C[iElem]'
        dN_mat = dB * mesh.C[iElem]'/Jac_ref_par

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1, curNodes]
        wgts = mesh.weights[curNodes]
        localMass = zeros(numNodes, numNodes)
        for iGauss = 1:length(gauss_rule.nodes)
            #compute the rational basis
            RR = N_mat[iGauss,:].* wgts
            dR = dN_mat[iGauss,:].* wgts
            w_sum = sum(RR)
            dw_xi = sum(dR)
            dR = dR/w_sum - RR*dw_xi/w_sum^2

            #compute the Jacobian of the transformation from parameter to physical space
            dxdxi = dR' * cpts
            Jac_par_phys = det(dxdxi)

            RR /= w_sum
            phys_pt = RR'*cpts
            localMass += a1(phys_pt) * (RR*RR') * Jac_par_phys * Jac_ref_par * gauss_rule.weights[iGauss]
            domainLength += Jac_par_phys * Jac_ref_par * gauss_rule.weights[iGauss]

        end
        #@show localMass
        #readline(stdin)
        mass[curNodes, curNodes] += localMass
    end
    @show domainLength
    return mass
end

"""
Assembles the RHS vector corresponding to the body force
RHS[i]=∫_Ω ϕ_i(x)*f(x) dΩ
"""
function assemble_rhs(mesh::Mesh, f::Function, gauss_rule)
    B, dB = bernsteinBasis(gauss_rule.nodes, degP)
    domainLength = 0
    rhs = zeros(Complex{Float64}, mesh.numBasis)
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 2]
        Jac_ref_par = (uMax-uMin)/2
        #@show Jac_ref_par
        #compute the (B-)spline basis functions and derivatives with Bezier extraction
        N_mat = B * mesh.C[iElem]'
        dN_mat = dB * mesh.C[iElem]'/Jac_ref_par

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1, curNodes]
        wgts = mesh.weights[curNodes]
        localRhs = zeros(numNodes)
        for iGauss = 1:length(gauss_rule.nodes)
            #compute the rational basis
            RR = N_mat[iGauss,:].* wgts
            dR = dN_mat[iGauss,:].* wgts
            w_sum = sum(RR)
            dw_xi = sum(dR)
            dR = dR/w_sum - RR*dw_xi/w_sum^2

            #compute the derivatives w.r.t to the physical space
            dxdxi = dR' * cpts
            Jac_par_phys = det(dxdxi)
            RR /= w_sum

            phys_pt = RR'*cpts
            localRhs += Jac_par_phys * Jac_ref_par * f(phys_pt) * RR * gauss_rule.weights[iGauss]
            domainLength += Jac_par_phys * Jac_ref_par * gauss_rule.weights[iGauss]
        end
        rhs[curNodes] += localRhs
    end
    @show domainLength
    return rhs
end


function uniquetol(arr, tol=1e-10)
    sort!(arr)
    delete_index = Array{Int64,1}(undef, 0)
    for i=1:length(arr)-1
        if abs(arr[i+1]-arr[i])<tol
            push!(delete_index, i+1)
        end
    end
    deleteat!(arr, delete_index)
end

"""
Helper function to compute the setdifference within a tolerance
"""
function setdifftol(arrA, arrB, tol=1e-10)
    delete_index = Array{Int64,1}(undef, 0)
    for i=1:length(arrA)
        for j=1:length(arrB)
            if abs(arrA[i]-arrB[j])<tol
                push!(delete_index, i)
                break
            end
        end
    end
    deleteat!(arrA, delete_index)
end

"""
Apply the boundary conditions to a linear system for a 1D
nurbs geometry
"""
function applyBCnurbs(bound_cond, stiff, mass, rhs, nrb)
    #collect the dofs and values corresponding to the boundary
    lhs = stiff+mass
    bcdof = Array{Int64,1}(undef, 0)
    bcval = Array{Float64,1}(undef, 0)
    for i=1:length(bound_cond)

        if bound_cond[i].type=="Dirichlet"
            bcdof = vcat(bcdof, findall(nrb.coefs[1,:].==bound_cond[i].x_val))
            bcval = vcat(bcval, bound_cond[i].op_val)
        elseif bound_cond[i].type=="Neumann"
            bcdof_neu = findall(nrb.coefs[1,:].==bound_cond[i].x_val)
            rhs[bcdof_neu] .+= bound_cond[i].op_val
        elseif bound_cond[i].type=="Robin"
            #Assume that the basis functions have value 1 on the Robin boundary
            bcdof_rob = findall(nrb.coefs[1,:].==bound_cond[i].x_val)
            lhs[bcdof_rob, bcdof_rob] .+= bound_cond[i].α
            rhs[bcdof_rob] .+= bound_cond[i].op_val
        else
            error("Boundary condition type not implemented")
        end
    end
    rhs = rhs - lhs[:,bcdof]*bcval
    rhs[bcdof] = bcval
    lhs[bcdof, :] .= 0.
    lhs[:, bcdof] .= 0.
    lhs[bcdof, bcdof] = sparse(I, length(bcdof), length(bcdof))
    return lhs, rhs
end
