
using SparseArrays
using FastGaussQuadrature

include("SplineModel.jl")
include("MeshUtils.jl")

struct Boundary2D
    type::String #can be "Dirichlet", "Neumann" or "Robin"
    side::String #can be "Down" (v=0), "Right" (u=1), "Up" (v=1) or "Left" (u=0)
    α::Number #parameter value for Robin boundary condition (αu+u')
    op_val::Function #boundary operator value (as a function of (x,y))
    #To Do: Partial boundary conditions, Lagrange multipliers (?)
end

struct BoundaryIndex2D
    down::Array{Int64, 1}
    up::Array{Int64, 1}
    left::Array{Int64, 1}
    right::Array{Int64, 1}
end

struct Material2D
    Emod::Float64
    nu::Float64
    Cmat::Array{Float64,2}
end

"""
Constructor for Boundary2D which sets α=0
"""
function Boundary2D(type::String, side::String, op_val::Function)
    return Boundary2D(type, side, 0., op_val)
end

struct Problem2D
    rhs_f::Function
    a0::Function
    a1::Function
    boundary_cond::Array{Boundary2D}
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
Generates the 2D Bernstein polynomials at points (pts_u, pts_v)
"""
function bernsteinBasis2D(pts_u::Array{Float64,1}, pts_v::Array{Float64,1},
                            degP::Array{Int64,1})
    B_u, dB_u = bernsteinBasis(pts_u, degP[1])
    B_v, dB_v = bernsteinBasis(pts_v, degP[2])

    numPtsU = length(pts_u)
    numPtsV = length(pts_v)
    basisCounter = 0
    Buv = zeros(numPtsU, numPtsV, (degP[1]+1)*(degP[2]+1))
    dBdu = zeros(numPtsU, numPtsV, (degP[1]+1)*(degP[2]+1))
    dBdv = zeros(numPtsU, numPtsV, (degP[1]+1)*(degP[2]+1))
    for j=1:degP[2]+1
        for i=1:degP[1]+1
            basisCounter += 1
            Buv[:,:,basisCounter] = B_u[:,i]*B_v[:,j]'
            dBdu[:,:,basisCounter] = dB_u[:,i]*B_v[:,j]'
            dBdv[:,:,basisCounter] = B_u[:,i]*dB_v[:,j]'
        end
    end
    return Buv, dBdu, dBdv
end

"""
Aseembles the stiffness matrix K_ij = ∫ a0(x)∇ϕ_i(x,y)∇ϕ_j(x) dΩ
"""
function assemble_stiff2D(mesh::Mesh, a0::Function, gauss_rule::Array{GaussQuad,1})
    Buv, dBdu, dBdv = bernsteinBasis2D(gauss_rule[1].nodes, gauss_rule[2].nodes, mesh.degP)
    numGaussU = length(gauss_rule[1].nodes)
    numGaussV = length(gauss_rule[2].nodes)
    #allocate memory for the triplet arrays
    indexCounter = 0
    for iElem=1:mesh.numElem
        numNodes = length(mesh.elemNode[iElem])
        indexCounter += numNodes^2
    end
    II = zeros(Int64, indexCounter)
    JJ = zeros(Int64, indexCounter)
    S = zeros(indexCounter)

    indexCounter = 0
    domainArea = 0
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        Jac_ref_par = (uMax-uMin)*(vMax-vMin)/4

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localStiff = zeros(numNodes, numNodes)

        for jGauss = 1:numGaussV
            for iGauss = 1:numGaussU
                #compute the (B-)spline basis functions and derivatives with Bezier extraction
                N_mat = mesh.C[iElem] * Buv[iGauss, jGauss, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, jGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, jGauss, :] * 2/(vMax-vMin)

                #compute the rational basis
                RR = N_mat.* wgts
                dRdu = dN_du.* wgts
                dRdv = dN_dv.* wgts
                w_sum = sum(RR)
                dw_xi = sum(dRdu)
                dw_eta = sum(dRdv)

                dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
                dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2
                #compute the derivatives w.r.t to the physical space
                dR = [dRdu'; dRdv']

                dxdxi = dR * cpts'
                dR = dxdxi\dR
                Jac_par_phys = det(dxdxi)
                RR /= w_sum
                phys_pt = cpts*RR

                localArea = Jac_par_phys * Jac_ref_par *
                    gauss_rule[1].weights[iGauss] * gauss_rule[2].weights[jGauss]
                localStiff += localArea * a0(phys_pt[1], phys_pt[2]) * (dR'*dR)
                domainArea += localArea
            end
        end
        II[indexCounter+1:indexCounter+numNodes^2] = repeat(curNodes, numNodes)
        JJ[indexCounter+1:indexCounter+numNodes^2] = reshape(repeat(curNodes', numNodes, 1), numNodes^2)
        S[indexCounter+1:indexCounter+numNodes^2] = reshape(localStiff, numNodes^2)
        indexCounter += numNodes^2

    end
    @show domainArea
    stiff = sparse(II, JJ, S, mesh.numBasis, mesh.numBasis)
    return stiff
end

"""
Aseembles the stiffness matrix for 2D Elasticity Kₑ = ∫ B^T*C*B dΩ
"""
function assemble_stiff_elast2D(mesh::Mesh, Cmat::Array{Float64,2}, gauss_rule::Array{GaussQuad,1})
    Buv, dBdu, dBdv = bernsteinBasis2D(gauss_rule[1].nodes, gauss_rule[2].nodes, mesh.degP)
    numGaussU = length(gauss_rule[1].nodes)
    numGaussV = length(gauss_rule[2].nodes)
    #allocate memory for the triplet arrays
    indexCounter = 0
    for iElem=1:mesh.numElem
        numNodes = length(mesh.elemNode[iElem])
        indexCounter += 4*numNodes^2
    end
    II = zeros(Int64, indexCounter)
    JJ = zeros(Int64, indexCounter)
    S = zeros(indexCounter)

    indexCounter = 0
    domainArea = 0
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        Jac_ref_par = (uMax-uMin)*(vMax-vMin)/4

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        curNodesXY = reshape(hcat(2*curNodes.-1, 2*curNodes)', 2*numNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localStiff = zeros(2*numNodes, 2*numNodes)

        for jGauss = 1:numGaussV
            for iGauss = 1:numGaussU
                #compute the (B-)spline basis functions and derivatives with Bezier extraction
                N_mat = mesh.C[iElem] * Buv[iGauss, jGauss, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, jGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, jGauss, :] * 2/(vMax-vMin)

                #compute the rational basis
                RR = N_mat.* wgts
                dRdu = dN_du.* wgts
                dRdv = dN_dv.* wgts
                w_sum = sum(RR)
                dw_xi = sum(dRdu)
                dw_eta = sum(dRdv)
                dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
                dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2
                
                #compute the derivatives w.r.t to the physical space
                dR = [dRdu'; dRdv']
                dxdxi = dR * cpts'
                dR = dxdxi\dR
                Jac_par_phys = det(dxdxi)
                RR /= w_sum
                phys_pt = cpts*RR

                B = zeros(2*numNodes,3);
                B[1:2:2*numNodes-1,1] = dR[1,:]
                B[2:2:2*numNodes,2] = dR[2,:]
                B[1:2:2*numNodes-1,3] = dR[2,:]
                B[2:2:2*numNodes,3] = dR[1,:]

                localArea = Jac_par_phys * Jac_ref_par *
                    gauss_rule[1].weights[iGauss] * gauss_rule[2].weights[jGauss]
                localStiff += localArea *  (B*Cmat*B')
                domainArea += localArea
            end
        end
        II[indexCounter+1:indexCounter+4*numNodes^2] = repeat(curNodesXY, 2*numNodes)
        JJ[indexCounter+1:indexCounter+4*numNodes^2] = reshape(repeat(curNodesXY',
                                                    2*numNodes, 1), 4*numNodes^2)
        S[indexCounter+1:indexCounter+4*numNodes^2] = reshape(localStiff, 4*numNodes^2)
        indexCounter += 4*numNodes^2

    end
    @show domainArea
    stiff = sparse(II, JJ, S, 2*mesh.numBasis, 2*mesh.numBasis)
    return stiff
end

"""
Aseembles the mass matrix M_ij = ∫ a1(x)ϕ_i(x,y)ϕ_j(x) dΩ
"""
function assemble_mass2D(mesh::Mesh, a1::Function, gauss_rule::Array{GaussQuad,1})
    Buv, dBdu, dBdv = bernsteinBasis2D(gauss_rule[1].nodes, gauss_rule[2].nodes, mesh.degP)
    numGaussU = length(gauss_rule[1].nodes)
    numGaussV = length(gauss_rule[2].nodes)
    #allocate memory for the triplet arrays
    indexCounter = 0
    for iElem=1:mesh.numElem
        numNodes = length(mesh.elemNode[iElem])
        indexCounter += numNodes^2
    end
    II = zeros(Int64, indexCounter)
    JJ = zeros(Int64, indexCounter)
    S = zeros(indexCounter)

    indexCounter = 0
    domainArea = 0
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        Jac_ref_par = (uMax-uMin)*(vMax-vMin)/4

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localMass = zeros(numNodes, numNodes)

        for jGauss = 1:numGaussV
            for iGauss = 1:numGaussU
                #compute the (B-)spline basis functions and derivatives with Bezier extraction
                N_mat = mesh.C[iElem] * Buv[iGauss, jGauss, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, jGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, jGauss, :] * 2/(vMax-vMin)

                #compute the rational basis
                RR = N_mat.* wgts
                dRdu = dN_du.* wgts
                dRdv = dN_dv.* wgts
                w_sum = sum(RR)
                dw_xi = sum(dRdu)
                dw_eta = sum(dRdv)

                dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
                dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2
                #compute the derivatives w.r.t to the physical space
                dR = [dRdu'; dRdv']

                dxdxi = dR * cpts'
                dR = dxdxi\dR
                Jac_par_phys = det(dxdxi)
                RR /= w_sum
                phys_pt = cpts*RR

                localArea = Jac_par_phys * Jac_ref_par *
                    gauss_rule[1].weights[iGauss] * gauss_rule[2].weights[jGauss]
                localMass += localArea * a1(phys_pt[1], phys_pt[2]) * (RR*RR')
                domainArea += localArea
            end
        end
        II[indexCounter+1:indexCounter+numNodes^2] = repeat(curNodes, numNodes)
        JJ[indexCounter+1:indexCounter+numNodes^2] = reshape(repeat(curNodes', numNodes, 1), numNodes^2)
        S[indexCounter+1:indexCounter+numNodes^2] = reshape(localMass, numNodes^2)
        indexCounter += numNodes^2

    end
    @show domainArea
    stiff = sparse(II, JJ, S, mesh.numBasis, mesh.numBasis)
    return stiff
end

"""
Assembles the RHS vector corresponding to the body force
RHS[i]=∫_Ω ϕ_i(x,y)*f(x,y) dΩ
"""
function assemble_rhs2D(mesh::Mesh, f::Function,  gauss_rule::Array{GaussQuad,1})
    Buv, dBdu, dBdv = bernsteinBasis2D(gauss_rule[1].nodes, gauss_rule[2].nodes, mesh.degP)
    numGaussU = length(gauss_rule[1].nodes)
    numGaussV = length(gauss_rule[2].nodes)
    domainArea = 0
    rhs = zeros(Complex{Float64}, mesh.numBasis)
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        Jac_ref_par = (uMax-uMin)*(vMax-vMin)/4

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localRhs = zeros(numNodes)
        for jGauss = 1:numGaussV
            for iGauss = 1:numGaussU
                #compute the (B-)spline basis functions and derivatives with Bezier extraction
                N_mat = mesh.C[iElem] * Buv[iGauss, jGauss, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, jGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, jGauss, :] * 2/(vMax-vMin)

                #compute the rational basis
                RR = N_mat.* wgts
                dRdu = dN_du.* wgts
                dRdv = dN_dv.* wgts
                w_sum = sum(RR)
                dw_xi = sum(dRdu)
                dw_eta = sum(dRdv)

                dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
                dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2
                #compute the derivatives w.r.t to the physical space
                dR = [dRdu'; dRdv']

                dxdxi = dR * cpts'
                dR = dxdxi\dR
                Jac_par_phys = det(dxdxi)

                RR /= w_sum
                phys_pt = cpts*RR
                localArea = Jac_par_phys * Jac_ref_par *
                    gauss_rule[1].weights[iGauss] * gauss_rule[2].weights[jGauss]
                localRhs += localArea * f(phys_pt[1], phys_pt[2]) * RR
                domainArea += localArea
            end
        end
        rhs[curNodes] += localRhs
    end
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
Returns the boundary (down, right, up, left) indices for a matrix of dimension
(p+1)*(q+1)
"""
function getBoundaryIndices(degP)
    indexDown = 1:degP[1]+1
    indexRight = degP[1]+1:degP[1]+1:prod(degP.+1)
    indexUp = 1+(degP[1]+1)*degP[2]:prod(degP.+1)
    indexLeft = 1:degP[1]+1:1+(degP[1]+1)*degP[2]
    index_all = BoundaryIndex2D(indexDown, indexUp, indexLeft, indexRight)
    return index_all
end

"""
Classifies the boundary nodes and boundary elements according to the side in the
parameter space (i.e. bcdof_down and elem_down for v=0, bcdof_up and elem_up
for v=1, bcdof_left and elem_left for u=0, bcdof_right and elem_right for u=1)
"""
function classifyBoundary2D(mesh::Mesh)
    bcdof_down = Array{Int64,1}(undef, 0)
    elem_down = Array{Int64, 1}(undef, 0)
    bcdof_up = Array{Int64,1}(undef, 0)
    elem_up = Array{Int64,1}(undef, 0)
    bcdof_left = Array{Int64,1}(undef, 0)
    elem_left = Array{Int64,1}(undef, 0)
    bcdof_right = Array{Int64,1}(undef, 0)
    elem_right = Array{Int64,1}(undef, 0)
    index_all = getBoundaryIndices(mesh.degP)
    tolEq = 1e-10

    for iElem = 1:mesh.numElem
        if abs(mesh.elemVertex[iElem,1])<tolEq
            #u_min = 0
            bcdof_left = vcat(bcdof_left, mesh.elemNode[iElem][index_all.left])
            elem_left = vcat(elem_left, iElem)
        end
        if abs(mesh.elemVertex[iElem,2])<tolEq
            #v_min = 0
            bcdof_down = vcat(bcdof_down, mesh.elemNode[iElem][index_all.down])
            elem_down = vcat(elem_down, iElem)
        end
        if abs(mesh.elemVertex[iElem,3]-1)<tolEq
            #u_max = 1
            bcdof_right = vcat(bcdof_right, mesh.elemNode[iElem][index_all.right])
            elem_right = vcat(elem_right, iElem)
        end
        if abs(mesh.elemVertex[iElem,4]-1)<tolEq
            #v_max = 1
            bcdof_up = vcat(bcdof_up, mesh.elemNode[iElem][index_all.up])
            elem_up = vcat(elem_up, iElem)
        end
    end
    bcdof_down = unique(bcdof_down)
    bcdof_up = unique(bcdof_up)
    bcdof_left = unique(bcdof_left)
    bcdof_right = unique(bcdof_right)
    bcdof_all = BoundaryIndex2D(bcdof_down, bcdof_up, bcdof_left, bcdof_right)
    elem_all = BoundaryIndex2D(elem_down, elem_up, elem_left, elem_right)
    return bcdof_all, elem_all
end

"""
Apply the Neumann boundary conditions to the RHS for a 2D NURBS geometry, scalar
solution field
"""
function applyNeumannScalar2D(mesh::Mesh, rhs, elem_list, direction::String,
                                quad_rule_1d::GaussQuad, op_val::Function)
    #Select the Gauss/edge points for evaluating the  1D the basis functions
    #select the scatter vector corresponding to the edge DOFs
    index_all = getBoundaryIndices(mesh.degP)
    numNodesEdge = length(quad_rule_1d.nodes)
    if direction == "right"
        pts_u = [1.]
        pts_v = quad_rule_1d.nodes
        sctr = index_all.right
    elseif direction == "left"
        pts_u = [-1.]
        pts_v = quad_rule_1d.nodes
        sctr = index_all.left
    elseif direction == "down"
        pts_u = quad_rule_1d.nodes
        pts_v = [-1.]
        sctr = index_all.down
    elseif direction == "up"
        pts_u = quad_rule_1d.nodes
        pts_v = [1.]
        sctr = index_all.up
    else
        error("Wrong direction given")
    end

    #Form the 2D tensor product of the basis functions
    Buv, dBdu, dBdv = bernsteinBasis2D(pts_u, pts_v, mesh.degP)

    #Evaluate the Neumann integral on each element
    domainLength = 0.
    for iElem in elem_list
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        if direction == "right" || direction == "left"
            Jac_ref_par = (vMax-vMin)/2
        else
            Jac_ref_par = (uMax-uMin)/2
        end

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localRhs = zeros(numNodesEdge)
        for iGauss = 1:length(quad_rule_1d.nodes)
            #compute the (B-)spline basis functions and derivatives with Bezier extraction
            if direction == "right" || direction == "left"
                N_mat = mesh.C[iElem] * Buv[1, iGauss, :]
                dN_du = mesh.C[iElem] * dBdu[1, iGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[1, iGauss, :] * 2/(vMax-vMin)
            else
                N_mat = mesh.C[iElem] * Buv[iGauss, 1, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, 1, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, 1, :] * 2/(vMax-vMin)
            end

            #compute the rational basis
            RR = N_mat.* wgts
            dRdu = dN_du.* wgts
            dRdv = dN_dv.* wgts
            w_sum = sum(RR)
            dw_xi = sum(dRdu)
            dw_eta = sum(dRdv)
            dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
            dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2

            #compute the derivatives w.r.t to the physical space
            dR = [dRdu'; dRdv']
            dxdxi = dR * cpts'

            #Jacobian of face mapping
            if direction=="right" || direction=="left"
                eJac = dxdxi[2,1]^2 + dxdxi[2,2]^2
            else
                eJac = dxdxi[1,1]^2 + dxdxi[1,2]^2
            end
            Jac_par_phys = sqrt(eJac)
            dR = dxdxi\dR
            RR /= w_sum
            phys_pt = cpts*RR
            gFunc = op_val(phys_pt[1], phys_pt[2])
            localLength = Jac_par_phys * Jac_ref_par * quad_rule_1d.weights[iGauss]
            localRhs += RR[sctr].*gFunc.*localLength
            domainLength += localLength
        end
        rhs[curNodes[sctr]] += localRhs
    end
    @show domainLength
    return rhs
end

"""
Apply the Neumann boundary conditions to the RHS for a 2D NURBS geometry, elasticity
solution field
"""
function applyNeumannElast2D(mesh::Mesh, rhs, elem_list, direction::String,
                                quad_rule_1d::GaussQuad, op_val::Function)
    #Select the Gauss/edge points for evaluating the  1D the basis functions
    #select the scatter vector corresponding to the edge DOFs
    index_all = getBoundaryIndices(mesh.degP)
    numNodesEdge = length(quad_rule_1d.nodes)
    if direction == "right"
        pts_u = [1.]
        pts_v = quad_rule_1d.nodes
        sctr = index_all.right
    elseif direction == "left"
        pts_u = [-1.]
        pts_v = quad_rule_1d.nodes
        sctr = index_all.left
    elseif direction == "down"
        pts_u = quad_rule_1d.nodes
        pts_v = [-1.]
        sctr = index_all.down
    elseif direction == "up"
        pts_u = quad_rule_1d.nodes
        pts_v = [1.]
        sctr = index_all.up
    else
        error("Wrong direction given")
    end

    #Form the 2D tensor product of the basis functions
    Buv, dBdu, dBdv = bernsteinBasis2D(pts_u, pts_v, mesh.degP)

    #Evaluate the Neumann integral on each element
    domainLength = 0.
    for iElem in elem_list
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        if direction == "right" || direction == "left"
            Jac_ref_par = (vMax-vMin)/2
        else
            Jac_ref_par = (uMax-uMin)/2
        end

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        curNodesXY = reshape(hcat(2*curNodes[sctr].-1, 2*curNodes[sctr])', 2*numNodesEdge)

        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localRhs = zeros(2*numNodesEdge)
        for iGauss = 1:length(quad_rule_1d.nodes)
            #compute the (B-)spline basis functions and derivatives with Bezier extraction
            if direction == "right" || direction == "left"
                N_mat = mesh.C[iElem] * Buv[1, iGauss, :]
                dN_du = mesh.C[iElem] * dBdu[1, iGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[1, iGauss, :] * 2/(vMax-vMin)
            else
                N_mat = mesh.C[iElem] * Buv[iGauss, 1, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, 1, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, 1, :] * 2/(vMax-vMin)
            end

            #compute the rational basis
            RR = N_mat.* wgts
            dRdu = dN_du.* wgts
            dRdv = dN_dv.* wgts
            w_sum = sum(RR)
            dw_xi = sum(dRdu)
            dw_eta = sum(dRdv)
            dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
            dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2

            #compute the derivatives w.r.t to the physical space
            dR = [dRdu'; dRdv']
            dxdxi = dR * cpts'

            #Jacobian of face mapping
            if direction=="right" || direction=="left"
                eJac = dxdxi[2,1]^2 + dxdxi[2,2]^2
            else
                eJac = dxdxi[1,1]^2 + dxdxi[1,2]^2
            end
            Jac_par_phys = sqrt(eJac)

            dR = dxdxi\dR
            RR /= w_sum
            phys_pt = cpts*RR

            gFunc = op_val(phys_pt[1], phys_pt[2])
            localLength = Jac_par_phys * Jac_ref_par * quad_rule_1d.weights[iGauss]
            localRhs[1:2:end-1] += RR[sctr].*gFunc[1].*localLength
            localRhs[2:2:end] += RR[sctr].*gFunc[2].*localLength
            domainLength += localLength

        end
        rhs[curNodesXY] += localRhs
    end
    @show domainLength
    return rhs
end

"""
Apply the Robin boundary conditions to the RHS for a 2D NURBS geometry, scalar
solution field
"""
function applyRobinScalar2D(mesh::Mesh, lhs, rhs, elem_list, direction::String,
                                quad_rule_1d::GaussQuad, α::Number, op_val::Function)
    #Select the Gauss/edge points for evaluating the  1D the basis functions
    #select the scatter vector corresponding to the edge DOFs
    index_all = getBoundaryIndices(mesh.degP)
    numNodesEdge = length(quad_rule_1d.nodes)
    if direction == "right"
        pts_u = [1.]
        pts_v = quad_rule_1d.nodes
        sctr = index_all.right
    elseif direction == "left"
        pts_u = [-1.]
        pts_v = quad_rule_1d.nodes
        sctr = index_all.left
    elseif direction == "down"
        pts_u = quad_rule_1d.nodes
        pts_v = [-1.]
        sctr = index_all.down
    elseif direction == "up"
        pts_u = quad_rule_1d.nodes
        pts_v = [1.]
        sctr = index_all.up
    else
        error("Wrong direction given")
    end

    #allocate memory for the triplet arrays
    indexCounter = length(elem_list)*length(sctr)^2
    II = zeros(Int64, indexCounter)
    JJ = zeros(Int64, indexCounter)
    S = zeros(Complex{Float64}, indexCounter)
    indexCounter = 0

    #Form the 2D tensor product of the basis functions
    Buv, dBdu, dBdv = bernsteinBasis2D(pts_u, pts_v, mesh.degP)

    #Evaluate the Neumann integral on each element
    domainLength = 0.
    for iElem in elem_list
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        if direction == "right" || direction == "left"
            Jac_ref_par = (vMax-vMin)/2
        else
            Jac_ref_par = (uMax-uMin)/2
        end

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(sctr)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localBndMass = zeros(numNodesEdge,numNodesEdge)
        for iGauss = 1:length(quad_rule_1d.nodes)
            #compute the (B-)spline basis functions and derivatives with Bezier extraction
            if direction == "right" || direction == "left"
                N_mat = mesh.C[iElem] * Buv[1, iGauss, :]
                dN_du = mesh.C[iElem] * dBdu[1, iGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[1, iGauss, :] * 2/(vMax-vMin)
            else
                N_mat = mesh.C[iElem] * Buv[iGauss, 1, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, 1, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, 1, :] * 2/(vMax-vMin)
            end

            #compute the rational basis
            RR = N_mat.* wgts
            dRdu = dN_du.* wgts
            dRdv = dN_dv.* wgts
            w_sum = sum(RR)
            dw_xi = sum(dRdu)
            dw_eta = sum(dRdv)
            dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
            dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2

            #compute the derivatives w.r.t to the physical space
            dR = [dRdu'; dRdv']
            dxdxi = dR * cpts'

            #Jacobian of face mapping
            if direction=="right" || direction=="left"
                eJac = dxdxi[1,2]^2 + dxdxi[2,2]^2
            else
                eJac = dxdxi[1,1]^2 + dxdxi[2,1]^2
            end
            Jac_par_phys = sqrt(eJac)
            dR = dxdxi\dR
            RR /= w_sum
            phys_pt = cpts*RR

            localBndMass += (RR[sctr]*RR[sctr]').*α.*Jac_ref_par.*Jac_par_phys.*quad_rule_1d.weights[iGauss]
            domainLength += Jac_par_phys * Jac_ref_par * quad_rule_1d.weights[iGauss]
        end
        II[indexCounter+1:indexCounter+numNodes^2] = repeat(curNodes[sctr], numNodes)
        JJ[indexCounter+1:indexCounter+numNodes^2] = reshape(repeat(curNodes[sctr]', numNodes, 1), numNodes^2)
        S[indexCounter+1:indexCounter+numNodes^2] = reshape(localBndMass, numNodes^2)
        indexCounter += numNodes^2
    end
    @show domainLength
    bnd_mass = sparse(II, JJ, S, mesh.numBasis, mesh.numBasis)
    lhs += bnd_mass
    return lhs, rhs
end

"""
Apply the boundary conditions to a linear system for a 2D
NURBS geometry and elasticity problems
"""
function applyBCnurbsElast(mesh::Mesh, bound_cond::Array{Boundary2D, 1}, lhs, rhs,
                            bcdof_all, elem_all, quad_rule)
    #collect the dofs and values corresponding to the boundary
    bcdof = Array{Int64,1}(undef, 0)
    bcval = Array{Float64,1}(undef, 0)
    #assume homogeneous boundary conditions
    evalPt = [0.,0.]
    for i=1:length(bound_cond)
        if bound_cond[i].type=="Dirichlet"
            if bound_cond[i].side=="Down"
                #check x-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[1]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.down.-1)
                end
                #check y-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[2]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.down)
                end
            end
            if bound_cond[i].side=="Up"
                #check x-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[1]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.up.-1)
                end
                #check y-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[2]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.up)
                end
            end
            if bound_cond[i].side=="Left"
                #check x-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[1]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.left.-1)
                end
                #check y-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[2]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.left)
                end
            end
            if bound_cond[i].side=="Right"
                #check x-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[1]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.right.-1)
                end
                #check y-direction
                if bound_cond[i].op_val(evalPt[1], evalPt[2])[2]!=undef
                    bcdof = vcat(bcdof, 2*bcdof_all.right)
                end
            end
        elseif bound_cond[i].type=="Neumann"
            if bound_cond[i].side=="Down"
                rhs = applyNeumannElast2D(mesh, rhs, elem_all.down, "down",
                                        quad_rule[1], bound_cond[i].op_val)
            end
            if bound_cond[i].side=="Up"
                rhs = applyNeumannElast2D(mesh, rhs, elem_all.up, "up",
                                        quad_rule[1], bound_cond[i].op_val)
            end
            if bound_cond[i].side=="Left"
                rhs = applyNeumannElast2D(mesh, rhs, elem_all.left, "left",
                                        quad_rule[2], bound_cond[i].op_val)
            end
            if bound_cond[i].side=="Right"
                rhs = applyNeumannElast2D(mesh, rhs, elem_all.right, "right",
                                        quad_rule[2], bound_cond[i].op_val)
            end
        else
            error("Boundary condition type not implemented")
        end
    end
    unique!(bcdof)
    bcval = zero(bcdof)
    rhs = rhs - lhs[:,bcdof]*bcval
    rhs[bcdof] = bcval
    lhs[bcdof, :] .= 0.
    lhs[:, bcdof] .= 0.
    lhs[bcdof, bcdof] = sparse(I, length(bcdof), length(bcdof))
    return lhs, rhs
end

function applyBCnurbs(mesh::Mesh, bound_cond::Array{Boundary2D, 1}, lhs, rhs,
                        bcdof_all, elem_all, quad_rule)
    #collect the dofs and values corresponding to the boundary
    bcdof = Array{Int64,1}(undef, 0)
    bcval = Array{Float64,1}(undef, 0)
    for i=1:length(bound_cond)
        if bound_cond[i].type=="Dirichlet"
            if bound_cond[i].side=="Down"
                bcdof = vcat(bcdof, bcdof_all.down)
            end
            if bound_cond[i].side=="Up"
                bcdof = vcat(bcdof, bcdof_all.up)
            end
            if bound_cond[i].side=="Left"
                bcdof = vcat(bcdof, bcdof_all.left)
            end
            if bound_cond[i].side=="Right"
                bcdof = vcat(bcdof, bcdof_all.right)
            end
        elseif bound_cond[i].type=="Neumann"
            if bound_cond[i].side=="Down"
                rhs = applyNeumannScalar2D(mesh, rhs, elem_all.down, "down",
                                        quad_rule[1], bound_cond[i].op_val)
            end
            if bound_cond[i].side=="Up"
                rhs = applyNeumannScalar2D(mesh, rhs, elem_all.up, "up",
                                        quad_rule[1], bound_cond[i].op_val)
            end
            if bound_cond[i].side=="Left"
                rhs = applyNeumannScalar2D(mesh, rhs, elem_all.left, "left",
                                        quad_rule[2], bound_cond[i].op_val)
            end
            if bound_cond[i].side=="Right"
                rhs = applyNeumannScalar2D(mesh, rhs, elem_all.right, "right",
                                        quad_rule[2], bound_cond[i].op_val)
            end
        elseif bound_cond[i].type=="Robin"
            α = bound_cond[i].α
            op_val = bound_cond[i].op_val
            if bound_cond[i].side=="Down"
                lhs, rhs = applyRobinScalar2D(mesh, lhs, rhs, elem_all.down, "down",
                                        quad_rule[1], α, op_val)
            end
            if bound_cond[i].side=="Up"
                lhs, rhs = applyRobinScalar2D(mesh, lhs, rhs, elem_all.up, "up",
                                        quad_rule[1], α, op_val)
            end
            if bound_cond[i].side=="Left"
                lhs, rhs = applyRobinScalar2D(mesh, lhs, rhs, elem_all.left, "left",
                                        quad_rule[2], α, op_val)
            end
            if bound_cond[i].side=="Right"
                lhs, rhs = applyRobinScalar2D(mesh, lhs, rhs, elem_all.right, "right",
                                        quad_rule[2], α, op_val)
            end
        else
            error("Boundary condition type not implemented")
        end
    end
    unique!(bcdof)
    bcval = zero(bcdof)
    rhs = rhs - lhs[:,bcdof]*bcval
    rhs[bcdof] = bcval
    lhs[bcdof, :] .= 0.
    lhs[:, bcdof] .= 0.
    lhs[bcdof, bcdof] = sparse(I, length(bcdof), length(bcdof))
    return lhs, rhs
end

function compErrorNorm(mesh::Mesh, sol0, exactSol::Function, derivExactSol::Function,
                        a0::Function, gauss_rule::Array{GaussQuad,1})
    B_u, dB_u = bernsteinBasis(gauss_rule[1].nodes, mesh.degP[1])
    B_v, dB_v = bernsteinBasis(gauss_rule[2].nodes, mesh.degP[2])
    domainArea = 0

    numGaussU = length(gauss_rule[1].nodes)
    numGaussV = length(gauss_rule[2].nodes)
    basisCounter = 0
    Buv = zeros(numGaussU, numGaussV, (degP[1]+1)*(degP[2]+1))
    dBdu = zeros(numGaussU, numGaussV, (degP[1]+1)*(degP[2]+1))
    dBdv = zeros(numGaussU, numGaussV, (degP[1]+1)*(degP[2]+1))

    for j=1:degP[2]+1
        for i=1:degP[1]+1
            basisCounter += 1
            Buv[:,:,basisCounter] = B_u[:,i]*B_v[:,j]'
            dBdu[:,:,basisCounter] = dB_u[:,i]*B_v[:,j]'
            dBdv[:,:,basisCounter] = B_u[:,i]*dB_v[:,j]'
        end
    end
    l2NormErr = 0.
    l2NormSol = 0.
    h1NormErr = 0.
    h1NormSol = 0.
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        Jac_ref_par = (uMax-uMin)*(vMax-vMin)/4

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localRhs = zeros(numNodes)
        for jGauss = 1:numGaussV
            for iGauss = 1:numGaussU
                #compute the (B-)spline basis functions and derivatives with Bezier extraction
                N_mat = mesh.C[iElem] * Buv[iGauss, jGauss, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, jGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, jGauss, :] * 2/(vMax-vMin)

                #compute the rational basis
                RR = N_mat.* wgts
                dRdu = dN_du.* wgts
                dRdv = dN_dv.* wgts
                w_sum = sum(RR)
                dw_xi = sum(dRdu)
                dw_eta = sum(dRdv)

                dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
                dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2
                #compute the derivatives w.r.t. the physical space
                dR = [dRdu'; dRdv']

                dxdxi = dR * cpts'
                dR = dxdxi\dR
                Jac_par_phys = det(dxdxi)

                RR /= w_sum
                phys_pt = cpts*RR
                solVal = RR'*sol0[curNodes]
                dXsolVal = dR*sol0[curNodes]
                exSolVal = real(exactSol(phys_pt[1], phys_pt[2]))
                dXexSolVal = real(derivExactSol(phys_pt[1], phys_pt[2]))
                localArea = Jac_par_phys * Jac_ref_par *
                    gauss_rule[1].weights[iGauss] * gauss_rule[2].weights[jGauss]
                domainArea += localArea
                l2NormErr += localArea * (exSolVal-solVal)^2
                l2NormSol += localArea*exSolVal^2
                funValue = a0(phys_pt[1], phys_pt[2])
                h1NormErr += localArea * funValue *sum((dXexSolVal.-dXsolVal).^2)
                h1NormSol += localArea * funValue *sum((dXexSolVal).^2)
            end
        end
    end
    relL2Err = sqrt(l2NormErr/l2NormSol)
    relH1Err = sqrt(h1NormErr/h1NormSol)
    return real(relL2Err), real(relH1Err)
end

function compErrorNormElast(mesh::Mesh, Cmat, sol0, exactSolDisp::Function, exactSolStress::Function,
                         gauss_rule::Array{GaussQuad,1})
    invC = inv(Cmat)
    B_u, dB_u = bernsteinBasis(gauss_rule[1].nodes, mesh.degP[1])
    B_v, dB_v = bernsteinBasis(gauss_rule[2].nodes, mesh.degP[2])
    domainArea = 0

    numGaussU = length(gauss_rule[1].nodes)
    numGaussV = length(gauss_rule[2].nodes)
    basisCounter = 0
    Buv = zeros(numGaussU, numGaussV, (degP[1]+1)*(degP[2]+1))
    dBdu = zeros(numGaussU, numGaussV, (degP[1]+1)*(degP[2]+1))
    dBdv = zeros(numGaussU, numGaussV, (degP[1]+1)*(degP[2]+1))

    for j=1:degP[2]+1
        for i=1:degP[1]+1
            basisCounter += 1
            Buv[:,:,basisCounter] = B_u[:,i]*B_v[:,j]'
            dBdu[:,:,basisCounter] = dB_u[:,i]*B_v[:,j]'
            dBdv[:,:,basisCounter] = B_u[:,i]*dB_v[:,j]'
        end
    end
    l2NormErr = 0.
    l2NormSol = 0.
    h1NormErr = 0.
    h1NormSol = 0.
    for iElem = 1:mesh.numElem
        uMin = mesh.elemVertex[iElem, 1]
        uMax = mesh.elemVertex[iElem, 3]
        vMin = mesh.elemVertex[iElem, 2]
        vMax = mesh.elemVertex[iElem, 4]
        Jac_ref_par = (uMax-uMin)*(vMax-vMin)/4

        #compute the rational spline basis
        curNodes = mesh.elemNode[iElem]
        numNodes = length(curNodes)
        curNodesXY = reshape(hcat(2*curNodes.-1, 2*curNodes)', 2*numNodes)
        cpts = mesh.controlPoints[1:2, curNodes]
        wgts = mesh.weights[curNodes]
        localRhs = zeros(numNodes)
        for jGauss = 1:numGaussV
            for iGauss = 1:numGaussU
                #compute the (B-)spline basis functions and derivatives with Bezier extraction
                N_mat = mesh.C[iElem] * Buv[iGauss, jGauss, :]
                dN_du = mesh.C[iElem] * dBdu[iGauss, jGauss, :] * 2/(uMax-uMin)
                dN_dv = mesh.C[iElem] * dBdv[iGauss, jGauss, :] * 2/(vMax-vMin)

                #compute the rational basis
                RR = N_mat.* wgts
                dRdu = dN_du.* wgts
                dRdv = dN_dv.* wgts
                w_sum = sum(RR)
                dw_xi = sum(dRdu)
                dw_eta = sum(dRdv)

                dRdu = dRdu/w_sum - RR*dw_xi/w_sum^2
                dRdv = dRdv/w_sum - RR*dw_eta/w_sum^2
                #compute the derivatives w.r.t to the physical space
                dR = [dRdu'; dRdv']

                dxdxi = dR * cpts'
                if abs(det(dxdxi))<1e-12
                    @warn "Singularity in mapping at $phys_pt"
                    dR = pinv(dxdxi)*dR
                else
                    dR = dxdxi\dR
                end
                Jac_par_phys = det(dxdxi)
                RR /= w_sum
                phys_pt = cpts*RR

                B = zeros(2*numNodes,3);
                B[1:2:2*numNodes-1,1] = dR[1,:]
                B[2:2:2*numNodes,2] = dR[2,:]
                B[1:2:2*numNodes-1,3] = dR[2,:]
                B[2:2:2*numNodes,3] = dR[1,:]
                solValX = RR'*sol0[2*curNodes.-1]
                solValY = RR'*sol0[2*curNodes]
                stressVect = Cmat*B'*sol0[curNodesXY]   
               
                exSolVal = exactSolDisp(phys_pt[1], phys_pt[2])
                exStressVal = exactSolStress(phys_pt[1], phys_pt[2])
                localArea = Jac_par_phys * Jac_ref_par *
                    gauss_rule[1].weights[iGauss] * gauss_rule[2].weights[jGauss]
                domainArea += localArea
                l2NormErr += localArea * ((exSolVal[1]-solValX)^2 + (exSolVal[2]-solValY)^2)
                l2NormSol += localArea * (exSolVal[1]^2 + exSolVal[2]^2)
                
                h1NormErr += localArea * (exStressVal-stressVect)'*invC*(exStressVal-stressVect)
                h1NormSol += localArea * exStressVal'*invC*exStressVal
            end
        end
    end
    relL2Err = sqrt(l2NormErr/l2NormSol)
    relH1Err = sqrt(h1NormErr/h1NormSol)
    return real(relL2Err), real(relH1Err)
end