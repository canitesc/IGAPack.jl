include("bernstein.jl")
using Plots
using LaTeXStrings
using WriteVTK
plotlyjs()

mutable struct Mesh
    elemVertex::Array{Float64,2}
    elemNode::Array{Array{Int64,1},1}
    degP::Array{Int64,1}
    C::Array{Array{Float64,2},1}
    numBasis::Int64
    numElem::Int64
    #localKnots::Array{Array{Float64,1},1}
    controlPoints::Array{Float64, 2}
    weights::Array{Float64, 1}
end

mutable struct Neighbors
    left::Array{Int64,1}
    right::Array{Int64,1}
end

"""
create the IEN (node index to element) array for a given knot vector and p and
elementVertex array

INPUT: knotVector - the given knot vector
        numElem - number of non-zero knot-spans
        p - polynomial degree
OUTPUT: IEN - array with nb rows and p+1 columns, where each row indicates the
                global node indices for a knot-span
        elemVertex - array with nb rows and 2 columns, where each row indicates
        the left and right knot of a non-empty knot-span
"""
function makeIEN(knotVector, numElem, p)
    IEN  = zeros(Int64, numElem, p+1)
    tolEq = 1e-10
    elemVertex = zeros(numElem, 2)
    elementCounter = 0
    for indexKnot = 1:length(knotVector)-1
        if knotVector[indexKnot+1]>knotVector[indexKnot]+tolEq
            elementCounter += 1
            IEN[elementCounter,:] = indexKnot-p:indexKnot
            elemVertex[elementCounter, :] = [knotVector[indexKnot], knotVector[indexKnot+1]]
        end
    end
    @assert numElem==elementCounter "Wrong number of elements passed"
    return IEN, elemVertex
end

function makeIEN(knotVectors::Array{Array{Float64,1},1}, numElem::Integer, p::Array{Int64, 1})
    numDim = length(knotVectors)
    numEntries = prod(p.+1)
    IEN = zeros(Int64, numElem, numEntries)
    tolEq = 1e-10
    elementVertex = zeros(numElem, 2^numDim)
    elementCounter = 0
    if numDim==2
        knotU = knotVectors[1]
        knotV = knotVectors[2]
        lenU = length(knotU)-p[1]-1
        for j=1:length(knotV)-1
            for i=1:length(knotU)-1
                if (knotU[i+1]>knotU[i]) && (knotV[j+1]>knotV[j])
                    elementCounter += 1
                    elementVertex[elementCounter, :] = [knotU[i], knotV[j], knotU[i+1], knotV[j+1]]

                    # now we add the nodes from i-p,..., i in the u direction
                    # j-q,..., j in the v direction
                    tcount = 0
                    for t2=j-degP[2]:j
                        for t1=i-degP[1]:i
                            tcount += 1
                            IEN[elementCounter, tcount] = t1+(t2-1)*lenU
                        end
                    end
                end
            end
        end
    end
    @assert numElem==elementCounter "Wrong number of elements passed"
    return IEN, elementVertex
end



"""
Initializes the elementNode, and elementChildren arrays
INPUT: knotVector - the given knot vector
        numElem - number of non-zero knot-spans
OUTPUT:
        elementNeighbors - array with nb elements of type Neighbors
        elementChildren - array with nb elements of type Children
"""
function initNeighborsChildren(knotVector, numElem)
    tolEq = 1e-10
    elementNeighbors = Array{Neighbors,1}(undef, numElem)
    elementChildren = Array{Array{Int64,1},1}(undef, numElem)
    elementCounter = 0
    for indexKnot = 1:length(knotVector)-1
        if knotVector[indexKnot+1]>knotVector[indexKnot]+tolEq
            elementCounter = elementCounter + 1
            elementNeighbors[elementCounter] = Neighbors(Array{Int64,1}(undef,1),
                                                         Array{Int64,1}(undef,1))

            if elementCounter>1
                elementNeighbors[elementCounter].left = [elementCounter-1]
            else
                elementNeighbors[elementCounter].left = []
            end
            if elementCounter<numElem
                elementNeighbors[elementCounter].right = [elementCounter+1]
            else
                elementNeighbors[elementCounter].right = []
            end
            elementChildren[elementCounter] = []
        end
    end
    @assert numElem==elementCounter "Wrong number of elements passed"
    return elementNeighbors, elementChildren
end

"""
Initialize the local knots corresponding the basis functions determined by the
knotVector with degree degP.
INPUT:
    knotVector - the given knot vector
    degP - polynomial degree
OUTPUT:
    localKnots: array of local knots vectors
"""
function initLocalKnots(knotVector, degP)
    numBasis = length(knotVector) - degP - 1
    localKnots = Array{Array{Float64,1}}(undef, numBasis)
    splineChildren = Array{Array{Float64,1}}(undef, numBasis)
    indexKnot = 1
    for i = 1:numBasis
        localKnots[i]= knotVector[i:i+degP+1]
    end
    return localKnots
end

"""
Initializes an IGA mesh from a NURBS object
"""
function genMesh(nurbs::NURBS)
    if length(nurbs.knots)==1
        #1D mesh
        knotU = nurbs.knots[1]
        degP = nurbs.order.-1
        Cmat, numElem = bezierExtraction(knotU, degP[1])
        C = Array{Array{Float64,2},1}(undef, numElem)
        for i=1:numElem
            C[i] = Cmat[:,:,i]
        end
        numBasis = length(knotU) - degP[1] - 1
        IEN, elemVertex = makeIEN(knotU, numElem, degP[1])

        #localKnots = initLocalKnots(knotU, degP)

    elseif length(nurbs.knots)==2
        #2D mesh
        knotU = nurbs.knots[1]
        knotV = nurbs.knots[2]
        degP = nurbs.order.-1
        CmatU, numElemU = bezierExtraction(knotU, degP[1])
        CmatV, numElemV = bezierExtraction(knotV, degP[2])
        numElem = numElemU * numElemV
        C = Array{Array{Float64,2},1}(undef, numElem)
        indexMatrix = permutedims(reshape(1:numElem, numElemU, numElemV), [2,1])
        for j=1:numElemV
            for i=1:numElemU
                elementIndex = indexMatrix[j,i]
                C[elementIndex] = kron(CmatV[:,:,j], CmatU[:,:,i])
            end
        end
        IEN, elemVertex = makeIEN([knotU, knotV], numElem, degP)
        numBasis = (length(knotU) - degP[1] - 1)*(length(knotV) - degP[2] - 1)
    end
    cpts = reshape(nrb.coefs[1:3,:,:,:], 3, numBasis)
    wgts = reshape(nrb.coefs[4,:,:,:], numBasis)
    for i=1:3
        cpts[i,:] ./= wgts
    end
    elemNode = Array{Array{Int64,1},1}(undef, numElem)
    for i=1:numElem
        elemNode[i]=IEN[i,:]
    end
    IGAmesh = Mesh(elemVertex, elemNode, degP, C, numBasis, numElem, cpts, wgts)
    return IGAmesh
end
"""
Plots the basis functions of a mesh in parameter space
"""
function plotBasisParam(mesh::Mesh)
    if length(degP)==1
        #1D plot
        colorList = ["blue", "red", "green", "black", "magenta"]
        graph=Plots.plot(title="B-Splines of degree $(mesh.degP[1])")
        numPtsElem = 11
        evalPts = LinRange(-1, 1, numPtsElem)
        B, dB = bernsteinBasis(evalPts, degP[1])

        for iBasis = 1:mesh.numBasis
            colorIndex = ((iBasis-1) % length(colorList))+1
            for iElem = 1:mesh.numElem
                localIndex = findall(isequal(iBasis), mesh.elemNode[iElem])
                if length(localIndex)>0
                    uMin = mesh.elemVertex[iElem, 1]
                    uMax = mesh.elemVertex[iElem, 2]
                    plotPts = LinRange(uMin, uMax, numPtsElem)
                    plotVal = B*(mesh.C[iElem][localIndex,:])'
                    graph = plot!(plotPts, plotVal, color=colorList[colorIndex], leg=false, line=2)
                    graph = scatter!([uMin], [0], color="red", leg=false, markersize=5)
                    graph = scatter!([uMax], [0], color="red", leg=false, markersize=5)
                end
            end
        end
    elseif length(degP)==2
        #2D plot
        graph=plot(title="B-Splines of degree $(mesh.degP)")
        numPtsElem = 11
        evalPtsU = LinRange(-1, 1, numPtsElem)
        evalPtsV = LinRange(-1, 1, numPtsElem)
        Bu, dBu = bernsteinBasis(evalPtsU, degP[1])
        Bv, dBv = bernsteinBasis(evalPtsV, degP[2])
        basisCounter = 0
        Buv = zeros(numPtsElem, numPtsElem, (degP[1]+1)*(degP[2]+1))
        for j=1:degP[2]+1
            for i=1:degP[1]+1
                basisCounter += 1
                Buv[:,:,basisCounter] = Bu[:,i]*Bv[:,j]'
            end
        end

        for iBasis = 1:mesh.numBasis
            for iElem = 1:mesh.numElem
                localIndex = findall(isequal(iBasis), mesh.elemNode[iElem])
                if length(localIndex)>0
                    uMin = mesh.elemVertex[iElem, 1]
                    uMax = mesh.elemVertex[iElem, 3]
                    vMin = mesh.elemVertex[iElem, 2]
                    vMax = mesh.elemVertex[iElem, 4]
                    plotPtsU = LinRange(uMin, uMax, numPtsElem)
                    plotPtsV = LinRange(vMin, vMax, numPtsElem)
                    cR = zeros(numPtsElem, numPtsElem)
                    for j=1:numPtsElem
                        for i=1:numPtsElem
                            cR[i,j] = (mesh.C[iElem][localIndex,:]*Buv[i,j,:])[1]
                        end
                    end


                end
            end
        end
    end

    display(graph)
end

"""
Plots the computed solution
"""
function plotSol(mesh::Mesh, sol0, fileName::String)

    if length(degP)==1
        graph=plot(title="Approximate solution")
        numPtsElem = 11
        evalPts = LinRange(-1, 1, numPtsElem)
        B, dB = bernsteinBasis(evalPts, degP[1])

        for iElem in 1:mesh.numElem
            curNodes = mesh.elemNode[iElem]
            uMin = mesh.elemVertex[iElem, 1]
            uMax = mesh.elemVertex[iElem, 2]
            plotPts = LinRange(uMin, uMax, numPtsElem)
            splineVal = B*(mesh.C[iElem])'
            cpts = mesh.controlPoints[1, curNodes]
            wgts = mesh.weights[curNodes]
            basisVal = zero(splineVal)
            for iPlotPt = 1:numPtsElem
                RR = splineVal[iPlotPt,:].* wgts
                w_sum = sum(RR)
                RR /= w_sum
                basisVal[iPlotPt,:] = RR
            end
            solVal = basisVal*sol0[curNodes]
            graph = plot!(plotPts, solVal, color="blue", leg=false, line=2)
            graph = plot!(plotPts, zeros(length(plotPts)), color="black", leg=false, line=1)
            graph = scatter!([uMin], [0], color="red", leg=false, markersize=5)
            graph = scatter!([uMax], [0], color="red", leg=false, markersize=5)
        end
        display(graph)
    elseif length(degP)==2
        #Set the order of points for VTK_QUADRATIC_QUAD  w.r.t. a 3x3 grid from Figure 3 of
        # https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
        VTKOrder = [1, 3, 9, 7, 2, 6, 8, 4]
        numPtsElem = 3
        evalPtsU = LinRange(-1, 1, numPtsElem)
        evalPtsV = LinRange(-1, 1, numPtsElem)
        Bu, dBu = bernsteinBasis(evalPtsU, degP[1])
        Bv, dBv = bernsteinBasis(evalPtsV, degP[2])
        basisCounter = 0
        Buv = zeros(numPtsElem, numPtsElem, (degP[1]+1)*(degP[2]+1))
        for j=1:degP[2]+1
            for i=1:degP[1]+1
                basisCounter += 1
                Buv[:,:,basisCounter] = Bu[:,i]*Bv[:,j]'
            end
        end
        pointsVTK = zeros(2, mesh.numElem*8)
        pointValVTK = zeros(1, mesh.numElem*8)
        localPoints = zeros(2, 9)
        localPointVal = zeros(9)
        meshCells = Array{MeshCell, 1}(undef, mesh.numElem)
        for iElem in 1:mesh.numElem
            curNodes = mesh.elemNode[iElem]
            uMin = mesh.elemVertex[iElem, 1]
            uMax = mesh.elemVertex[iElem, 3]
            vMin = mesh.elemVertex[iElem, 2]
            vMax = mesh.elemVertex[iElem, 4]
            cpts = mesh.controlPoints[1:2, curNodes]
            wgts = mesh.weights[curNodes]
            localPtCounter = 0
            for jPt = 1:numPtsElem
                for iPt = 1:numPtsElem
                    #compute the (B-)spline basis functions and derivatives with Bezier extraction
                    N_mat = mesh.C[iElem] * Buv[iPt, jPt, :]
                    #compute the rational basis
                    RR = N_mat.* wgts
                    w_sum = sum(RR)
                    RR /= w_sum
                    phys_pt = cpts*RR
                    localPtCounter += 1
                    localPoints[:, localPtCounter ]=phys_pt
                    solVal = RR'*sol0[curNodes]
                    localPointVal[localPtCounter]=solVal
                end
            end
            ptRange = (iElem-1)*8+1:iElem*8
            pointsVTK[:, ptRange] = localPoints[:, VTKOrder]
            pointValVTK[:, ptRange] = localPointVal[VTKOrder]
            meshCells[iElem] = MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD, collect(ptRange))
        end
        #Output data to VTK
        vtkfile = vtk_grid(fileName, pointsVTK, meshCells)
        vtkfile["U", VTKPointData()] = pointValVTK
        vtk_save(vtkfile)
        println("Output written to $fileName.vtk")
    end

end

"""
Plots the computed solution for elasticity
"""
function plotSolElast(mesh::Mesh, Cmat::Array{Float64}, sol0, fileName::String)    
    if length(degP)==2
        #Set the order of points for VTK_QUADRATIC_QUAD  w.r.t. a 3x3 grid from Figure 3 of
        # https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
        VTKOrder = [1, 3, 9, 7, 2, 6, 8, 4]
        numPtsElem = 3
        evalPtsU = LinRange(-1, 1, numPtsElem)
        evalPtsV = LinRange(-1, 1, numPtsElem)
        Buv, dBdu, dBdv = bernsteinBasis2D(collect(evalPtsU), collect(evalPtsV), mesh.degP)
        pointsVTK = zeros(2, mesh.numElem*8)
        pointDispVTK = zeros(2, mesh.numElem*8)
        pointStressVTK = zeros(3, mesh.numElem*8)
        pointStressVMVTK = zeros(mesh.numElem*8)
        localPoints = zeros(2, 9)
        localPointDisp = zeros(2, 9)
        localPointStress = zeros(3, 9)
        localPointStressVM = zeros(1, 9)
        meshCells = Array{MeshCell, 1}(undef, mesh.numElem)
        for iElem in 1:mesh.numElem
            curNodes = mesh.elemNode[iElem]
            numNodes = length(curNodes)
            curNodesXY = reshape(hcat(2*curNodes.-1, 2*curNodes)', 2*numNodes)
            uMin = mesh.elemVertex[iElem, 1]
            uMax = mesh.elemVertex[iElem, 3]
            vMin = mesh.elemVertex[iElem, 2]
            vMax = mesh.elemVertex[iElem, 4]
            cpts = mesh.controlPoints[1:2, curNodes]
            wgts = mesh.weights[curNodes]
            localPtCounter = 0
            B = zeros(2*numNodes,3)
            for jPt = 1:numPtsElem
                for iPt = 1:numPtsElem
                    #compute the (B-)spline basis functions and derivatives with Bezier extraction
                    N_mat = mesh.C[iElem] * Buv[iPt, jPt, :]
                    dN_du = mesh.C[iElem] * dBdu[iPt, jPt, :] * 2/(uMax-uMin)
                    dN_dv = mesh.C[iElem] * dBdv[iPt, jPt, :] * 2/(vMax-vMin)
                    
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
                    RR /= w_sum
                    phys_pt = cpts*RR
                    if abs(det(dxdxi))<1e-12
                        @warn "Singularity in mapping at $phys_pt"
                        dR = pinv(dxdxi)*dR
                    else
                        dR = dxdxi\dR
                    end
                    Jac_par_phys = det(dxdxi)
                                        
                    B[1:2:2*numNodes-1,1] = dR[1,:]
                    B[2:2:2*numNodes,2] = dR[2,:]
                    B[1:2:2*numNodes-1,3] = dR[2,:]
                    B[2:2:2*numNodes,3] = dR[1,:]
                    localPtCounter += 1
                    localPoints[:, localPtCounter ]=phys_pt
                    solValX = RR'*sol0[2*curNodes.-1]
                    solValY = RR'*sol0[2*curNodes]
                    localPointDisp[:, localPtCounter]=[solValX, solValY]
                    stressVect = Cmat*B'*sol0[curNodesXY]                    
                    localPointStress[:, localPtCounter] = stressVect
                    stressVM = sqrt(stressVect[1]^2 - stressVect[1]*stressVect[2] + stressVect[2]^2 
                                    +3*stressVect[3]^2)                                   
                    localPointStressVM[localPtCounter] = stressVM
                end
            end
            ptRange = (iElem-1)*8+1:iElem*8
            pointsVTK[:, ptRange] = localPoints[:, VTKOrder]
            pointDispVTK[:, ptRange] = localPointDisp[:, VTKOrder]
            pointStressVTK[:, ptRange] = localPointStress[:, VTKOrder]
            pointStressVMVTK[ptRange] = localPointStressVM[VTKOrder]
            meshCells[iElem] = MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD, collect(ptRange))
        end
        #Output data to VTK
        vtkfile = vtk_grid(fileName, pointsVTK, meshCells)      
        vtkfile["U", VTKPointData()] = pointDispVTK
        vtkfile["Stress", VTKPointData()] = pointStressVTK
        vtkfile["StressVM"] = pointStressVMVTK        
        vtk_save(vtkfile)
        println("Output written to $fileName.vtk")
    end

end

"""
Plots the error in the approximate solution
"""
function plotSolError(mesh::Mesh, sol0, exactSol::Function, fileName::String)
    if length(degP)==1
        graph=plot(title="Error \$u-u_h\$")
        numPtsElem = 11
        evalPts = LinRange(-1, 1, numPtsElem)
        B, dB = bernsteinBasis(evalPts, degP[1])

        for iElem in 1:mesh.numElem
            curNodes = mesh.elemNode[iElem]
            uMin = mesh.elemVertex[iElem, 1]
            uMax = mesh.elemVertex[iElem, 2]
            plotPts = LinRange(uMin, uMax, numPtsElem)
            physPts = zeros(numPtsElem)
            splineVal = B*(mesh.C[iElem])'
            cpts = mesh.controlPoints[1, curNodes]
            wgts = mesh.weights[curNodes]
            basisVal = zero(splineVal)
            for iPlotPt = 1:numPtsElem
                RR = splineVal[iPlotPt,:].* wgts
                w_sum = sum(RR)
                RR /= w_sum
                physPts[iPlotPt] = RR'*cpts
                basisVal[iPlotPt,:] = RR
            end
            exSolVal = real(exactSol.(physPts))
            solVal = basisVal*sol0[curNodes]
            graph = plot!(plotPts, exSolVal-solVal, color="blue", leg=false, line=2)
            graph = plot!(plotPts, zeros(length(plotPts)), color="black", leg=false, line=1)
            graph = scatter!([uMin], [0], color="red", leg=false, markersize=5)
            graph = scatter!([uMax], [0], color="red", leg=false, markersize=5)
        end
        display(graph)
    elseif length(degP)==2
        #Set the order of points for VTK_QUADRATIC_QUAD  w.r.t. a 3x3 grid from Figure 3 of
        # https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
        VTKOrder = [1, 3, 9, 7, 2, 6, 8, 4]
        numPtsElem = 3
        evalPtsU = LinRange(-1, 1, numPtsElem)
        evalPtsV = LinRange(-1, 1, numPtsElem)
        Bu, dBu = bernsteinBasis(evalPtsU, degP[1])
        Bv, dBv = bernsteinBasis(evalPtsV, degP[2])
        basisCounter = 0
        Buv = zeros(numPtsElem, numPtsElem, (degP[1]+1)*(degP[2]+1))
        for j=1:degP[2]+1
            for i=1:degP[1]+1
                basisCounter += 1
                Buv[:,:,basisCounter] = Bu[:,i]*Bv[:,j]'
            end
        end
        pointsVTK = zeros(2, mesh.numElem*8)
        pointValVTK = zeros(1, mesh.numElem*8)
        localPoints = zeros(2, 9)
        localPointVal = zeros(9)
        meshCells = Array{MeshCell, 1}(undef, mesh.numElem)
        for iElem in 1:mesh.numElem
            curNodes = mesh.elemNode[iElem]
            uMin = mesh.elemVertex[iElem, 1]
            uMax = mesh.elemVertex[iElem, 3]
            vMin = mesh.elemVertex[iElem, 2]
            vMax = mesh.elemVertex[iElem, 4]
            cpts = mesh.controlPoints[1:2, curNodes]
            wgts = mesh.weights[curNodes]
            localPtCounter = 0
            for jPt = 1:numPtsElem
                for iPt = 1:numPtsElem
                    #compute the (B-)spline basis functions and derivatives with Bezier extraction
                    N_mat = mesh.C[iElem] * Buv[iPt, jPt, :]
                    #compute the rational basis
                    RR = N_mat.* wgts
                    w_sum = sum(RR)
                    RR /= w_sum
                    phys_pt = cpts*RR
                    localPtCounter += 1
                    localPoints[:, localPtCounter ]=phys_pt

                    solVal = RR'*sol0[curNodes]
                    exSolVal = real(exactSol(phys_pt[1], phys_pt[2]))
                    localPointVal[localPtCounter]=exSolVal-solVal
                end
            end
            ptRange = (iElem-1)*8+1:iElem*8
            pointsVTK[:, ptRange] = localPoints[:, VTKOrder]
            pointValVTK[:, ptRange] = localPointVal[VTKOrder]
            meshCells[iElem] = MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD, collect(ptRange))
        end
        #Output data to VTK
        vtkfile = vtk_grid(fileName, pointsVTK, meshCells)
        vtkfile["U-U_h", VTKPointData()] = pointValVTK
        vtk_save(vtkfile)
        println("Output written to $fileName.vtk")
    end
end

"""
Plots the error in the approximate solution
"""
function plotSolErrorElast(mesh::Mesh, Cmat::Array{Float64}, sol0, exactSolDisp::Function, exactSolStress::Function, fileName::String)
    if length(degP)==2
         #Set the order of points for VTK_QUADRATIC_QUAD  w.r.t. a 3x3 grid from Figure 3 of
        # https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
        VTKOrder = [1, 3, 9, 7, 2, 6, 8, 4]
        numPtsElem = 3
        evalPtsU = LinRange(-1, 1, numPtsElem)
        evalPtsV = LinRange(-1, 1, numPtsElem)
        Buv, dBdu, dBdv = bernsteinBasis2D(collect(evalPtsU), collect(evalPtsV), mesh.degP)
        pointsVTK = zeros(2, mesh.numElem*8)
        pointDispVTK = zeros(2, mesh.numElem*8)
        pointStressVTK = zeros(3, mesh.numElem*8)
        pointStressVMVTK = zeros(mesh.numElem*8)
        localPoints = zeros(2, 9)
        localPointDisp = zeros(2, 9)
        localPointStress = zeros(3, 9)
        localPointStressVM = zeros(1, 9)
        meshCells = Array{MeshCell, 1}(undef, mesh.numElem)
        for iElem in 1:mesh.numElem
            curNodes = mesh.elemNode[iElem]
            numNodes = length(curNodes)
            curNodesXY = reshape(hcat(2*curNodes.-1, 2*curNodes)', 2*numNodes)
            uMin = mesh.elemVertex[iElem, 1]
            uMax = mesh.elemVertex[iElem, 3]
            vMin = mesh.elemVertex[iElem, 2]
            vMax = mesh.elemVertex[iElem, 4]
            cpts = mesh.controlPoints[1:2, curNodes]
            wgts = mesh.weights[curNodes]
            localPtCounter = 0
            for jPt = 1:numPtsElem
                for iPt = 1:numPtsElem
                    #compute the (B-)spline basis functions and derivatives with Bezier extraction
                    N_mat = mesh.C[iElem] * Buv[iPt, jPt, :]
                    dN_du = mesh.C[iElem] * dBdu[iPt, jPt, :] * 2/(uMax-uMin)
                    dN_dv = mesh.C[iElem] * dBdv[iPt, jPt, :] * 2/(vMax-vMin)
                    
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
                    RR /= w_sum
                    phys_pt = cpts*RR
                    if abs(det(dxdxi))<1e-12
                        @warn "Singularity in mapping at $phys_pt"
                        dR = pinv(dxdxi)*dR
                    else
                        dR = dxdxi\dR
                    end
                    Jac_par_phys = det(dxdxi)
                    
                    B = zeros(2*numNodes,3);
                    B[1:2:2*numNodes-1,1] = dR[1,:]
                    B[2:2:2*numNodes,2] = dR[2,:]
                    B[1:2:2*numNodes-1,3] = dR[2,:]
                    B[2:2:2*numNodes,3] = dR[1,:]
                    localPtCounter += 1
                    localPoints[:, localPtCounter ]=phys_pt
                    solValX = RR'*sol0[2*curNodes.-1]
                    solValY = RR'*sol0[2*curNodes]

                    #evaluate the exact displacement and stresses
                    exSolVal = exactSolDisp(phys_pt[1], phys_pt[2])                    
                    exSolStressVal = exactSolStress(phys_pt[1], phys_pt[2])
                    exStressVM = sqrt(exSolStressVal[1]^2 - exSolStressVal[1]*exSolStressVal[2] + exSolStressVal[2]^2 
                                    +3*exSolStressVal[3]^2)

                    #compute the errors (exact - approximate)
                    localPointDisp[:, localPtCounter]=[exSolVal[1]-solValX, exSolVal[2]-solValY]
                    stressVect = Cmat*B'*sol0[curNodesXY]                    
                    localPointStress[:, localPtCounter] = exSolStressVal - stressVect
                    stressVM = sqrt(stressVect[1]^2 - stressVect[1]*stressVect[2] + stressVect[2]^2 
                                    +3*stressVect[3]^2)                                   
                    localPointStressVM[localPtCounter] = exStressVM - stressVM
                end
            end
            ptRange = (iElem-1)*8+1:iElem*8
            pointsVTK[:, ptRange] = localPoints[:, VTKOrder]
            pointDispVTK[:, ptRange] = localPointDisp[:, VTKOrder]
            pointStressVTK[:, ptRange] = localPointStress[:, VTKOrder]
            pointStressVMVTK[ptRange] = localPointStressVM[VTKOrder]
            meshCells[iElem] = MeshCell(VTKCellTypes.VTK_QUADRATIC_QUAD, collect(ptRange))
        end
        #Output data to VTK
        vtkfile = vtk_grid(fileName, pointsVTK, meshCells)      
        vtkfile["error U", VTKPointData()] = pointDispVTK
        vtkfile["error Stress", VTKPointData()] = pointStressVTK
        vtkfile["error StressVM"] = pointStressVMVTK        
        vtk_save(vtkfile)
        println("Output written to $fileName.vtk")
    end
end