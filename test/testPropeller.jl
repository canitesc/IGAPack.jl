# script for generating the patches for a propeller model
using LinearAlgebra
include("../src/SplinePlotting.jl")
plotly()

function getLineProj(R, theta, theta1)
    #gets the coordinates of the points on the tangent to the arc
    pointCenter = zeros(2)
    pointRight = zeros(2)
    pointLeft = zeros(2)

    pointCenter[1] = R*cos(theta)
    pointCenter[2] = R*sin(theta)

    pointRight[1] = (R*sin(theta)+R*cot(theta)*cos(theta))/(tan(theta+theta1)+cot(theta))
    pointRight[2] = tan(theta+theta1)*pointRight[1]

    pointLeft[1] = (R*sin(theta)+R*cot(theta)*cos(theta))/(tan(theta-theta1)+cot(theta))
    pointLeft[2] = tan(theta-theta1)*pointLeft[1]
    return pointLeft, pointRight, pointCenter
end

#dimensions of the hub
radInt = 5 #inner radius
radExt = 6 #outer radius
hubHeight = 10 #height
centerOffset = 0.5 #amount by which the blade interface is lifted from the hub
P_D = 1.4 # P/D pitch distribution

numBlades = 8  #number of blades


bladeThicknessRatio = 0.02 #
bladeHeightRatio = 0.85 #
bladeRotAngle = [-0.8*(pi/2-atan(P_D/pi)),-0.88*(pi/2-atan(P_D/pi)),-0.96*(pi/2-atan(P_D/pi)), -(pi/2-atan(P_D/pi)),-(pi/2-atan(P_D/pi)), -(pi/2-atan(P_D/pi)), -(pi/2-atan(P_D/pi))]
###########BLADE SHAPE DATA#############
coeff1 = 0.2
coeff2 = 3.0
coeff3 = 0.1
coeff4 = 1.2
#top and bottom surface of the blade
coefs_ext = zeros(4,7,3)
coefs_int = zeros(4,7,3)
coefs_ext[:,:,1] = [  0   30   60   90   120  165  170
                      0    20*coeff4        29*coeff4         30*coeff4         24*coeff4        0*coeff4         -45
                      -0.75*coeff1  -0.65*coeff1  -0.50*coeff1  -0.35*coeff1  -0.20*coeff1  0.0*coeff1  0.0*coeff1
                      1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000]

coefs_ext[:,:,2] = [  0   30   60   90   120  165  170
                      -15  -13 -12  -12 -14  -40  -50
                       -0.75*coeff2  -0.65*coeff2  -0.50*coeff2  -0.35*coeff2  -0.20*coeff2  0.0*coeff2  0.0*coeff2
                      1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000]

coefs_ext[:,:,3] = [  0   30   60   90   120  165  170
                      -60  -67 -74  -81 -88  -89  -58
                       -0.75*coeff3  -0.65*coeff3  -0.50*coeff3  -0.35*coeff3  -0.20*coeff3  0.0*coeff3  0.0*coeff3
                      1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000]
#midsurface of the blade
coefs_int[:,:,1] = [ 0   30   60   90   120  165  170
                     0+2  20*coeff4+2  29*coeff4+2  30*coeff4+2  24*coeff4+2  0*coeff4+2  -45+2
                    -0.75*coeff1/2  -0.65*coeff1/2  -0.50*coeff1/2  -0.35*coeff1/2  -0.20*coeff1/2  0.0*coeff1/2  0.0*coeff1/2
                     1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000]

coefs_int[:,:,2] = [ 0   30   60   90   120  165  175
                      -15  -13 -12  -12 -14  -40  -50
                     -0.75*coeff2/2  -0.65*coeff2/2  -0.50*coeff2/2  -0.35*coeff2/2  -0.20*coeff2/2  0.0*coeff2/2  0.0*coeff2/2
                     1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];

coefs_int[:,:,3] = [ 0   30   60   90   120  165  170
                     -60-2  -67-2 -74-2  -81-2 -88-2  -89-2  -58-2
                     -0.75*coeff3/2  -0.65*coeff3/2  -0.50*coeff3/2  -0.35*coeff3/2  -0.20*coeff3/2  0.0*coeff3/2  0.0*coeff3/2
                     1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
#####END BLADE SHAPE DATA######



arcBlade = 2*pi/numBlades
centerBottom  = [0, 0, 0]

bladeHeightBottom = hubHeight*(1-bladeHeightRatio)/2
bladeHeightMiddle = hubHeight*bladeHeightRatio
bladeHeightTop = bladeHeightBottom
bladeBaseArcRatio = 0.5
bladeBaseBottom = bladeHeightBottom*(1-bladeBaseArcRatio)
bladeBaseTop = bladeHeightBottom + bladeHeightMiddle + bladeBaseArcRatio*bladeHeightTop
bladeCenterHeight = hubHeight/2


arcLeft = arcBlade*(1-bladeThicknessRatio)/2
arcRight = arcLeft
arcMiddle = arcBlade*bladeThicknessRatio

#define the arcs
nrbArcLeftExt = nrbcirc(radExt, centerBottom, 0, arcLeft)
nrbArcMiddleExt = nrbcirc(radExt, centerBottom, arcLeft, arcLeft+arcMiddle)
nrbArcRightExt = nrbcirc(radExt, centerBottom, arcLeft+arcMiddle, arcLeft+arcMiddle+arcRight)

nrbArcLeftInt = nrbcirc(radInt, centerBottom, 0, arcLeft)
nrbArcMiddleInt = nrbcirc(radInt, centerBottom, arcLeft, arcLeft+arcMiddle)
nrbArcRightInt = nrbcirc(radInt, centerBottom, arcLeft+arcMiddle, arcLeft+arcMiddle+arcRight)

knotExt = [0,0,0,1/3,1/3,2/3,2/3,1,1,1]
coefsExt = hcat(nrbArcLeftExt.coefs, nrbArcMiddleExt.coefs[:,2:end], nrbArcRightExt.coefs[:,2:end])
nrbArcExt = nrbmak(coefsExt, [knotExt])

knotInt = [0,0,0,1/3,1/3,2/3,2/3,1,1,1]
coefsInt = hcat(nrbArcLeftInt.coefs, nrbArcMiddleInt.coefs[:,2:end], nrbArcRightInt.coefs[:,2:end])
nrbArcInt = nrbmak(coefsInt, [knotInt])

knotThickness = [0,0,1,1]
coefsAnnulus = zeros(4,7,2)
coefsAnnulus[:,:,1] = coefsInt
coefsAnnulus[:,:,2] = coefsExt
nrbAnnulus = nrbmak(coefsAnnulus,[knotInt,knotThickness])

# extrude the annuli along the height
layerSet = [bladeBaseBottom, bladeHeightBottom, bladeCenterHeight, bladeHeightBottom+bladeHeightMiddle,  bladeBaseTop, hubHeight]
knotHeight = [0,0,0, 1/3, 1/3, 2/3, 2/3, 1, 1,1]

# ajust the z-coordinate for each layer
coefsHub = zeros(size(coefsAnnulus)...,length(layerSet)+1...)
coefsHub[:,:,:,1] = coefsAnnulus
for indexSlice = 1:length(layerSet)
    currentSlice = layerSet[indexSlice]
    coefsHub[:,:,:,indexSlice+1] = coefsAnnulus
    coefsHub[3,:,:,indexSlice+1] = currentSlice.*coefsAnnulus[4,:,:]
end
nrbHub = nrbmak(coefsHub,[knotInt, knotThickness, knotHeight])

# % hold on
theta = arcBlade/2
theta1 = theta - arcLeft


pointLeft, pointRight, pointCenter = getLineProj(radExt+centerOffset, theta, theta1)

#define the rotation axis and the rotation matrix
ux = cos(theta)
uy = sin(theta)
uz = 0
rotMatrix = [cos(bladeRotAngle[1])+ux^2*(1-cos(bladeRotAngle[1])) ux*uy*(1-cos(bladeRotAngle[1]))-uz*sin(bladeRotAngle[1]) ux*uz*(1-cos(bladeRotAngle[1]))+uy*sin(bladeRotAngle[1])
    uy*ux*(1-cos(bladeRotAngle[1]))+uz*sin(bladeRotAngle[1]) cos(bladeRotAngle[1])+uy^2*(1-cos(bladeRotAngle[1])) uy*uz*(1-cos(bladeRotAngle[1]))-ux*sin(bladeRotAngle[1])
    uz*ux*(1-cos(bladeRotAngle[1]))-uy*sin(bladeRotAngle[1]) uz*uy*(1-cos(bladeRotAngle[1]))+ux*sin(bladeRotAngle[1]) cos(bladeRotAngle[1])+uz^2*(1-cos(bladeRotAngle[1]))]

bladeInterfacePoints = [pointLeft pointCenter pointRight pointLeft pointCenter pointRight pointLeft pointCenter pointRight
                        ones(1,3)*layerSet[2] ones(1,3)*layerSet[3] ones(1,3)*layerSet[4]]
centerPoint = [pointCenter' layerSet[3]]
bladeInterfacePoints[1,:] = bladeInterfacePoints[1,:] .- centerPoint[1]
bladeInterfacePoints[2,:] = bladeInterfacePoints[2,:] .- centerPoint[2]
bladeInterfacePoints[3,:] = bladeInterfacePoints[3,:] .- centerPoint[3]
bladeInterfacePoints = rotMatrix*bladeInterfacePoints
bladeInterfacePoints[1,:] = bladeInterfacePoints[1,:] .+ centerPoint[1]
bladeInterfacePoints[2,:] = bladeInterfacePoints[2,:] .+ centerPoint[2]
bladeInterfacePoints[3,:] = bladeInterfacePoints[3,:] .+ centerPoint[3]

coefsNewHub = coefsHub

newCoefsLow = [bladeInterfacePoints[:,1:3]; ones(1,3)]
newCoefsMid = [bladeInterfacePoints[:,4:6]; ones(1,3)]
newCoefsTop = [bladeInterfacePoints[:,7:9]; ones(1,3)]
coefsNewHub[:,3:5,2,3] = newCoefsLow
coefsNewHub[:,3:5,2,4] = newCoefsMid
coefsNewHub[:,3:5,2,5] = newCoefsTop

nrbNewHub = nrbmak(coefsNewHub, [knotInt, knotThickness, knotHeight])
# Plots.surface(reuse = false)
# plt = nrbctrlplot(nrbNewHub)
# display(plt)

bladeWidth = norm(pointLeft-pointRight)
bladeHeight = layerSet[4]-layerSet[2]

coefs = coefs_ext
coefs_blade = zeros(size(coefs)...,3...)
scaleFactor = bladeHeight/60
coefs[1:2,:,:] = scaleFactor*coefs[1:2,:,:]
coefs_blade[:,:,:,1] = coefs
coefs[3,:,:] .= bladeWidth
coefs_blade[:,:,:,3] = coefs
coefs = coefs_int
coefs[1:2,:,:] = scaleFactor*coefs[1:2,:,:]
coefs[3,:,:] .= bladeWidth/2
coefs_blade[:,:,:,2] = coefs

numix = 5
numiy = 1
p=2
q=2


# initialize the knot vector
knotU = vcat(zeros(p),collect(LinRange(0,1,numix+1)),ones(p))
knotV = vcat(zeros(q),collect(LinRange(0,1,numiy+1)),ones(q))
knotW = vcat(zeros(q+1), ones(q+1))

# rotate about the x-axis by pi/2
thetaX = pi/2
rotMatrixX = [1  0  0; 0  cos(thetaX)  -sin(thetaX); 0  sin(thetaX)  cos(thetaX)]
for i=1:size(coefs_blade,3)
    for j=1:size(coefs_blade,4)
        coefs_blade[1:3,:,i,j] = rotMatrixX*coefs_blade[1:3,:,i,j]
    end
end

# translate back (in the positive direction of y-axis) by bladeWidth/2
coefs_blade[2,:,:,:] = coefs_blade[2,:,:,:].+bladeWidth/2
# rotate about the z-axis by theta
rotMatrixZ = [cos(theta)  -sin(theta)  0; sin(theta)  cos(theta)  0; 0  0  1]
for i=1:size(coefs_blade,3)
    for j=1:size(coefs_blade,4)
        coefs_blade[1:3,:,i,j] = rotMatrixZ*coefs_blade[1:3,:,i,j]
    end
end


# shift to the center of the target blade location
coefs_blade[3,:,:,:] = coefs_blade[3,:,:,:] .+ bladeHeight/2;
# rotate the blade by bladeRot
for i=1:size(coefs_blade,3)
    for j=1:size(coefs_blade,4)
        for k=1:size(coefs_blade,2)
            global rotMatrix
            R11 = cos(bladeRotAngle[k])+ux^2*(1-cos(bladeRotAngle[k]))
            R12 = ux*uy*(1-cos(bladeRotAngle[k]))-uz*sin(bladeRotAngle[k])
            R13 = ux*uz*(1-cos(bladeRotAngle[k]))+uy*sin(bladeRotAngle[k])
            R21 = uy*ux*(1-cos(bladeRotAngle[k]))+uz*sin(bladeRotAngle[k])
            R22 = cos(bladeRotAngle[k])+uy^2*(1-cos(bladeRotAngle[k]))
            R23 = uy*uz*(1-cos(bladeRotAngle[k]))-ux*sin(bladeRotAngle[k])
            R31 = uz*ux*(1-cos(bladeRotAngle[k]))-uy*sin(bladeRotAngle[k])
            R32 = uz*uy*(1-cos(bladeRotAngle[k]))+ux*sin(bladeRotAngle[k])
            R33 = cos(bladeRotAngle[k])+uz^2*(1-cos(bladeRotAngle[k]))
            rotMatrix = [ R11 R12 R13; R21 R22 R23; R31 R32 R33]
            coefs_blade[1:3,k,i,j] = rotMatrix*coefs_blade[1:3,k,i,j]
        end
    end
end

# shift to the target location
coefs_blade[1,:,:,:] = coefs_blade[1,:,:,:] .+ centerPoint[1]
coefs_blade[2,:,:,:] = coefs_blade[2,:,:,:] .+ centerPoint[2]
coefs_blade[3,:,:,:] = coefs_blade[3,:,:,:] .+ centerPoint[3]
nrbBlade = nrbmak(coefs_blade, [knotU, knotV,knotW]);

# Plots.surface(reuse = false)
# plt = nrbctrlplot(nrbBlade)
# display(plt)

#rotate nrbNewHub and nrbBlade to get the full propeller
Plots.surface(reuse = false)

hubPatches = Array{NURBS, 1}(undef, numBlades)
bladePatches = Array{NURBS, 1}(undef, numBlades)
for i=1:numBlades
    global plt
    global rotMatrixZ
    rotAngle = arcBlade*(i-1);
    rotMatrixZ = [cos(rotAngle) -sin(rotAngle) 0; sin(rotAngle) cos(rotAngle) 0; 0 0 1];
    #rotate and plot the blade
    coefs_blade_current = coefs_blade;
    for j=1:size(coefs_blade,3)
        for k=1:size(coefs_blade,4)
            coefs_blade_current[1:3,:,j,k] = rotMatrixZ*coefs_blade[1:3,:,j,k]
        end
    end
    bladePatches[i] = nrbmak(coefs_blade_current, [knotU, knotV, knotW])
    plt = nrbkntplot(bladePatches[i])
    #plt = nrbplot(bladePatches{i},[20,20,2])

    #rotate and plot the hub
    coefsNewHubCurrent = coefsNewHub;
    for j=1:size(coefsNewHub,3)
        for k=1:size(coefsNewHub,4)
            coefsNewHubCurrent[1:3,:,j,k] = rotMatrixZ*coefsNewHub[1:3,:,j,k]
        end
    end
    hubPatches[i] = nrbmak(coefsNewHubCurrent, [knotInt, knotThickness, knotHeight])
    plt = nrbkntplot(hubPatches[i])
    #plt = nrbplot(hubPatches[i],[20,2,20])
end
display(plt)

