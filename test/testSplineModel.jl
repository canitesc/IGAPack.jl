#tests for SplineModel.jl module
using Test
include("../src/SplineModel.jl")

#####tests for nrbmak()###############
#based on Demonstration 1
pnts = [0.5 1.5 4.5 3.0 7.5 6.0 8.5;
         3.0 5.5 5.5 1.5 1.5 4.0 4.5;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0]
knot_crv = [0, 0, 0, 1/4, 1/2, 3/4, 3/4, 1, 1, 1];
crv = nrbmak(pnts, [knot_crv])
@test crv.number == [7]
@test crv.coefs == [pnts; ones(1, crv.number[1])]
@test crv.knots == [knot_crv]
@test crv.order == [3]

#based on Demonstration 2
pnts = zeros(3,5,5);
pnts[:,:,1] = [ 0.0  3.0  5.0  8.0 10.0;
             0.0  0.0  0.0  0.0  0.0;
             2.0  2.0  7.0  7.0  8.0];
pnts[:,:,2] = [ 0.0  3.0  5.0  8.0 10.0;
             3.0  3.0  3.0  3.0  3.0;
             0.0  0.0  5.0  5.0  7.0];
pnts[:,:,3] = [ 0.0  3.0  5.0  8.0 10.0;
             5.0  5.0  5.0  5.0  5.0;
             0.0  0.0  5.0  5.0  7.0];
pnts[:,:,4] = [ 0.0  3.0  5.0  8.0 10.0;
             8.0  8.0  8.0  8.0  8.0;
             5.0  5.0  8.0  8.0 10.0];
pnts[:,:,5] = [ 0.0  3.0  5.0  8.0 10.0;
            10.0 10.0 10.0 10.0 10.0;
             5.0  5.0  8.0  8.0 10.0];
uknots = [0, 0, 0, 1/3, 2/3, 1, 1, 1]
vknots = [0, 0, 0, 1/3, 2/3, 1, 1, 1]
knots = [uknots, vknots]
srf = nrbmak(pnts,knots);
@test srf.number == [5, 5]
@test srf.coefs == [pnts; ones(1,5,5)]
@test srf.knots == [uknots, vknots]
@test srf.order == [3, 3]

#based on Demonstration 3
coefs =[ 6.0  0.0  6.0  1;
        -5.5  0.5  5.5  1;
        -5.0  1.0 -5.0  1;
        4.5  1.5 -4.5  1;
        4.0  2.0  4.0  1;
        -3.5  2.5  3.5  1;
        -3.0  3.0 -3.0  1;
        2.5  3.5 -2.5  1;
        2.0  4.0  2.0  1;
        -1.5  4.5  1.5  1;
        -1.0  5.0 -1.0  1;
        0.5  5.5 -0.5  1;
        0.0  6.0  0.0  1]';
knots = [0, 0, 0, 0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1, 1, 1];

crv3d = nrbmak(coefs,[knots]);
@test crv3d.number == [13]
@test crv3d.coefs == coefs
@test crv3d.knots ==[knots]
@test crv3d.order == [4]

##############tests for findspan()###############

n = 3
U = [0, 0, 0, 1/2, 1, 1, 1]
u = LinRange(0, 1, 10)
s = findspan(n, u, U)
@test s == [2*ones(5); 3*ones(5)]

p = 2
m = 7
n = m - p - 1
U = [zeros(p); LinRange(0,1,m+1-2*p); ones(p)]
u = [0, 0.11880, 0.55118, 0.93141, 0.40068, 0.35492, 0.44392, 0.88360, 0.35414, 0.92186, 0.83085, 1]
s = [2, 2, 3, 4, 3, 3, 3, 4, 3, 4, 4, 4]
@test findspan(n, u, U)==s

#############test for basisfun()##################

n = 3
U = [0, 0, 0, 1/2, 1, 1, 1]
p = 2;
u = LinRange(0, 1, 10);
s = findspan(n, u, U);
Bref = [1.000000000000000                   0                   0
        0.604938271604938   0.370370370370370   0.024691358024691
        0.308641975308642   0.592592592592593   0.098765432098765
        0.111111111111111   0.666666666666667   0.222222222222222
        0.012345679012346   0.592592592592593   0.395061728395062
        0.395061728395062   0.592592592592593   0.012345679012346
        0.222222222222222   0.666666666666667   0.111111111111111
        0.098765432098765   0.592592592592593   0.308641975308642
        0.024691358024691   0.370370370370370   0.604938271604938
        0                   0   1.000000000000000]
B = basisfun(s, u, p, U)
@test B≈Bref
############test for nrbeval()###################
uknots = [0, 0, 0, 1, 1, 1]
vknots = [0, 0, 0, .5, 1, 1, 1]
wknots = [0, 0, 0, 0, 1, 1, 1, 1]
knots = [uknots, vknots, wknots]
cx = [0 0.5 1]
nx = length(cx)
cy = [0 0.25 0.75 1]
ny = length(cy)
cz = [0 1/3 2/3 1]
nz = length(cz)
coefs = zeros(4, nx, ny, nz)
coefs[1,:,:,:] = repeat(reshape(cx,nx,1,1),1, ny, nz)
coefs[2,:,:,:] = repeat(reshape(cy,1,ny,1), nx, 1, nz)
coefs[3,:,:,:] = repeat(reshape(cz,1,1,nz), nx, ny, 1)
coefs[4,:,:,:] .= 1
nurbs = nrbmak(coefs, knots)
x = rand(5,1)
y = rand(5,1)
z = rand(5,1)
tt = [x y z]'
points = nrbeval(nurbs,tt)
@test points ≈ tt

uknots = [0, 0, 0, 1, 1, 1]
vknots = [0, 0, 0, 0, 1, 1, 1, 1]
wknots = [0, 0, 1, 1]
knots = [uknots, vknots, wknots]
cx = [0 0 1]
nx = length(cx)
cy = [0 0 0 1]
ny = length(cy)
cz = [0 1]
nz = length(cz)
coefs = zeros(4, nx, ny, nz)
coefs[1,:,:,:] = repeat(reshape(cx,nx,1,1), 1, ny, nz);
coefs[2,:,:,:] = repeat(reshape(cy,1,ny,1), nx, 1, nz);
coefs[3,:,:,:] = repeat(reshape(cz,1,1,nz), nx, ny, 1);
coefs[4,:,:,:] .= 1
nurbs = nrbmak(coefs, knots)
x = rand(5,1)
y = rand(5,1)
z = rand(5,1)
tt = [x y z]'
points = nrbeval(nurbs,tt)
@test points ≈ [x.^2 y.^3 z]'

uknots = [0, 0, 0, 1, 1, 1]
vknots = [0, 0, 0, 0, 1, 1, 1, 1]
wknots = [0, 0, 1, 1]
knots = [uknots, vknots, wknots]
cx = [0 0 1]
nx = length(cx)
cy = [0 0 0 1]
ny = length(cy)
cz = [0 1]
nz = length(cz)
coefs = zeros(4, nx, ny, nz)
coefs[1,:,:,:] = repeat(reshape(cx,nx,1,1), 1, ny, nz)
coefs[2,:,:,:] = repeat(reshape(cy,1,ny,1), nx, 1, nz)
coefs[3,:,:,:] = repeat(reshape(cz,1,1,nz), nx, ny, 1)
coefs[4,:,:,:] .= 1
coefs = coefs[[2, 1, 3, 4],:,:,:]
nurbs = nrbmak(coefs, knots)
x = rand(5,1)
y = rand(5,1)
z = rand(5,1)
tt = [x y z]'
points = nrbeval(nurbs,tt)
[y.^3 x.^2 z]'
@test points ≈ [y.^3 x.^2 z]'

#test kntrefine()
knots = [[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
coefs = zeros(4, 2, 3)
coefs[1,:,:] = [1 sqrt(2)/2 0; 2 sqrt(2) 0]
coefs[2,:,:] = [0 sqrt(2)/2 1; 0 sqrt(2) 2]
coefs[4,:,:] = [1 sqrt(2)/2 1; 1 sqrt(2)/2 1]
nrbs = nrbmak(coefs, knots)
nrbs = nrbkntins(nrbs, [[], [0.5, 0.6, 0.6]])
nrbs = nrbdegelev(nrbs, [0, 1])
nrbs = nrbkntins(nrbs, [[], [0.4]])
rknots, _, _ = kntrefine(nrbs.knots, [1, 1], [1, 1], [0, 0])
@test (rknots[1] == [0, 0, 0.5, 1, 1])
@test (rknots[2] == [0, 0, 0.2, 0.4, 0.45, 0.5, 0.55, 0.6, 0.8, 1, 1])

rknots, _, _ = kntrefine(nrbs.knots, [1, 1], [3, 3], [0, 0])
@test (rknots[1] == [0, 0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1, 1])
@test (rknots[2] == [0, 0, 0, 0, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.45, 0.45, 0.45, 0.5,
                0.5, 0.5, 0.55, 0.55, 0.55, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 1, 1, 1, 1])

rknots, _, _ = kntrefine(nrbs.knots, [1, 1], [3, 3], [2, 2])
@test (rknots[1] == [0, 0, 0, 0, 0.5, 1, 1, 1, 1])
@test (rknots[2] == [0, 0, 0, 0, 0.2, 0.4, 0.45, 0.5, 0.5, 0.55, 0.6, 0.6, 0.6, 0.8, 1, 1, 1, 1])

#test nrbsquare()
srf = nrbsquare([], 1, 2, 2, 4)
@test (srf.order == [3, 3])
knt = [0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1]
@test (srf.knots == [knt, knt])
x = LinRange(0, 1, 100)
X = [i for i in x, j in x]
Y = [j for i in x, j in x]
vals = nrbeval(srf, [collect(x), collect(x)])
@test (vals[1,:,:] ≈ X)
@test (vals[2,:,:] ≈ 2*Y)
