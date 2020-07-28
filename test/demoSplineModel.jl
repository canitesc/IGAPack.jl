#demos for SplineModel.jl module
include("../src/DemoGeometries.jl")
plotly()
#pyplot()
# ###### demo for nrbtestsrf()
# Plots.surface(reuse = false)
# srf = nrbtestsrf()
# plt = nrbplot(srf, [20, 30])
# display(plt)
#
# ##### demo for nrbtestcrv()
# Plots.surface(reuse = false)
# crv = nrbtestcrv()
# #plt = nrbplot(crv, [100])
# plt = nrbctrlplot(crv)
# display(plt)
#
# icrv = nrbkntins(crv, [[0.125, 0.375, 0.625, 0.875]])
# Plots.surface(reuse = false)
# plt =nrbctrlplot(icrv)
# display(plt)
#
# icrv = nrbdegelev(crv, [1])
# Plots.surface(reuse = false)
# plt =nrbctrlplot(icrv)
# display(plt)
#
#
# ###### demo for nrbtestcrv3d()
# Plots.surface(reuse = false)
# crv3d = nrbtestcrv3d()
# #plt = nrbplot(crv3d, [100])
# plt = nrbctrlplot(crv3d)
# display(plt)
#
# ##### demo for nrbtestvol()
# Plots.surface(reuse = false)
# horseshoe = nrbtestvol()
# plt = nrbplot(horseshoe, [6, 6, 50])
# display(plt)
#
# ##### demo for nrbline()
# Plots.surface(reuse = false)
# crv = nrbline([0.0 0.0 0.0]',[5.0 4.0 2.0]')
# #crv2 = nrbline()
# plt = nrbplot(crv, [1])
# #plt = nrbplot(crv2, [1])
# display(plt)
#
# ##### demo for nrbrect(w,h)
# Plots.surface(reuse = false)
# crv = nrbrect(2,1)
# plt = nrbplot(crv, [4])
# display(plt)
#
#
# ##### demo for nrbcirc()
# Plots.surface(reuse = false)
# for r=1:9
#     global crv = nrbcirc(r, [], 45*pi/180, 315*pi/180)
#     #global plt=nrbplot(crv, [50])
#     global plt=nrbkntplot(crv);
# end
# display(plt)
#
##### demo for nrbkntplot()
Plots.surface(reuse = false)
crv = nrbtestcrv()
plt = nrbkntplot(crv)
display(plt)

##### demo for nrbkntplot()
Plots.surface(reuse = false)
crv = nrbtestcrv3d()
plt = nrbkntplot(crv)
display(plt)

###### demo for nrbtestsrf()
Plots.surface(reuse = false)
srf = nrbtestsrf()
#plt = nrbkntplot(srf)
plt = nrbctrlplot(srf)
display(plt)
#
# ##### demo for nrbtestvol()
# Plots.surface(reuse = false)
# horseshoe = nrbtestvol()
# #plt = nrbkntplot(horseshoe)
# plt = nrbctrlplot(horseshoe)
# display(plt)
#
#
# ######### demo for nrbruled()
# Plots.surface(reuse = false)
# crv1 = nrbtestcrv()
# crv2 = nrbtform(nrbcirc(4,[4.5,0],pi,0.0), vectrans([0.0, 4.0, -4.0]));
# srf = nrbruled(crv1, crv2);
# plt = nrbplot(srf,[40 20]);
# #plt = nrbctrlplot(srf);
# display(plt)
# # title ('Ruled surface construction from two NURBS curves.');
# # hold off
#
# ######### demo for nrb4surf()
# Plots.surface(reuse = false)
# srf = nrb4surf([0.0, 0.0, 0.5],[1.0, 0.0, -0.5],[0.0, 1.0, -0.5],[1.0, 1.0, 0.5]);
# plt = nrbplot(srf,[10,10]);
# display(plt)
#
#
# ##### demo for nrbtransp()
# Plots.surface(reuse = false)
# srf = nrb4surf([0 0 0], [1 0 1], [0 1 1], [1 1 2]);
# plt = nrbplot(srf,[20 5]);
# coefs = srf.coefs
# coefs[3,:,:] = coefs[3,:,:] .+ 10
# srf = nrbmak(coefs, srf.knots);
# srf = nrbtransp(srf);
# plt = nrbplot(srf,[20 5]);
# display(plt)
#
# ##### demo for nrbcoons()
# Plots.surface(reuse = false)
# pnts = [ 0.0  3.0  4.5  6.5 8.0 10.0;
#          0.0  0.0  0.0  0.0 0.0  0.0;
#          2.0  2.0  7.0  4.0 7.0  9.0];
# crv1 = nrbmak(pnts, [[0, 0, 0, 1/3, 0.5, 2/3, 1, 1, 1]]);
#
# pnts= [ 0.0  3.0  5.0  8.0 10.0;
#         10.0 10.0 10.0 10.0 10.0;
#         3.0  5.0  8.0  6.0 10.0];
# crv2 = nrbmak(pnts, [[0, 0, 0, 1/3, 2/3, 1, 1, 1]]);
#
# pnts= [ 0.0 0.0 0.0 0.0;
#         0.0 3.0 8.0 10.0;
#         2.0 0.0 5.0 3.0];
# crv3 = nrbmak(pnts, [[0, 0, 0, 0.5, 1, 1, 1]]);
#
# pnts= [ 10.0 10.0 10.0 10.0 10.0;
#         0.0   3.0  5.0  8.0 10.0;
#         9.0   7.0  7.0 10.0 10.0];
# crv4 = nrbmak(pnts, [[0, 0, 0, 0.25, 0.75, 1, 1, 1]]);
# srf = nrbcoons(crv1, crv2, crv3, crv4);
#
# plt = nrbplot(srf,[20 20]);
# display(plt)

#### demo for nrbrevolve()
Plots.surface(reuse = false)
sphere = nrbrevolve(nrbcirc(1,[],0.0,pi),[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]);
plt = nrbplot(sphere,[40, 40]);
torus = nrbrevolve(nrbcirc(0.2,[0.9, 1.0]),[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]);
plt = nrbplot(torus,[40, 40]);
plt = nrbplot(nrbtform(torus,vectrans([-1.8])),[20, 10]);
display(plt)

#### demo for nrbrevolve()
Plots.surface(reuse = false)
pnts = [3.0 5.5 5.5 1.5 1.5 4.0 4.5;
        0.0 0.0 0.0 0.0 0.0 0.0 0.0;
        0.5 1.5 4.5 3.0 7.5 6.0 8.5]
crv = nrbmak(pnts,[[0, 0, 0, 1/4, 1/2, 3/4, 3/4, 1, 1, 1]])
xx = vecrotz(25*pi/180)*vecroty(15*pi/180)*vecrotx(20*pi/180)
nrb = nrbtform(crv,vectrans([5, 5])*xx)
pnt = [5, 5, 0]
vect = xx*[0, 0, 1, 1]
srf = nrbrevolve(nrb,pnt,vect[1:3])
plt = nrbplot(srf, [40, 40])
display(plt)

#### demo for nrbrevolve()
Plots.surface(reuse = false)
crv1 = nrbcirc(1,[0, 0],0, pi/2);
crv2 = nrbcirc(2,[0, 0],0, pi/2);
srf = nrbruled(crv1, crv2);
srf = nrbtform(srf, [1 0 0 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]);
vol = nrbrevolve(srf, [0, 0, 0], [1, 0, 0], pi/2);
plt = nrbplot(vol, [30, 30, 30])
display(plt)
