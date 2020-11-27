%% Plot vector graph for 2D CFD
rescale = 5;

N = 64;
inner_radius = 5;
outer_radius = 100;
[x, y, r, dtheta] = gen_donut_mesh2(inner_radius, outer_radius, N);

load("/Users/chenjenchen/Mirror/cluster/pinst_compressible/build/example/U2.txt")
load("/Users/chenjenchen/Mirror/cluster/pinst_compressible/build/example/U3.txt")

U2 = U2.*rescale;
U3 = U3.*rescale;

figure
quiver(x,y,U2',U3','AutoScale','off');