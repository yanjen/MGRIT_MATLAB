%% Plot 2D compressible in MATLAB
% This colormap comes from paraview
b = [0.231373 0.298039 0.752941];
w = [0.865003 0.865003 0.865003];
r = [0.705882 0.0156863 0.14902];
Nch = 64;

cmap = zeros(2*Nch+1, 3);
for i = 1:3
    cmap(1:Nch+1, i) = linspace(b(i), w(i), Nch+1);
    cmap(Nch+1:end, i) = linspace(w(i), r(i), Nch+1);
end

N = 64;
inner_radius = 5;
outer_radius = 100;
[x, y, r, dtheta] = gen_donut_mesh2(inner_radius, outer_radius, N);

fcfd = figure;
hold on;
vertices = [x' y'];
faces = zeros(N*(length(r) - 1), 4);
for i = 1:length(r) - 1
    for j = 1:N
        if j == N
            faces((i-1)*N+j,1) = (i-1)*N+j;
            faces((i-1)*N+j,2) = i*N+j;
            faces((i-1)*N+j,3) = (i-1)*N+j+1;
            faces((i-1)*N+j,4) = (i-2)*N+j+1;
        else
            faces((i-1)*N+j,1) = (i-1)*N+j;
            faces((i-1)*N+j,2) = i*N+j;
            faces((i-1)*N+j,3) = i*N+j+1;
            faces((i-1)*N+j,4) = (i-1)*N+j+1;
        end
    end
end
data = [zeros(N,1);ones(N*(length(r) - 1),1)];
patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', data, 'FaceColor', 'interp');
colormap(cmap);
colorbar
axis equal;
