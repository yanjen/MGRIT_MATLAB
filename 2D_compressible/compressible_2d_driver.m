function compressible_2d_driver()
%% Parameters
N = 32;
dt = 0.01;
inner_radius = 1;
outer_radius = 15;
[x, y, dr, dtheta] = gen_donut_mesh(inner_radius, outer_radius, N + 1, N);
l = length(x);
density = 1;
x_velocity = 0.1;
y_velocity = 0;
gamma = 1.4;
P = 1/gamma * ones(1,l);
max_iter = 300;

%% Initialization
rho = density * ones(1,l);
rhou = [zeros(1,N), x_velocity * density * ones(1,l-N)];
rhov = y_velocity * density * ones(1,l);
E = P ./ (gamma - 1) .* ones(1,l) + 0.5 .* (rhou .* rhou + rhov .* rhov) ./ rho;
u = rhou ./ rho;
v = rhov ./ rho;

%% Plotting
% subplot(2,2,1)
% scatter3(x,y,rho)
% zlabel("Density")
% subplot(2,2,2)
% scatter3(x,y,u)
% zlabel("x-Velocity")
% subplot(2,2,3)
% scatter3(x,y,v)
% zlabel("y-Velocity")
% subplot(2,2,4)
% scatter3(x,y,E)
% zlabel("Energy")
xx = x(abs(y) < 1e-12);
uu = u(abs(y) < 1e-12);
[xx, order] = sort(xx);
uu = uu(order);
plot(xx,uu)

%% Time marching
for iter = 1:max_iter
    [rho, rhou, rhov, E, P] = Lax_Wendroff_vecv3d(rho, rhou, rhov, E, P, x, y, N, gamma, dr, dtheta, dt);
    
    [rho, rhou, rhov, E] = extrapolate_boundary3d(rho, rhou, rhov, E, P, N, gamma);
    
%     [rhou, E] = viscosity_fix(rho, rhou, E, P, gamma, dx, dt);
    
    u = rhou ./ rho;
    v = rhov ./ rho;
    
    %% Plotting
%     subplot(2,2,1)
%     scatter3(x,y,rho)
%     zlabel("Density")
%     subplot(2,2,2)
%     scatter3(x,y,u)
%     zlabel("x-Velocity")
%     subplot(2,2,3)
%     scatter3(x,y,v)
%     zlabel("y-Velocity")
%     subplot(2,2,4)
%     scatter3(x,y,E)
%     zlabel("Energy")
    xx = x(abs(y) < 1e-12);
    uu = u(abs(y) < 1e-12);
    [xx, order] = sort(xx);
    uu = uu(order);
    plot(xx,uu)
    
    pause(0.1)
end

end