function compressible_1d_driver()
%% Parameters
N = 65;
dt = 0.01;
x = linspace(1,15,N);
dx = x(2) - x(1);
density = 1;
x_velocity = 0.1;
gamma = 1.4;
P = 1/gamma * ones(1,N);
max_iter = 600;

%% Initialization
rho = density * ones(1,N);
rhou = [0, x_velocity * density * ones(1,N-1)];
E = P ./ (gamma - 1) .* ones(1,N) + 0.5 .* rhou .* rhou ./ rho;
u = rhou ./ rho;

%% Plotting
subplot(3,1,1)
plot(x,rho)
ylabel("Density")
subplot(3,1,2)
plot(x,u)
ylabel("Velocity")
subplot(3,1,3)
plot(x,E)
ylabel("Energy")

%% Time marching
for iter = 1:max_iter
    [rho, rhou, E, P] = Lax_Wendroff_vecv(rho, rhou, E, P, gamma, dt);
    
    [rho, rhou, E] = extrapolate_boundary(rho, rhou, E);
    
    [rhou, E] = viscosity_fix(rho, rhou, E, P, gamma, dx, dt);
    
    u = rhou ./ rho;
    
    %% Plotting
    subplot(3,1,1)
    plot(x,rho)
    ylabel("Density")
    subplot(3,1,2)
    plot(x,u)
    ylabel("Velocity")
    subplot(3,1,3)
    plot(x,E)
    ylabel("Energy")
    
    pause(0.01)
end

end