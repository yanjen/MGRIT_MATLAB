function compressible_1d_driver()
%% Parameters
N = 32;
dt = 0.01;
x = linspace(1,15,N);
density = 1;
x_velocity = 0.1;
gamma = 1.4;
P = 1/gamma * ones(1,N);
max_iter = 600;

%% Initialization
rho = density * ones(1,N);
rhou = [0, x_velocity * density * ones(1,N-1)];
E = P ./ (gamma - 1) .* ones(1,N) + 0.5 .* rhou .* rhou ./ rho;

%% Plotting
subplot(3,1,1)
plot(x,rho)
subplot(3,1,2)
plot(x,rhou)
subplot(3,1,3)
plot(x,E)

%% Time marching
for iter = 1:max_iter
    u = rhou ./ rho;
    F1 = rhou;
    F2 = rhou .* u + P;
    F3 = (E + P) .* u;
    P = (E - 0.5 .* rho .* u .* u) .* (gamma - 1);
    
    rho = Lax_Wendroff_vecv(rho, F1, dt, 2);
    rhou = Lax_Wendroff_vecv(rhou, F2, dt, 0);
    E = Lax_Wendroff_vecv(E, F3, dt, 2);
    
    %% Plotting
    subplot(3,1,1)
    plot(x,rho)
    subplot(3,1,2)
    plot(x,rhou)
    subplot(3,1,3)
    plot(x,E)
    
    pause(0.05)
end

end