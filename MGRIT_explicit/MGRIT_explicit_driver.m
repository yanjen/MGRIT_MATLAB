function MGRIT_explicit_driver()

Q = 1;
h = 0.1;
dt = 0.001;
x = 0:h:1;
l = length(x);
N = 40;
T = zeros(1, l*N);
fc_ratio = 2;
level = 2;
iter = 20;

figure
subplot(2,1,1)
title('Result at time T')
hold on
subplot(2,1,2)
title('Results for all time steps')
hold on

%% Plot MGRIT explicit result by iteration
for i = 1:iter
    T = MGRIT_explicit_method1(T, Q, h, dt, N, fc_ratio, level);
    subplot(2,1,1)
    plot(x, T(l*(N-1) + 1:l*N));
    subplot(2,1,2)
    plot(linspace(0, N, l*N), T);
    pause(0.5)
end

%% Plot explicit time marching result
T_ETM = zeros(size(x));
for i = 1:N-1
    T_ETM = explicit_time_marching(T_ETM, Q, h, dt);
end
subplot(2,1,1)
p_ETM = plot(x, T_ETM); M1 = 'Explicit time marching';
pause(0.5)

%% Plot with 2x coarse x-grid
h = 0.2;
x = 0:h:1;
T_2 = zeros(size(x));
for i = 1:N-1
    T_2 = explicit_time_marching(T_2, Q, h, dt);
end
p_2 = plot(x, T_2); M2 = 'Explicit with 2x coarser x-grids';
legend([p_ETM; p_2], M1, M2);

end