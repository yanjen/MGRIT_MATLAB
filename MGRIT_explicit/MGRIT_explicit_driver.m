function MGRIT_explicit_driver()

Q = 1;
h = 2.5;
dt = 0.0075;
x = 0:h:300;
l = length(x);
N = 40;
Y = zeros(1, l);
for i = 1:length(x)
    if x(i) > 49 && x(i) < 111
        Y(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
Y = repmat(Y,1,N);
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
% axis tight manual

%% Plot MGRIT explicit result by iteration
for i = 1:iter
    Y = MGRIT_explicit_method1(Y, Q, h, dt, N, fc_ratio, level);
    subplot(2,1,1)
    plot(x, Y(l*(N-1) + 1:l*N));
    subplot(2,1,2)
    plot(linspace(0, N, l*N), Y);
    pause(0.5)
end

%% Plot explicit time marching result
Y_ETM = zeros(1, l);
for i = 1:length(x)
    if x(i) > 49 && x(i) < 111
        Y_ETM(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
for i = 1:N-1
    Y_ETM = Lax_Wendroff(Y_ETM, h, dt);
end
subplot(2,1,1)
p_ETM = plot(x, Y_ETM); M1 = 'Explicit time marching';
legend(p_ETM, M1);
pause(0.5)

% %% Plot with 2x coarse x-grid
% h = 0.2;
% x = 0:h:1;
% T_2 = zeros(size(x));
% for i = 1:N-1
%     T_2 = explicit_time_marching(T_2, Q, h, dt);
% end
% p_2 = plot(x, T_2); M2 = 'Explicit with 2x coarser x-grids';
% legend([p_ETM; p_2], M1, M2);

end