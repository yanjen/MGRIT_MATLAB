function MGRIT_explicit_driver()
%% Initialization
option = 'NoCoarseDirect';
h = 1.25;
dt = 1;
Max = 1000;
x = 0:h:Max;
l = length(x);
N = 600;
Y = zeros(1, l);
for i = 1:length(x)
    if x(i) > 49 && x(i) < 111
        Y(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
Y = repmat(Y,1,N);
fc_ratio = 4;
level = 2;
iter = 43;

%% Initial different y values ====
% ratio = fc_ratio^(level - 1);
% for i = 2:N
%     if mod(i,fc_ratio) == 0
%         Y((i-1)*l + 1:i*l) = Lax_Wendroff(Y((i-2)*l + 1:(i-1)*l), h*ratio, dt*ratio);
%     else
%         Y((i-1)*l + 1:i*l) = Y((i-2)*l + 1:(i-1)*l);
%     end
% end
%=================================

f = figure;
% subplot(2,1,1)
title('Explicit MGRIT','FontSize',20)
hold on
% subplot(2,1,2)
% title('Results for all time steps')
% hold on
% axis tight manual

%% Plot theoretical answer
Z = zeros(1, l);
for i = 1:length(x)
    if x(i) > 49 + N && x(i) < 111 + N - 1
        Z(i) = 100*sin(pi*(x(i)-(50 + N-1))/60);
    end
end
% subplot(2,1,1)
p_Z = plot(x, Z, 'LineWidth', 1.5); M0 = 'Exact solution';
ylim([-20 120])
ylabel('value','FontSize',20)
xlabel('x','FontSize',20)
% title('Explicit MGRIT without coarest explicit and with fine update','FontSize',15)

%% Plot MGRIT explicit result by iteration
for i = 1:iter
    Y = MGRIT_explicit_method1(Y, h, dt, N, fc_ratio, level, option);
%     subplot(2,1,1)
    plot(x, Y(l*(N-1) + 1:l*N));
%     subplot(2,1,2)
%     plot(linspace(0, N, l*N), Y);
    pause(0.05)
end
plot(x, Y(l*(N-1) + 1:l*N), 'LineWidth', 2);

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
% subplot(2,1,1)
p_ETM = plot(x, Y_ETM,'LineWidth',2); M1 = 'Lax-Wendroff';
legend([p_Z,p_ETM], M0, M1,'FontSize',15);
pause(0.5)

%% Plot implicit time marching result
% Y_ITM = zeros(1, l);
% for i = 1:length(x)
%     if x(i) > 49 && x(i) < 111
%         Y_ITM(i) = 100*sin(pi*(x(i)-50)/60);
%     end
% end
% for i = 1:N-1
%     Y_ITM = backward_euler(Y_ITM, h, dt);
% end
% subplot(2,1,1)
% p_ITM = plot(x, Y_ITM,'LineWidth',2); M2 = 'Backward Euler';
% legend([p_Z,p_ITM], M0, M2,'FontSize',15);
% pause(0.5)

%% Plot with coarse grids
ratio = fc_ratio^(level - 1);
h = h * ratio;
dt = dt * ratio;
x = 0:h:Max;
Y_2 = zeros(size(x));
for i = 1:length(x)
    if x(i) > 49 && x(i) < 111
        Y_2(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
for i = 1:N/ratio - 1
    Y_2 = Lax_Wendroff(Y_2, h, dt);
end
Y_2 = Lax_Wendroff(Y_2, h, dt/ratio);
Y_2 = Lax_Wendroff(Y_2, h, dt/ratio);
Y_2 = Lax_Wendroff(Y_2, h, dt/ratio);
% subplot(2,1,1)
p_2 = plot(x, Y_2, '-x'); M3 = 'Explicit with coarsest grids';
legend([p_ETM; p_2], M1, M3);

end