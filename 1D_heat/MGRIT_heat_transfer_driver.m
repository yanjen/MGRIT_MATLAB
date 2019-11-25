function MGRIT_heat_transfer_driver()

Q = 1;
h = 0.2;
dt = 0.001;
x = 0:h:1;
N = 40;
T = zeros(1, length(x)*N);
fc_ratio = 2;
level = 2;
iter = 20;

figure
% subplot(2,1,1)
hold on
% subplot(2,1,2)
% hold on

for i = 1:iter
    T = MGRIT_heat_explicit(T, Q, h, dt, N, fc_ratio, level);
%     subplot(2,1,1)
    plot(x, T(N*6-5:N*6))
%     subplot(2,1,2)
%     plot(1:N*6, T)
    pause(0.1)
end

T_ETM = zeros(1, length(x));
for i = 1:N-1
    T_ETM = explicit_time_marching(T_ETM, Q, h, dt);
end
% subplot(2,1,1)
plot(x, T_ETM)

end