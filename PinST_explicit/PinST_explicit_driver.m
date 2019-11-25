function Output = PinST_explicit_driver()
%% Parameters
h = 1.25;
dt = 1;
Max = 300;
x = 0:h:Max;
l = length(x);
N = 128;
Y = zeros(1, l);
for i = 1:length(x)
    if x(i) >= 50 && x(i) <= 110
        Y(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
fc_ratio = 2;
level = 3;
num_processors = 8;
tol = 0.2;
Y = repmat(Y,1,num_processors + 1);

figure;
Z = reshape(Y, l, num_processors + 1);
surf(Z)
pause(0.5);

%% Initialization
Y = PinST_initialization(Y, h, dt, N, num_processors, fc_ratio, level);
Z = reshape(Y, l, num_processors + 1);
surf(Z)
pause(0.5);

%% Parallel update by explicit PinST
last_time_step = Z(:, num_processors + 1);
iter_num = num_processors + 1;
for i = level - 1 : -1 : 0
    for j = 1:iter_num
        Y = PinST_relaxation(Y, h, dt, N, fc_ratio, i, level, num_processors);
        Z = reshape(Y, l, num_processors + 1);
        this_time_step = Z(:, num_processors + 1);
        error = sum(abs(last_time_step - this_time_step)) / l;
        disp(strcat("error = ", num2str(error)));
        if error < tol
            disp(strcat("In level ", int2str(i), ", iteration number is ", int2str(j)));
            break;
        end
        last_time_step = this_time_step;
        surf(Z)
        pause(0.5);
    end
end

%% Plot result comparison
% Plot exact result
figure;
A = zeros(1, l);
for i = 1:length(x)
    if x(i) >= 50+N && x(i) <= 110+N
        A(i) = 100*sin(pi*(x(i)-50-N)/60);
    end
end
plot(x, A, "LineWidth", 2);
hold on

% Plot Lax-Wendroff result
B = zeros(1, l);
for i = 1:length(x)
    if x(i) >= 50 && x(i) <= 110
        B(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
for i = 1:N
    B = Lax_Wendroff(B, h, dt);
end
plot(x, B, "LineWidth", 2);

% Plot explicit PinST result
plot(x, Z(:, num_processors + 1), "LineWidth", 2);

title(strcat("PinST result (tol = ", num2str(tol), ")"), "FontSize", 20);
legend("exact", "Lax-Wendroff", "explicit PinST", "Location", "northwest", "FontSize", 15);

%% Output last result
Output = Z(:,num_processors + 1)';

end