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
    if x(i) > 49 && x(i) < 111
        Y(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
fc_ratio = 2;
level = 3;
num_processors = 8;
Y = repmat(Y,1,num_processors + 1);

Z = reshape(Y, l, num_processors + 1);
surf(Z)
pause(0.5);

%% Initialization
Y = PinST_initialization(Y, h, dt, N, num_processors, fc_ratio, level);
Z = reshape(Y, l, num_processors + 1);
surf(Z)
pause(0.5);

iter_num = 1;
for i = level - 1 : -1 : 0
    for j = 1:iter_num
        Y = PinST_relaxation(Y, h, dt, N, fc_ratio, i, level, num_processors);
        Z = reshape(Y, l, num_processors + 1);
        surf(Z)
        pause(0.5);
    end
end

Output = Z(:,num_processors + 1)';

end