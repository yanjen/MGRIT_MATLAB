function PinST_explicit_driver()
%% Parameters
h = 1.25;
dt = 1;
Max = 300;
x = 0:h:Max;
l = length(x);
N = 60;
Y = zeros(1, l);
for i = 1:length(x)
    if x(i) > 49 && x(i) < 111
        Y(i) = 100*sin(pi*(x(i)-50)/60);
    end
end
Y = repmat(Y,1,N);
fc_ratio = 2;
level = 4;
num_processors = 4;

Z = reshape(Y, l, N);
surf(Z)
pause(1);

%% Initialization
Y = PinST_initialization(Y, h, dt, N, fc_ratio, level);
Z = reshape(Y, l, N);
surf(Z)
pause(1);

for i = level - 1 : -1 : 1
    Y = PinST_relaxation(Y, h, dt, N, fc_ratio, i, num_processors);
    Z = reshape(Y, l, N);
    surf(Z)
    pause(1);
end

for j = 1:2
    Y = PinST_relaxation(Y, h, dt, N, fc_ratio, 1, num_processors);
    Z = reshape(Y, l, N);
    surf(Z)
    pause(1)
end

end