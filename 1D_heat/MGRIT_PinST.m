function T = MGRIT_PinST(T, Q, h, dt, N, iter)

% Define fine-coarse ratio
fc_ratio = 4;

%% MGRIT Parallel-in-Space/Time method
[A, x, b] = gen_1d_heat_transfer(T, Q, h, dt, N);

for i = 1:iter
    x = MGRIT_1d(A, x, b, N, fc_ratio);
end

l = length(x)/N + 1;
T = zeros(1,l*N);
for i = 1:N
    T((i-1)*l + 2:i*l) = x((i-1)*(l-1) + 1:i*(l-1));
end

end