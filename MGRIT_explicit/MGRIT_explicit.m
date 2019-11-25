function T = MGRIT_explicit(T, Q, h, dt, N, fc_ratio, level)

if level == 0
    return
end

l = length(T) / N;
R = T;
for i = 1:N-1
    R(i*l + 1:(i+1)*l) = Lax_Wendroff(T((i-1)*l + 1:i*l), h, dt);
end

N2 = ceil(N / fc_ratio);
T2 = zeros(1, l*N2);
for i = 1:N2
    T2((i-1)*l + 1:i*l) = R((i-1)*fc_ratio*l + 1:((i-1)*fc_ratio + 1)*l);
end
R2 = MGRIT_explicit(T2, Q, h, dt*fc_ratio, N2, fc_ratio, level - 1);
for i = 1:N2
    R((i-1)*fc_ratio*l + 1:((i-1)*fc_ratio + 1)*l) = R2((i-1)*l + 1:i*l);
end

S = R;
for i = 1:N-1
    S(i*l + 1:(i+1)*l) = Lax_Wendroff(R((i-1)*l + 1:i*l), h, dt);
end

T = S;

end