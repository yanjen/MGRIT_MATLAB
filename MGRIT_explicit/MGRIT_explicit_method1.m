function T = MGRIT_explicit_method1(T, Q, h, dt, N, fc_ratio, level)

if level == 0
    return
end

l = length(T) / N;
R = T;
for i = 1:N-1
    R(i*l + 1:(i+1)*l) = Lax_Wendroff(T((i-1)*l + 1:i*l), h, dt);
end

N2 = ceil(N / fc_ratio);
l2 = ceil(l / fc_ratio);
T2 = zeros(1, l2*N2);
for i = 1:N2
    for j = 1:l2
        T2((i-1)*l2 + (j-1) + 1) = R((i-1)*fc_ratio*l + fc_ratio*(j-1) + 1);
    end
end
R2 = MGRIT_explicit_method1(T2, Q, h*fc_ratio, dt*fc_ratio, N2, fc_ratio, level - 1);
for i = 1:N2
    for j = 1:l2
        for k = fc_ratio:-1:1
            if ((i-1)*fc_ratio)*l + fc_ratio*(j-1) + k < length(T)+1
                R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + k) = R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + k) + R2((i-1)*l2 + (j-1) + 1) - R((i-1)*fc_ratio*l + fc_ratio*(j-1) + 1);
            end
        end
    end
end

S = R;
for i = 1:N-1
    S(i*l + 1:(i+1)*l) = Lax_Wendroff(R((i-1)*l + 1:i*l), h, dt);
end

T = S;

end