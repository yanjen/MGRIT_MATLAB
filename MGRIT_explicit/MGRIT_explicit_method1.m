function T = MGRIT_explicit_method1(T, h, dt, N, fc_ratio, level, option)
l = length(T) / N;


if level == 0
    if strcmp(option,'CoarseDirect') == 1
        for i = 1:N-1
            T(i*l + 1:(i+1)*l) = Lax_Wendroff(T((i-1)*l + 1:i*l), h, dt);
        end
    end
    return
end

R = T;
for i = 1:N-1
    R(i*l + 1:(i+1)*l) = Lax_Wendroff(T((i-1)*l + 1:i*l), h, dt);
end
T = R;

N2 = ceil(N / fc_ratio);
l2 = ceil(l / fc_ratio);
T2 = zeros(1, l2*N2);
for i = 1:N2
    for j = 1:l2
        T2((i-1)*l2 + (j-1) + 1) = R((i-1)*fc_ratio*l + fc_ratio*(j-1) + 1);
    end
end
R2 = MGRIT_explicit_method1(T2, h*fc_ratio, dt*fc_ratio, N2, fc_ratio, level - 1, option);
for i = 1:N2
    for j = 1:l2
%% Update same amount
%         for k = fc_ratio:-1:1
%             if ((i-1)*fc_ratio)*l + fc_ratio*(j-1) + k < length(T)+1
%                 R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + k) = R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + k) + R2((i-1)*l2 + (j-1) + 1) - R((i-1)*fc_ratio*l + fc_ratio*(j-1) + 1);
%             end
%         end
%% Update only on coarse
%         R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + 1) = R2((i-1)*l2 + (j-1) + 1);
%% Update the average of front and back update
        R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + 1) = R2((i-1)*l2 + (j-1) + 1);
    end
    for j = 1:l2-1
            front = R2((i-1)*l2 + (j-1) + 1) - T((i-1)*fc_ratio*l + fc_ratio*(j-1) + 1);
%             back = R2((i-1)*l2 + (j-1) + 2) - T((i-1)*fc_ratio*l + fc_ratio*(j) + 1);
            R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + 2) = R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + 2) + front;
    end
%% Update by interpolation
%     for j = 1:l2 - 1
%         R(((i-1)*fc_ratio)*l + fc_ratio*(j-1) + 2) = (R2((i-1)*l2 + (j-1) + 1) + R2((i-1)*l2 + j + 1))/2;
%     end
end

S = R;
for i = 1:N-1
    if mod(i,2) == 1
        S(i*l + 1:(i+1)*l) = Lax_Wendroff(R((i-1)*l + 1:i*l), h, dt);
    end
end

T = S;

end