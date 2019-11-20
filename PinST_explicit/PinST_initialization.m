function Y = PinST_initialization(Y, h, dt, N, nprocs, fc_ratio, level)

ratio = fc_ratio^level;
h_coarse = h * ratio;
dt = dt * ratio;

L = length(Y)/(nprocs + 1);
y = interp1(0:h:300, Y(1:L), 0:h_coarse:300);

for proc = 1:nprocs
    n_start = floor(N * (proc - 1) / nprocs) + 1;
    n_start = floor((n_start - 1) / ratio) * ratio + 1;
    n_end = floor((N * proc) / nprocs) + 1;
    if proc ~= nprocs
        n_end = floor((n_end - 1) / ratio) * ratio + 1;
    end
    for j = 1:(n_end - n_start)/ratio
        y = Lax_Wendroff(y, h_coarse, dt);
    end
    Y(proc*L + 1:(proc + 1)*L) = interp1(0:h_coarse:300, y, 0:h:300);
end
% L = length(Y)/N;
% l = (L - 1)/ratio + 1;
% 
% n = floor((N - 1)/ratio);
% for i = 1:n
%     last_time_step = (i - 1) * ratio + 1;
%     last_values = Y((last_time_step - 1) * L + 1 : ratio : last_time_step * L);
%     next_values = Lax_Wendroff(last_values, h, dt);
%     for j = i * ratio + 1 : min((i + 1) * ratio, N)
%         Y((j - 1) * L + 1 : j * L) = interp1(0:h:300, next_values, 0:h/ratio:300);
%     end
% end

end