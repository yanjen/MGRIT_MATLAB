function Y = PinST_relaxation(Y, h, dt, N, fc_ratio, level, num_processors)

ratio = fc_ratio ^ level;
h = h * ratio;
dt = dt * ratio;

L = length(Y)/N;
l = (L - 1) / ratio + 1;

n = floor((N - 1) / ratio);
npn = ceil(n / num_processors);

Z = Y;
last_proc_values = zeros(1,l);
for proc = 1:num_processors
    for i = (proc - 1) * npn + 1 : min(proc * npn, n)
        last_time_step = (i - 1) * ratio + 1;
        if i == (proc - 1) * npn + 1
            last_values = Z((last_time_step - 1) * L + 1 : ratio : last_time_step * L);
        else
            last_values = Y((last_time_step - 1) * L + 1 : ratio : last_time_step * L);
        end
        next_values = Lax_Wendroff(last_values, h, dt);
        next_values = next_values + last_proc_values;
        for j = i * ratio + 1 : min((i + 1) * ratio, N)
            Y((j - 1) * L + 1 : j * L) = interp1(0:h:300, next_values, 0:h/ratio:300);
        end
    end
    last_proc_values = next_values - Y((i * ratio + 1) * L + 1 : ratio : (i * ratio + 2) * L);
end

end