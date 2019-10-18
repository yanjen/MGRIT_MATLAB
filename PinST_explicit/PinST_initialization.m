function Y = PinST_initialization(Y, h, dt, N, fc_ratio, level)

ratio = fc_ratio^level;
h = h * ratio;
dt = dt * ratio;

L = length(Y)/N;
l = (L - 1)/ratio + 1;

n = floor((N - 1)/ratio);
for i = 1:n
    last_time_step = (i - 1) * ratio + 1;
    last_values = Y((last_time_step - 1) * L + 1 : ratio : last_time_step * L);
    next_values = Lax_Wendroff(last_values, h, dt);
    for j = i * ratio + 1 : min((i + 1) * ratio, N)
        Y((j - 1) * L + 1 : j * L) = interp1(0:h:300, next_values, 0:h/ratio:300);
    end
end

end