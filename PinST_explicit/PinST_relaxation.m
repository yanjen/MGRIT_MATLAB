function Y = PinST_relaxation(Y, h, dt, N, fc_ratio, level, total_level, nprocs)

coarse_ratio = fc_ratio ^ total_level;
fine_ratio = fc_ratio ^ level;
h_coarse = h * coarse_ratio;
h_fine = h * fine_ratio;
dt_coarse = dt * coarse_ratio;
dt_fine = dt * fine_ratio;

L = length(Y)/(nprocs + 1);

Z = Y;
%% Parallel Part
for proc = 1:nprocs
    n_start = floor(N * (proc - 1) / nprocs) + 1;
    n_start = floor((n_start - 1) / coarse_ratio) * coarse_ratio + 1;
    n_end = floor((N * proc) / nprocs) + 1;
    if proc ~= nprocs
        n_end = floor((n_end - 1) / coarse_ratio) * coarse_ratio + 1;
    end
    y_fine = interp1(0:h:300, Z((proc - 1) * L + 1:proc * L), 0:h_fine:300);
    for j = 1:(n_end - n_start)/fine_ratio
        y_fine = Lax_Wendroff(y_fine, h_fine, dt_fine);
    end
    y_coarse = interp1(0:h:300, Z((proc - 1) * L + 1:proc * L), 0:h_coarse:300);
    for j = 1:(n_end - n_start)/coarse_ratio
        y_coarse = Lax_Wendroff(y_coarse, h_coarse, dt_coarse);
    end
    Y(proc * L + 1:(proc + 1) * L) = interp1(0:h_fine:300, y_fine, 0:h:300) - interp1(0:h_coarse:300, y_coarse, 0:h:300);
end

%% Sequential Part
y = interp1(0:h:300, Y(1:L), 0:h_coarse:300);
for proc = 1:nprocs
    n_start = floor(N * (proc - 1) / nprocs) + 1;
    n_start = floor((n_start - 1) / coarse_ratio) * coarse_ratio + 1;
    n_end = floor((N * proc) / nprocs) + 1;
    if proc ~= nprocs
        n_end = floor((n_end - 1) / coarse_ratio) * coarse_ratio + 1;
    end
    for j = 1:(n_end - n_start)/coarse_ratio
        y = Lax_Wendroff(y, h_coarse, dt_coarse);
    end
    Y(proc * L + 1:(proc + 1) * L) = Y(proc * L + 1:(proc + 1) * L) + interp1(0:h_coarse:300, y, 0:h:300);
    y = interp1(0:h:300, Y(proc * L + 1:(proc + 1) * L), 0:h_coarse:300);
end

end