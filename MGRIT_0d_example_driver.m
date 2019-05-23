%% MGRIT driver

fc_ratio = 4;
[A, x, b] = gen_trivial(10);

[A, x, b] = reorganize_fine_coarse(A, x, b, fc_ratio);

x = MGRIT(A, b, x, fc_ratio);

x = reverse_permute_fine_coarse(x, fc_ratio);