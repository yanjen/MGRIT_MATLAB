function x = MGRIT_1d(A, x, b, N, fc_ratio)

l = length(b) / N;
n = l * N;
nc = l * ceil(N/fc_ratio);
nf = n - nc;

if nf == 0
    x = A\b;
else

%% Reorganize the matrix and 
[B, x, y] = reorganize_fine_coarse_1d(A, x, b, N, fc_ratio);

A_ff = B(1:nf, 1:nf);
A_fc = B(1:nf, nf+1:n);
A_cf = B(nf+1:n, 1:nf);

% Restriction Matrix
R = [-A_cf*(A_ff\speye(nf)) speye(nc)];
% Prolongation Matrix
P = [-A_ff\A_fc;speye(nc)];

R1 = [zeros(nc, nf) speye(nc)];

P1 = [zeros(nf, nc);speye(nc)];

%% Do FCF relaxization
x = FCF_relaxation(R, B, P, x, y);

%% Restriction
b2 = R1*(y - B*x);

%% Coarse grid operation
x2 = MGRIT_1d(R1*B*P1, b2, R1*x, nc/l, fc_ratio);

%% Correction
x = x + P*x2;

%% Reverse the permutation
x = reverse_permute_fine_coarse_1d(x, N, fc_ratio);

end

end