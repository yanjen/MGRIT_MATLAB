function x = MGRIT_0d(A, b, x, fc_ratio)

n = length(b);
nc = ceil(n/fc_ratio);
nf = n - nc;

if length(A) < 5
    x = A\b;
else
    
[A, x, b] = reorganize_fine_coarse(A, x, b, fc_ratio);

A_ff = A(1:nf, 1:nf);
A_fc = A(1:nf, nf+1:n);
A_cf = A(nf+1:n, 1:nf);

% Restriction Matrix
R = [-A_cf*(A_ff\eye(nf)) eye(nc)];
% Prolongation Matrix
P = [-A_ff\A_fc;eye(nc)];

R1 = [zeros(nc, nf) eye(nc)];

P1 = [zeros(nf, nc);eye(nc)];

%% Do FCF relaxization
x = FCF_relaxation(R, A, P, x, b);

%% Restriction
b2 = R1*(b - A*x);

%% Coarse grid operation
x2 = MGRIT_0d(R1*A*P1, b2, R1*x, fc_ratio);

%% Correction
x = x + P*x2;

x = reverse_permute_fine_coarse(x, fc_ratio);

end

end