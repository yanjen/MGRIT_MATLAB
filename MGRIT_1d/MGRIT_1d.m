function x = MGRIT_1d(A, x, b, N, fc_ratio)

l = length(b) / N;
n = l * N;
nc = l * ceil(N/fc_ratio);
nf = n - nc;

if length(A) < 5*l
    x = A\b;
else

A_ff = A(1:nf, 1:nf);
A_fc = A(1:nf, nf+1:n);
A_cf = A(nf+1:n, 1:nf);

% Restriction Matrix
R = [-A_cf*(A_ff\speye(nf)) speye(nc)];
% Prolongation Matrix
P = [-A_ff\A_fc;speye(nc)];

R1 = [zeros(nc, nf) speye(nc)];

%% Do FCF relaxization
x = FCF_relaxation(R, A, P, x, b);

%% Restriction
b2 = R1*(b - A*x);

%% Coarse grid operation
x2 = MGRIT_0d(R*A*P, b2, R*x, fc_ratio);

%% Correction
x = x + P*x2;

end

end