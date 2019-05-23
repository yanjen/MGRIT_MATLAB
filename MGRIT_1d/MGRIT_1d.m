function T_all = MGRIT_1d(T_init, Q, h, dt, N)

l = length(T_init);
e = ones((l-1)*N,1);
A = spdiags([-dt/h^2*e (1+2*dt/h^2)*e -dt/h^2*e], -1:1, (l-1), (l-1));
A(l-1,l-2) = A(l-1,l-2) - dt/h^2;
ACell = repmat({A}, 1, N);
A = blkdiag(ACell{:});
v = -ones(1,(l-1)*(N-1));
V = diag(v, -l+1);
A = A + sparse(V);

b = zeros(l-1,1);
b(l-1) = Q*dt;
b = repmat(b,N,1);

T_all = zeros((l-1)*N,1);

[A, ~, b] = reorganize_fine_coarse_1d(A,T_all,b,N,4);

T_all = (A\b)';

T_all = reverse_permute_fine_coarse_1d(T_all, N, 4);

end