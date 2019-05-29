function [A, x, b] = gen_1d_heat_transfer(T, Q, h, dt, N)

l = length(T);
e = ones((l-1)*N,1);
A = spdiags([-dt/h^2*e (1+2*dt/h^2)*e -dt/h^2*e], -1:1, (l-1), (l-1));
A(l-1,l-2) = A(l-1,l-2) - dt/h^2;
ACell = repmat({A}, 1, N);
A = blkdiag(ACell{:});
v = -ones((l-1)*(N-1), 1);
V = spdiags(v, -l+1, (l-1)*N, (l-1)*(N-1));
V = [V zeros((l-1)*N, l-1)];
A = A + sparse(V);

b = zeros(l-1,1);
b(l-1) = Q*dt;
b = repmat(b,N,1);

x = zeros((l-1)*N,1);

end