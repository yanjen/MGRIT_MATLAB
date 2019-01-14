function [A, x, b] = gen_trivial(n)

A = eye(n);
for i = 2:n
A(i,i-1) = -1;
end
A = sparse(A);

b = ones(n, 1);
x = randn(n, 1);
x(1) = 1;

end