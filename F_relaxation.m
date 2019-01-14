function x = F_relaxation(A, S, x, b)

n = size(S, 1);

x = (eye(n)-S*((S'*A*S)\S')*A)*x + b;

end