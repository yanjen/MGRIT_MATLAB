function x = FCF_relaxation(R, A, P, x, b)

n = size(R, 2);
nc = size(R, 1);
nf = n - nc;

I = eye(nc);
R1 = [zeros(nc, nf) I];

r = b - A*x;
x = A\(P*(I-R*A*P)*R1*r+b);

end