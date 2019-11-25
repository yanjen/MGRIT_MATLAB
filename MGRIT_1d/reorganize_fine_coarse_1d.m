function [A, x, b] = reorganize_fine_coarse_1d(A, x, b, N, ratio)

Nc = ceil(N/ratio);
Nf = N - Nc;

x = permute_fine_coarse_1d(x, N, ratio);
b = permute_fine_coarse_1d(b, N, ratio);

l = length(b)/N;
B = spdiags([[diag(A,-1);0] diag(A) [0;diag(A,1)]], -1:1, l*N, l*N);
for i = 2:N
    if (mod(i-1, ratio) == 1)
        B((i - ceil(i/ratio) - 1)*l + 1:(i - ceil(i/ratio))*l, (Nf+ceil((i-1)/ratio) - 1)*l + 1:(Nf+ceil((i-1)/ratio))*l) = A((i-1)*l+1:i*l, (i-2)*l+1:(i-1)*l);
        continue;
    end
    if (mod(i, ratio) == 1)
        B((Nf+ceil(i/ratio)-1)*l+1:(Nf+ceil(i/ratio))*l, (i-ceil(i/ratio)-1)*l+1:(i-ceil(i/ratio))*l) = A((i-1)*l+1:i*l, (i-2)*l+1:(i-1)*l);
        continue;
    end
    B((i-ceil(i/ratio)-1)*l+1:(i-ceil(i/ratio))*l, (i-1-ceil(i/ratio)-1)*l+1:(i-1-ceil(i/ratio))*l) = A((i-1)*l+1:i*l, (i-2)*l+1:(i-1)*l);
end
A = B;

end