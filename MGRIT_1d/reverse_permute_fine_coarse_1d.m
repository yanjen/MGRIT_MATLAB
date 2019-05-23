function v = reverse_permute_fine_coarse_1d(v, N, ratio)

Nc = ceil(N/ratio);
Nf = N - Nc;
l = length(v)/N;
u = zeros(l*N, 1);

for i = 1:Nc
    u((i-1)*ratio*l+1:(1+(i-1)*ratio)*l) = v((Nf+i-1)*l+1:(Nf+i)*l);
    if N > 1 + (i-1)*ratio
        u((2+(i-1)*ratio-1)*l+1:i*ratio*l) = v((i-1)*(ratio-1)*l+1:i*(ratio-1)*l);
    end
end

v = u(1:l*N);

end