function v = reverse_permute_fine_coarse_1d(v, N, ratio)

Nc = ceil(N/ratio);
Nf = N - Nc;
l = length(v)/N;
u = zeros(l*N, 1);

for i = 1:Nc
    u((i-1)*ratio*l+1:(1+(i-1)*ratio)*l) = v((Nf+i-1)*l+1:(Nf+i)*l);
    if N > 1 + (i-1)*ratio
        u_end = min(l*N, i*ratio*l);
        v_end = (i-1)*(ratio-1)*l+1 + u_end - ((2+(i-1)*ratio-1)*l+1);
        u((2+(i-1)*ratio-1)*l+1:u_end) = v((i-1)*(ratio-1)*l+1:v_end);
    end
end

v = u(1:l*N);

end