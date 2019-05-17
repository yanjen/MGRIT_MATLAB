function v = reverse_permute_fine_coarse(v, ratio)

n = length(v);
nc = ceil(n/ratio);
nf = n - nc;
u = zeros(n, 1);

for i = 1:nc
    u(1+(i-1)*ratio) = v(nf+i);
    if n > 1 + (i-1)*ratio
        u(2+(i-1)*ratio:i*ratio) = v(1+(i-1)*(ratio-1):i*(ratio-1));
    end
end

v = u(1:n);
end