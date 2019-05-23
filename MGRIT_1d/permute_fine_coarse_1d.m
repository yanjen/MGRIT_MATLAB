function v = permute_fine_coarse_1d(v, N, ratio)

Nc = ceil(N/ratio);
l = length(v)/N;

vc = zeros(l*Nc,1);
for i = 1:Nc
    vc(1 + (i-1)*l:i*l) = v((i-1)*(ratio-1)*l + 1:(1+(i-1)*(ratio-1))*l);
    v((i-1)*(ratio-1)*l + 1:(1+(i-1)*(ratio-1))*l) = [];
end
v = [v;vc];

end