function v = permute_fine_coarse(v, ratio)

n = length(v);
nc = ceil(n/ratio);

vc = zeros(nc, 1);
for i = 1:nc
    vc(i) = v(1+(i-1)*(ratio-1));
    v(1+(i-1)*(ratio-1)) = [];
end
v = [v;vc];
end