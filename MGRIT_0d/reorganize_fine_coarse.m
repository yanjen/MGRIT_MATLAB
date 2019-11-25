function [A, u, g] = reorganize_fine_coarse(A, u, g, fc_ratio)
%% Reorganize fine, coarse grid
% For now, assume that the first element is C-point

n = size(A, 1);
nc = ceil(n/fc_ratio);
nf = n - nc;

u = permute_fine_coarse(u, fc_ratio);
g = permute_fine_coarse(g, fc_ratio);

B = eye(n);
for i = 2:n
    if (mod(i-1, fc_ratio) == 1)
        B(i - ceil(i/fc_ratio), nf+ceil((i-1)/fc_ratio)) = A(i, i-1);
        continue;
    end
    if (mod(i, fc_ratio) == 1)
        B(nf+ceil(i/fc_ratio), i-ceil(i/fc_ratio)) = A(i, i-1);
        continue;
    end
    B(i-ceil(i/fc_ratio), i-1-ceil(i/fc_ratio)) = A(i, i-1);
end
A = sparse(B);

end