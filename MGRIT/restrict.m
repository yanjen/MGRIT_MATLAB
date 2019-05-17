function R = restrict(m, n, fc_ratio)
%% Compute Restrict Matrix

R = spalloc(m, n, m);
for i = 1:m
    R(m, 1 + (m-1)*fc_ratio) = 1;
end

end