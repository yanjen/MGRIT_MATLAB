function T = explicit_time_marching(T, Q, h, dt)

l = length(T);
T = [T, T(l-1)];
dTdt = zeros(1, l);

for i = 1:l-1
    dTdt(i+1) = (T(i) - 2*T(i+1) + T(i+2))/h^2 + Q;
end
% dTdt(l) = dTdt(l) + Q;

T = T(1:l) + dTdt * dt;

end