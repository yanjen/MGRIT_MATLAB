function T = implicit_time_marching(T, Q, h, dt)

l = length(T);
A = conv2(eye(l-1), [-dt/h^2 (1+2*dt/h^2) -dt/h^2], 'same');
A(l-1,l-2) = A(l-1,l-2) - dt/h^2;

b = T(2:end)';
b(l-1) = b(l-1) + Q*dt;

x = A\b;
T = [0;x]';

end