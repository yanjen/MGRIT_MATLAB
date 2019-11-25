function u = backward_euler(u, h, dt)

l = length(u);
A = conv2(eye(l), [-dt/(2*h) 1 dt/(2*h)], 'same');
A(1,2) = 0;
A(l,l-1) = 0;

b = u';
x = A\b;

u = x';

end