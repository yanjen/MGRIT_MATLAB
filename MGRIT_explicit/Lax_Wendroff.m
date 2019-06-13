function S = Lax_Wendroff(T, h, dt)

l = length(T);
S = zeros(1,l);

C = 300;
CR = C*dt/h;
CX = 0.5*CR*CR;

for i = 2:l-1
    S(i) = (1-2*CX)*T(i) + (-0.5*CR+CX)*T(i+1) + (0.5*CR+CX)*T(i-1);
end

end