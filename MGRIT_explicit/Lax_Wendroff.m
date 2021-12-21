function S = Lax_Wendroff(U, g, h, dt)
%% U^{i+1} = \Psi(U^i)+g

l = length(U);
S = zeros(1,l);

C = 1;
CR = C*dt/h;
CX = 0.5*CR*CR;

for i = 2:l-1
    S(i) = (1-2*CX)*U(i) + (-0.5*CR+CX)*U(i+1) + (0.5*CR+CX)*U(i-1);
end

S(1) = (1-2*CX)*U(1) + (-0.5*CR+CX)*U(2) + (0.5*CR+CX)*U(l);
S(l) = (1-2*CX)*U(l) + (-0.5*CR+CX)*U(1) + (0.5*CR+CX)*U(l-1);

S = S + g;

end