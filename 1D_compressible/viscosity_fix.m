function [rhou, E] = viscosity_fix(rho, rhou, E, P, gamma, dx, dt)

C2ref = 0.51;
mu = 0.2;
PRL = 0.72;

u = rhou ./ rho;
l = length(rho);
drho_dx = zeros(1,l-1);
T = gamma .* P ./ rho;
dT_dx = zeros(1,l-1);
uc = zeros(1,l-1);

temp = gamma .* P ./ rho;
v = temp .^ 1.5 * (1 + C2ref) * mu ./ (temp + C2ref);
viscosity = zeros(1,l-1);

for i = 1:l-1
    drho_dx(i) = (rho(i+1) - rho(i)) / dx;
    dT_dx(i) = (T(i+1) - T(i)) / dx;
    uc(i) = (u(i+1) + u(i)) / 2;
    viscosity(i) = (v(i+1) + v(i)) / 2;
end

tau_xx = 4 / 3 .* viscosity .* drho_dx;
Qx = -viscosity ./ PRL ./ (gamma - 1) .* mu;

Grad_rhou = tau_xx;
Grad_E = -Qx;

for i = 1:l-1
    rhou(i) = rhou(i) + dt * Grad_rhou(i);
    rhou(i+1) = rhou(i+1) - dt * Grad_rhou(i);
    E(i) = E(i) + dt * Grad_E(i);
    E(i+1) = E(i+1) - dt * Grad_E(i);
end

end