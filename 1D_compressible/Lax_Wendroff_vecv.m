function [rho, rhou, E, P] = Lax_Wendroff_vecv(rho, rhou, E, P, gamma, dt)

l = length(rho);

u = rhou ./ rho;
F1 = rhou;
F2 = rhou .* u + P;
F3 = (E + P) .* u;
P = (E - 0.5 .* rho .* u .* u) .* (gamma - 1);

for i = 2:l-1
    rho(i) = rho(i) - (dt + 0.5 * dt * dt) * F1(i + 1);
    rho(i) = rho(i) - (-dt + 0.5 * dt * dt) * F1(i - 1);
    rho(i) = rho(i) + dt * dt * F1(i);
    rhou(i) = rhou(i) - (dt + 0.5 * dt * dt) * F2(i + 1);
    rhou(i) = rhou(i) - (-dt + 0.5 * dt * dt) * F2(i - 1);
    rhou(i) = rhou(i) + dt * dt * F2(i);
    E(i) = E(i) - (dt + 0.5 * dt * dt) * F3(i + 1);
    E(i) = E(i) - (-dt + 0.5 * dt * dt) * F3(i - 1);
    E(i) = E(i) + dt * dt * F3(i);
end

end