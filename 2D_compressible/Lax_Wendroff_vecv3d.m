function [rho, rhou, rhov, E, P] = Lax_Wendroff_vecv3d(rho, rhou, rhov, E, P, x, y, num_theta, gamma, dr, dtheta, dt)

l = length(rho);

u = rhou ./ rho;
v = rhov ./ rho;
P = (E - 0.5 .* rho .* u .* u) .* (gamma - 1);

F1 = rhou;
F2 = rhou .* u + P;
F3 = rhov .* u;
F4 = (E + P) .* u;
G1 = rhov;
G2 = rhou .* v;
G3 = rhov .* v + P;
G4 = (E + P) .* v;

for i = num_theta + 1 : l - num_theta
    r = sqrt(x(i) * x(i) + y(i) * y(i));
    mesh_area = dr * r * dtheta;
    i_out = i + num_theta;
    i_in = i - num_theta;
    i_cc = i + 1;
    i_c = i - 1;
    if mod(i, num_theta) == 0
        i_cc = i_cc - num_theta;
    end
    if mod(i_c, num_theta) == 0
        i_c = i_c + num_theta;
    end
    r_out = sqrt(x(i_out) * x(i_out) + y(i_out) * y(i_out));
    r_in = sqrt(x(i_in) * x(i_in) + y(i_in) * y(i_in));
    [nrx, nry] = get_norm(x(i), y(i), x(i_out), y(i_out));
    [ncx, ncy] = get_norm(x(i), y(i), x(i_c), y(i_c));
    [nccx, nccy] = get_norm(x(i), y(i), x(i_cc), y(i_cc));
    
    rho(i) = rho(i) - (dt + 0.5 * dt * dt / dr) / mesh_area * (F1(i_out) * nrx + G1(i_out) * nry) * 0.5 * (r + r_out) * dtheta;
    rho(i) = rho(i) - (dt - 0.5 * dt * dt / dr) / mesh_area * (F1(i_in) * -nrx + G1(i_in) * -nry) * 0.5 * (r + r_in) * dtheta;
    rho(i) = rho(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F1(i_c) * ncx + G1(i_c) * ncy) * dr;
    rho(i) = rho(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F1(i_cc) * nccx + G1(i_cc) * nccy) * dr;
    rho(i) = rho(i) + 0.5 * dt * dt / mesh_area * (F1(i) * (2 * nrx * r * dtheta / dr + ncx * dr / r / dtheta + nccx * dr / r / dtheta));
    rho(i) = rho(i) + 0.5 * dt * dt / mesh_area * (G1(i) * (2 * nry * r * dtheta / dr + ncy * dr / r / dtheta + nccy * dr / r / dtheta));
    rhou(i) = rhou(i) - (dt + 0.5 * dt * dt / dr) / mesh_area * (F2(i_out) * nrx + G2(i_out) * nry) * 0.5 * (r + r_out) * dtheta;
    rhou(i) = rhou(i) - (dt - 0.5 * dt * dt / dr) / mesh_area * (F2(i_in) * -nrx + G2(i_in) * -nry) * 0.5 * (r + r_in) * dtheta;
    rhou(i) = rhou(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F2(i_c) * ncx + G2(i_c) * ncy) * dr;
    rhou(i) = rhou(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F2(i_cc) * nccx + G2(i_cc) * nccy) * dr;
    rhou(i) = rhou(i) + 0.5 * dt * dt / mesh_area * (F2(i) * (2 * nrx * r * dtheta / dr + ncx * dr / r / dtheta + nccx * dr / r / dtheta));
    rhou(i) = rhou(i) + 0.5 * dt * dt / mesh_area * (G2(i) * (2 * nry * r * dtheta / dr + ncy * dr / r / dtheta + nccy * dr / r / dtheta));
    rhov(i) = rhov(i) - (dt + 0.5 * dt * dt / dr) / mesh_area * (F3(i_out) * nrx + G3(i_out) * nry) * 0.5 * (r + r_out) * dtheta;
    rhov(i) = rhov(i) - (dt - 0.5 * dt * dt / dr) / mesh_area * (F3(i_in) * -nrx + G3(i_in) * -nry) * 0.5 * (r + r_in) * dtheta;
    rhov(i) = rhov(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F3(i_c) * ncx + G3(i_c) * ncy) * dr;
    rhov(i) = rhov(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F3(i_cc) * nccx + G3(i_cc) * nccy) * dr;
    rhov(i) = rhov(i) + 0.5 * dt * dt / mesh_area * (F3(i) * (2 * nrx * r * dtheta / dr + ncx * dr / r / dtheta + nccx * dr / r / dtheta));
    rhov(i) = rhov(i) + 0.5 * dt * dt / mesh_area * (G3(i) * (2 * nry * r * dtheta / dr + ncy * dr / r / dtheta + nccy * dr / r / dtheta));
    E(i) = E(i) - (dt + 0.5 * dt * dt / dr) / mesh_area * (F4(i_out) * nrx + G4(i_out) * nry) * 0.5 * (r + r_out) * dtheta;
    E(i) = E(i) - (dt - 0.5 * dt * dt / dr) / mesh_area * (F4(i_in) * -nrx + G4(i_in) * -nry) * 0.5 * (r + r_in) * dtheta;
    E(i) = E(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F4(i_c) * ncx + G4(i_c) * ncy) * dr;
    E(i) = E(i) - (dt + 0.5 * dt * dt / r / dtheta) / mesh_area * (F4(i_cc) * nccx + G4(i_cc) * nccy) * dr;
    E(i) = E(i) + 0.5 * dt * dt / mesh_area * (F4(i) * (2 * nrx * r * dtheta / dr + ncx * dr / r / dtheta + nccx * dr / r / dtheta));
    E(i) = E(i) + 0.5 * dt * dt / mesh_area * (G4(i) * (2 * nry * r * dtheta / dr + ncy * dr / r / dtheta + nccy * dr / r / dtheta));
end

end