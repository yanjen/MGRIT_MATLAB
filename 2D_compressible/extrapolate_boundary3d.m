function [rho, rhou, rhov, E] = extrapolate_boundary3d(rho, rhou, rhov, E, P, N, gamma)

U_boundary = 0;

for i = 1:N
    i_out = i + N;
    i_cc = i + 1;
    i_c = i - 1;
    if mod(i, N) == 0
        i_cc = i_cc - N;
    end
    if mod(i_c, N) == 0
        i_c = i_c + N;
    end
    rho1 = (rho(i) + rho(i_out)) / 2;
    rho2 = (rho(i) + rho(i_c)) / 2;
    rho3 = (rho(i) + rho(i_cc)) / 2;
    P1 = (P(i) + P(i_out)) / 2;
    P2 = (P(i) + P(i_c)) / 2;
    P3 = (P(i) + P(i_cc)) / 2;
    T1 = gamma * P1 / rho1;
    T2 = gamma * P2 / rho2;
    T3 = gamma * P3 / rho3;
    T0 = (T1 + T2 + T3) / 3;
    P0 = (P1 + P2 + P3) / 3;
    rho(i) = gamma * P0 / T0;
    rhou(i) = rho(i) * U_boundary;
    rhov(i) = rho(i) * U_boundary;
    E(i) = P0 / (gamma - 1);
end

end