function [U1, U2, U3] = extrapolate_boundary(U1, U2, U3)

U_boundary = 0;

R = 2 * U1(2) - U1(3);
% P = 2 * P(2) - P(3);
% E = P / (gamma - 1);
E = 2 * U3(2) - U3(3);

U1(1) = R;
U2(1) = R * U_boundary;
U3(1) = E;

end