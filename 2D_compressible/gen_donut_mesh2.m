function [x, y, r, dtheta] = gen_donut_mesh2(inner_r, outer_r, num_theta)

theta = linspace(0,2*pi,num_theta + 1);
theta = theta(1:end - 1);
dtheta = theta(2) - theta(1);
num_r = floor((log(outer_r+1e-10) - log(inner_r)) / log(1+dtheta)) + 2;
r = zeros(1,num_r);
r(1) = inner_r;

x = zeros(num_theta, num_r);
y = zeros(num_theta, num_r);

for i = 1:num_r
    for j = 1:num_theta
        x(j, i) = r(i) * cos(theta(j));
        y(j, i) = r(i) * sin(theta(j));
    end
    if i == num_r - 1
        r(num_r) = outer_r;
    elseif i < num_r - 1
        r(i + 1) = r(i) * (1 + dtheta);
    end
end

end