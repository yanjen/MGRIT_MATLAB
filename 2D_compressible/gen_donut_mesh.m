function [x, y, dr, dtheta] = gen_donut_mesh(inner_r, outer_r, num_r, num_theta)

r = linspace(inner_r,outer_r,num_r);
theta = linspace(0,2*pi,num_theta + 1);
theta = theta(1:end - 1);
dr = r(2) - r(1);
dtheta = theta(2) - theta(1);

l = num_r * num_theta;
x = zeros(1,l);
y = zeros(1,l);

for i = 1:num_r
    for j = 1:num_theta
        x((i-1) * num_theta + j) = r(i) * cos(theta(j));
        y((i-1) * num_theta + j) = r(i) * sin(theta(j));
    end
end

end