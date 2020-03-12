function output = Lax_Wendroff_vecv(input, velocity, dt, neumann_boundary)
l = length(input);
output = input;

for i = 2:l-1
    output(i) = output(i) - (dt + 0.5 * dt * dt) * velocity(i + 1);
    output(i) = output(i) - (-dt + 0.5 * dt * dt) * velocity(i - 1);
    output(i) = output(i) + dt * dt * velocity(i);
end

if neumann_boundary == 1
    output(1) = output(1) - dt * dt * velocity(i + 1);
    output(1) = output(1) + dt * dt * velocity(i);
end

if neumann_boundary == 2
    output(1) = output(1) - (dt + 0.5 * dt * dt) * velocity(i + 1);
    output(1) = output(1) + dt * dt * velocity(i);
end

end