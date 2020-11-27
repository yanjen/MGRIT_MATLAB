function [nx, ny] = get_norm(x, y, x_ref, y_ref)

nx = x_ref - x;
ny = y_ref - y;
nrms = sqrt(nx * nx + ny * ny);
nx = nx / nrms;
ny = ny / nrms;

end