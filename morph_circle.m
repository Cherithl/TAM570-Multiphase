function [X1, Y1] = morph_circle(X, Y, Nx, Ny, N)

% X size (N, E, N)
% Y size (N, E, N)

X1 = zeros(size(X));
Y1 = zeros(size(Y));

% inner and outer radii
R1 = 1;
R2 = 3;

Theta1 = pi/2;
Theta2 = 5*pi/2;

dr = (R2 - R1) / Ny;
dtheta = (Theta1 - Theta2) / Nx;

E = Nx * Ny; % total number of elements

for j = 1:Ny
    for i = 1:Nx
        e_idx = i + (j-1)*Nx;
        
        r1 = (j-1) * dr + R1;
        r2 = j * dr + R1;
        
        theta1 = (i-1) * dtheta + Theta1;
        theta2 = i * dtheta + Theta1;

        R = (r2 - r1) * (Y(:, e_idx, :) - Y(:, e_idx, 1))./(Y(:, e_idx, end) - Y(:, e_idx, 1)) + r1;
        T = (theta2 - theta1) * (X(:, e_idx, :) - X(1, e_idx, :))./(X(end, e_idx, :) - X(1, e_idx, :)) + theta1;

        X1(:, e_idx, :) = R .* cos(T);
        Y1(:, e_idx, :) = R .* sin(T);
    end 
end


figure(41);
colors = jet(E);
hold on;
for e = 1:E
    x = reshape(X(:, e, :), [], 1);
    y = reshape(Y(:, e, :), [], 1);
    scatter(x, y, 36, colors(e, :), '*'); 
end
axis square
hold off;

figure(42);
colors = jet(E);
hold on;
for e = 1:E
    x = reshape(X1(:, e, :), [], 1);
    y = reshape(Y1(:, e, :), [], 1);
    scatter(x, y, 36, colors(e, :), '*'); 
end
axis square
hold off;


