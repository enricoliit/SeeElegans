function [planeCoeff, normalVector] = fitPlaneAndPlot_2(points)
    % At least 3 points
    if size(points, 1) < 3 || size(points, 2) ~= 3
        error('Deve essere fornito un array con almeno 3 punti e dimensione Nx3.');
    end

    % Get coordinates x, y, z
    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);

    % Linear regression
    % Model: z = Ax + By + C
    [coeff, ~, ~, ~, ~] = regress(z, [x, y, ones(length(x), 1)]);

    % Plane coefficients (A, B, C, D)
    A = coeff(1);
    B = coeff(2);
    C = -1;
    D = coeff(3);
    planeCoeff = [A, B, C, D];

    % Normal vector to the plane
    normalVector = [A, B, C];
    normalVector = normalVector/norm(normalVector);
    normalVector = normalVector';

    % Plotting
    [X, Y] = meshgrid(linspace(min(x), max(x), 10), linspace(min(y), max(y), 10));
    Z = A*X + B*Y + D;

    plot3(x, y, z, 'ro');
    hold on;
    surf(X, Y, Z, 'FaceAlpha', 0.5);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Piano di adattamento e punti dati');
    hold off;
end
