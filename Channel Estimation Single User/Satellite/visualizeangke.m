% Parameters for illustration
theta_max = 30 * pi / 180; % Max angle
num_points = 500; % Generate multiple points

% Generate random points
xi = zeros(num_points, 1);
yi = zeros(num_points, 1);
for i = 1:num_points
    theta_x = theta_max * rand();
    theta_y = theta_max * rand();
    x = sin(theta_y) * cos(theta_x);
    y = sin(theta_y) * sin(theta_x);
    if x^2 + y^2 <= sin(theta_max)^2
        xi(i) = x;
        yi(i) = y;
    end
end

% 1. Plot Points on the Circle Region
figure;
scatter(xi, yi, 10, 'filled');
title('Generated Points within Circular Region');
xlabel('x'); ylabel('y');
axis equal;

% 2. Plot AoA and AoD in Polar Coordinates
AoA_el = (pi * rand(1, num_points) - pi / 2);
AoD_el = acos(asin(yi));
figure;
polarplot(AoA_el, ones(1, num_points), 'r.', 'DisplayName', 'AoA');
hold on;
polarplot(AoD_el, ones(1, num_points), 'b.', 'DisplayName', 'AoD');
title('Angles of Arrival and Departure');
legend;

% 3. Path Loss Visualization
SP.Re = 6371e3; % Earth's radius (meters)
SP.Hight = 500e3; % Example height (meters)
SP.Rs = 6000e3; % Satellite radius (meters)
SP.fc = 2e9; % Carrier frequency (Hz)

alpha_new = acos((SP.Rs / SP.Re) * sin(AoD_el));
D = sqrt(SP.Re^2 .* (sin(alpha_new)).^2 + SP.Hight^2 + 2 * SP.Hight * SP.Re) - SP.Re .* sin(alpha_new);
Beta = abs(randn(1, num_points) * 3); % Random PDP
PL = 20 * log10(SP.fc) + 20 * log10(D) + 32.45 + Beta + 2;

figure;
bar(PL);
title('Path Loss (PL) for Generated Points');
xlabel('Point Index');
ylabel('Path Loss (dB)');

% 4. PDP Distribution
figure;
histogram(Beta);
title('Power Delay Profile (PDP)');
xlabel('PDP (\beta)');
ylabel('Frequency');

