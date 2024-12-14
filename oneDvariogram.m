function oneDvariogram(x,z)
%x = linspace(0, 10, 20);      % 1D spatial locations
%z = sin(x) + 0.1 * randn(1, 20); % Example output values with noise

% Compute pairwise distances and semivariances
% Compute pairwise distances and semivariances
n = length(x);
distances = zeros(n, n);
semivariances = zeros(n, n);

for i = 1:n
    for j = 1:n
        distances(i, j) = abs(x(i) - x(j)); % Pairwise distances
        semivariances(i, j) = 0.5 * (z(i) - z(j))^2; % Semivariance
    end
end

% Flatten the matrices
distances = distances(:);
semivariances = semivariances(:);

% Remove zero distances (self-comparisons)
nonzeroIdx = distances > 0;
distances = distances(nonzeroIdx);
semivariances = semivariances(nonzeroIdx);

% Plot the variogram
figure;
scatter(distances, semivariances, '. b');
xlabel('Distance (Lag)');
ylabel('Semivariance');
title('Cloud Variogram');
grid minor;
end