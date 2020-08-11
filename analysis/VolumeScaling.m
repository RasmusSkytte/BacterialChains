% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_S2', 'dir')
    mkdir ../figures/Figure_S2
end

if ~exist('../figures/Figure_S6', 'dir')
    mkdir ../figures/Figure_S6
end

% Time-steps to check
theta = 1:-0.1:0.2;
repeats = 9;
N = 100;

% Make arrays to store the fitted values inmat
angle   = nan(length(theta), repeats);
d_angle = nan(length(theta), repeats);

% Prepare figure
fh1 = figure(1); clf; hold on; box on;
ax1 = gca;
ax1.Position = [0.17 0.17 0.79 0.8];

% Allocate array to store the angles
angles  = nan((N-1) * repeats, numel(theta));

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Prepare figures
cc = lines(numel(theta) + 2);
cc = cc(end-2:-1:1, :);

fh2 = figure;
fh.Resize = 'off';
ax2 = gca;
ax2.Position = [0.17 0.17 0.79 0.8];
ax2.NextPlot = 'add';
ax2.Box = 'on';

% Loop over runs
for t = 1:numel(theta)

    % Allocate array to store the distances
    distances = nan(N, repeats);

    for r = 1:repeats

        % Load data from the run
        data = importdata(sprintf('%s/BendingAngleValidation/theta_%.3f_pi/repeat_%d/CellData.txt', path, theta(t), r-1));

        % Open the log to get the meta data
        log = importdata(sprintf('%s/BendingAngleValidation/theta_%.3f_pi/repeat_%d/log.txt', path, theta(t), r-1));

        % Allocate array to store the centers
        centers = nan(N, 3);

        % Loop over cells to compute the angles
        for n = 1:N
        P = data(n, 2:4);
        Q = data(n, 5:7);

        centers(n, :) = (P + Q) / 2;

        if n < N
            R = data(n + 1, 2:4);
            S = data(n + 1, 5:7);

            PQ = P-Q;
            RS = R-S;
            angles(n + (r-1) * (N-1), t) = real(acos(abs(PQ * RS') / (norm(PQ) * norm(RS))) / pi);
        end

        end

        % Adjust centers
        [~, I] = min(sum((centers - repmat(mean(centers), N, 1)).^2, 2));
        centers = centers - repmat(centers(I, :), N, 1);
        distances(:, r) = sqrt(sum(centers.^2, 2));
        s = linspace(0, 125, 20);

    end
    if mod(t + 1, 2) == 0
        errorbar(ax2, s, mean(sum(reshape(distances, N, 1, repeats)<s), 3) / N, std(sum(reshape(distances, N, 1, repeats)<s), [], 3) / (N * sqrt(repeats)), '-o', 'MarkerSize', 4, 'Color', cc(t, :), 'MarkerFaceColor', 'Auto', 'DisplayName', sprintf('\\Theta / 2 = %.1f \\pi', theta(t) / 2))
    end
end

% Update the volume scaling figure
s = linspace(0, 125, 100);
p1 = plot(ax2, s, min(1, 2 * s / (N * (2 * log.data(10) + log.data(11)))), '--k', 'LineWidth', 2, 'DisplayName', 'Straight line');
p2 = plot(ax2, s, min(1, 4 * pi / 3 * s.^3 / (N * (4 * pi / 3 * log.data(10).^3  +  pi * log.data(10)^2 * log.data(11)))), ':k', 'LineWidth', 2, 'DisplayName', 'Solid sphere');
uistack(p1, 'bottom')
uistack(p2, 'bottom')
xlabel(ax2, 'r ({\mu}m)')
ylabel(ax2, 'V(r) / V_{max}')
legend(ax2, 'location', 'NorthEastOutside')
ax2.LineWidth = 1.5;
ax2.FontSize = 16;
ax2.XLim = [0 125];
ax2.XTick = 0:25:125;

pause(0.1); fh2.Position = [10 50 750 420]; pause(0.1);
print(fh2, '../figures/Figure_S6/FigS6.tif', '-dtiff', '-r900')

% Plot the convergence
errorbar(ax1, theta, mean(angles(:, :)), std(angles(:, :)) / sqrt(repeats * (N-1)), 'kx', 'LineWidth', 2, 'MarkerSize', 2, 'MarkerFaceColor', 'Auto', 'DisplayName', sprintf('n = %d', repeats))

plot(ax1, [0 1], [0 0.5], '--k', 'DisplayName', '<\theta> = \Theta / 2')
xlabel(ax1, '\Theta / \pi')
ylabel(ax1, '<\theta> / \pi')
% legend(ax1, 'location', 'southeast')
ax1.LineWidth = 1.5;
ax1.FontSize = 16;

pause(0.1); fh1.Position = [10 50 560 420]; pause(0.1);
print(fh1, '../figures/Figure_S2/FigS2.tif', '-dtiff', '-r900')
