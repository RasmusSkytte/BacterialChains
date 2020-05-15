% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_S3', 'dir')
    mkdir ../figures/Figure_S3
end

% Define the simulations
N = 1000;
theta = 0.2:0.2:1;

repeats = 100;

% Make arrays to store the random walk distance
x = nan(numel(N), repeats, numel(theta), 2*max(N)-1);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for i = 1:numel(N)
    for t = 1:numel(theta)
        for j = 1:repeats

            % Define the path for the data
            dpath = sprintf('%s/Chain/N_%d/theta_%.3f_pi', path, N(i), theta(t));

            % Check if data exists
            if exist(sprintf('%s/repeat_%d/CellData.txt', dpath, j-1), 'file')
                data = importdata(sprintf('%s/repeat_%d/RandomWalkDistance.txt', dpath, j-1));
            else
                continue
            end

            % Store the fitted values
            x(i, j, t, 1:numel(data)) = data;

        end
    end
end

% Compute RMS
d = squeeze(mean(x, 2));
e = squeeze(std(x, [], 2)) / sqrt(repeats);

% Plot assuming RW scaling
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.75];
cc = lines(numel(theta));

for t = 1:numel(theta)
    xx = 1:(2*N(i)-1);
    yy = squeeze(d(t, :));
    ee = squeeze(e(t, :));
    h = errorbar(ax, xx.^(0.5), yy, ee, 's', 'MarkerFaceColor', 'Auto', 'Color', cc(t, :), 'LineWidth', 1.5);
    axes(ax);
    hs = shadedErrorBar(h.XData, h.YData, [h.YPositiveDelta; h.YNegativeDelta], 'lineProps', {'Color', h.Color, 'LineWidth', 1.5});
    delete(h);
    plot([xx(1), xx(end)].^(0.5), [yy(1), yy(end)], 'k', 'LineWidth', 1)
end

xlabel(ax, 'Links^{0.5}')
ylabel(ax, 'Distance ({\mu}m)')

ax.XLim = [0 xx(end)^0.5];

ax.LineWidth = 1.5;
ax.FontSize = 16;

ax.YLabel.Position(1) = 1.2*ax.YLabel.Position(1);

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
saveas(fh, '../figures/Figure_S3/FigS3a.png')

% Plot assuming Flory scaling
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.75];
cc = lines(numel(theta));

for t = 1:numel(theta)
    xx = (1:(2*N(i)-1)) / 2;
    yy = squeeze(d(t, :));
    ee = squeeze(e(t, :));
    h = errorbar(ax, xx.^(0.588), yy, ee, 's', 'MarkerFaceColor', 'Auto', 'Color', cc(t, :), 'LineWidth', 1.5);
    axes(ax);
    hs(t) = shadedErrorBar(h.XData, h.YData, [h.YPositiveDelta; h.YNegativeDelta], 'lineProps', {'Color', h.Color, 'LineWidth', 1.5});
    hs(t).mainLine.DisplayName = sprintf('\\Theta / 2 = %.1f \\pi', theta(t) / 2);
    delete(h);
    plot([xx(1), xx(end)].^(0.588), [yy(1), yy(end)], 'k', 'LineWidth', 1)
end

xlabel(ax, 'Links^{0.588}')
ylabel(ax, 'Distance ({\mu}m)')

ax.XLim = [0 xx(end)^0.588];

ax.LineWidth = 1.5;
ax.FontSize = 16;

ax.YLabel.Position(1) = 1.2*ax.YLabel.Position(1);

legend([hs.mainLine], 'Location', 'NorthEastOutside')
ax.Position = [0.13 0.17 0.613 0.75];

pause(0.1); fh.Position = [10 50 725.5 420]; pause(0.1);
saveas(fh, '../figures/Figure_S3/FigS3b.png')