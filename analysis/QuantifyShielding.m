% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_6', 'dir')
    mkdir ../figures/Figure_6
end

% Define the simulations
N = [3 10 32 1000 316 1000];
theta = 0.2:0.2:1;

repeats = 9 * ones(size(N));

% Make arrays to store the test
ratios = nan(numel(N), max(repeats), numel(theta) + 2);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for n = 1:numel(N)

    % Define null hypthesis CDF
    nullCDF = (1 / N(n)):(1 / N(n)):1;

    for t = 1:(numel(theta) + 2)
        hits = [];
        for j = 1:repeats(n)

            % Define the path for the data
            if t == numel(theta) + 1
                dpath = sprintf('%s/WellMixed/N_%d', path, N(n));
            elseif t == numel(theta) + 2
                dpath = sprintf('%s/SphericalColony/N_%d', path, N(n));
            else
                dpath = sprintf('%s/Chain/N_%d/theta_%.3f_pi', path, N(n), theta(t));
            end

            % Check if data exists
            if exist(sprintf('%s/repeat_%d/PhageLocation.txt', dpath, j-1), 'file')
                data = importdata(sprintf('%s/repeat_%d/PhageLocation.txt', dpath, j-1));
            else
                continue
            end

            % Compute hits
            hits = histcounts(data(:, end), 0.5:(N(n) + 0.5));

            % Store the extremal ratio of hits
            ratios(n, j, t) = min(hits) / max(hits);

            if ( (t == 2 || t == numel(theta) + 2 ) && n == 1 && j == 1)
                fh = figure();
                ax = axes;
                scatter3(ax, data(:, 1), data(:, 2), data(:, 3), '.r', 'SizeData', 10);
                axis equal;
                ax.Visible = 'off';
                if t == 2
                    ax.Position = [0 0 1 1];
                    view([3 1 1])
                elseif t == numel(theta) + 2
                    ax.Position = [0.2 0.2 0.6 0.6];
                    view([-0.2 -11 6.5])
                end

                pause(0.1); fh.Position(1:3) = [10 50 280]; pause(0.1);
                if t == numel(theta) + 2
                    print(fh, '../figures/Figure_6/Fig6c.tif', '-dtiff', '-r900')
                else
                    print(fh, '../figures/Figure_6/Fig6b.tif', '-dtiff', '-r900')
                end

            end
        end
    end
end

% Plot the degree of shielding
fh = figure;
fh.Resize = 'off';
ax1 = axes;
ax2 = axes;
ax3 = axes;

ax1.NextPlot = 'add';
ax2.NextPlot = 'add';
ax3.NextPlot = 'add';

ax1.Box = 'on';
ax2.Box = 'on';
ax3.Box = 'on';

ax1.Position = [0.16  0.15 0.1    0.82];
ax2.Position = [0.325 0.15 0.4875 0.82];
ax3.Position = [0.87  0.15 0.1    0.82];

cc = lines(numel(theta) + 2);
for n = [1 2 3]

    for t = 1:(numel(theta) + 2)

        tt = ratios(n, :, t);

        if t <= numel(theta)
            errorbar(ax2, theta(t) / 2, nanmean(tt, 2), nanstd(tt, [], 2) / sqrt(sum(~isnan(tt))), 's', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        elseif t == numel(theta) + 1
            errorbar(ax1, 0.5, nanmean(tt, 2), nanstd(tt, [], 2) / sqrt(sum(~isnan(tt))), '^', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        elseif t == numel(theta) + 2
            errorbar(ax3, 0.5, nanmean(tt, 2), nanstd(tt, [], 2) / sqrt(sum(~isnan(tt))), 'o', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        end

    end

end

ylabel(ax1, 'min(hits) / max(hits)')
xlabel(ax2, '\Theta / 2\pi')

ax1.XTickLabel = [];
ax3.XTickLabel = [];

text(ax1, 0.5, -0.09, {'Well', 'Mixed'}, 'HorizontalAlignment', 'Center', 'FontSize', 12);
text(ax3, 0.5, -0.09, {'Spherical', 'Colony'}, 'HorizontalAlignment', 'Center', 'FontSize', 12);

ax1.XTickLabelRotation = 45;
ax3.XTickLabelRotation = 45;

ax2.YTickLabel = [];
ax3.YTickLabel = [];

ax1.LineWidth = 1.5;
ax2.LineWidth = 1.5;
ax3.LineWidth = 1.5;

ax1.FontSize = 16;
ax2.FontSize = 16;
ax3.FontSize = 16;

ax1.XLim = [0 1];
ax2.XLim = [(theta(1)/2) 0.5] - [0.05 -0.05];
ax3.XLim = [0 1];

ax1.XTick = [0 1];
ax2.XTick = 0:0.1:0.5;
ax3.XTick = [0 1];

ax1.YLim = [0 1];
ax2.YLim = [0 1];
ax3.YLim = [0 1];

xtickformat(ax2, '%.1f')
ytickformat(ax1, '%.1f')

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
print(fh, '../figures/Figure_6/Fig6a.tif', '-dtiff', '-r900')