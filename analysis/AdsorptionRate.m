% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_5', 'dir')
    mkdir ../figures/Figure_5
end

if ~exist('../figures/Figure_S8', 'dir')
    mkdir ../figures/Figure_S8
end

% Define the simulations
N = [1 3 10 32 100 316 1000];
theta = 0.2:0.2:1;

repeats = 9 * ones(size(N));

% Make arrays to store the fitted values in
f_eta   = nan(numel(N), max(repeats), numel(theta) + 2);
hits    = nan(numel(N), max(repeats), numel(theta) + 2);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for i = 1:numel(N)
    for t = 1:(numel(theta) + 2)
        for j = 1:repeats(i)

            % Define the path for the data
            if t == numel(theta) + 1
                dpath = sprintf('%s/WellMixed/N_%d', path, N(i));
            elseif t == numel(theta) + 2
                dpath = sprintf('%s/SphericalColony/N_%d', path, N(i));
            else
                dpath = sprintf('%s/Chain/N_%d/theta_%.3f_pi', path, N(i), theta(t));
            end

            % Check if data exists
            if exist(sprintf('%s/repeat_%d/Completed.txt', dpath, j-1), 'file')
                data = importdata(sprintf('%s/repeat_%d/ColonySize.txt', dpath, j-1));
                flog = importdata(sprintf('%s/repeat_%d/log.txt', dpath, j-1), '=');
            else
                continue
            end

            % Define the simulation volume
            V = flog.data(5).^3;

            % Get time and phage numbers
            T = data(:, 1);
            P = data(:, 2);

            % Store the number of hits
            hits(i, j, t) = P(1) - P(end);

            % Perform fit
            I = ceil(size(data, 1) / 2):size(data, 1);
            f = fit(T(I) - T(I(1)), P(I) / P(I(1)), 'exp(-a*x)', 'StartPoint', 1e5 / V);

            % Store the fitted values
            f_eta(i, j, t) = V * f.a;
        end
    end
end

% Detect completed runs for these graphs
I = sum(squeeze(sum(~isnan(f_eta),2) >= 3), 2) > 0;
N = N(I);

% Get the eta's
eta = squeeze(nanmean(f_eta(I, :, :), 2))';
d_eta = squeeze(nanstd(f_eta(I, :, :), [], 2))' ./ squeeze(sqrt(sum(~isnan(f_eta(I, :, :)), 2)))';

% reference point
ref = 1;

% Compute relative eta and error
eta_0 = repmat(eta(:, ref), 1, numel(N));
d_eta_0 = repmat(d_eta(:, ref), 1, numel(N));

r_eta   = eta ./ eta_0;
d_r_eta = r_eta .* sqrt( (d_eta ./ eta).^2 + (d_eta_0 ./ eta_0).^2);

% Store the AdsorptionRate
save('AdsorptionRate.mat', 'eta', 'd_eta', 'f_eta')

% Test for comparability
fh = figure;
fh.Resize = 'off';
ax1 = axes;
ax2 = axes;
ax3 = axes;

ax2.NextPlot = 'add';

ax1.Box = 'on';
ax2.Box = 'on';
ax3.Box = 'on';

ax1.Position = [0.16  0.15 0.1    0.78];
ax2.Position = [0.325 0.15 0.4875 0.78];
ax3.Position = [0.87  0.15 0.1    0.78];

cc = lines(numel(theta) + 2);
for t = 1:(numel(theta) + 2)

    if t <= numel(theta)
        errorbar(ax2, theta(t) / 2, eta(t, 1), d_eta(t, 1), 's', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'auto')
    elseif t == numel(theta) + 1
        errorbar(ax1, 0.5,          eta(t, 1), d_eta(t, 1), '^', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'auto')
    elseif t == numel(theta) + 2
        errorbar(ax3, 0.5,          eta(t, 1), d_eta(t, 1), 'o', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'auto')
    end

end

ylabel(ax1, '\eta_1 ({\mu}m^3/h)')
xlabel(ax2, '\Theta / 2\pi')

ax1.XTickLabel = [];
ax3.XTickLabel = [];

oom = mode(10.^floor(log10(eta(:, 1))));
text(ax1, 0.5, (mean(round(eta(:, 1) / oom, 1))-0.0575)*oom, {'Well', 'Mixed'}, 'HorizontalAlignment', 'Center', 'FontSize', 11);
text(ax3, 0.5, (mean(round(eta(:, 1) / oom, 1))-0.0575)*oom, {'Spherical', 'Colony'}, 'HorizontalAlignment', 'Center', 'FontSize', 11);

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

ax1.YLim = (nanmean(round(eta(:, 1) / oom, 1)) + [-0.05 0.05])*oom;
ax2.YLim = (nanmean(round(eta(:, 1) / oom, 1)) + [-0.05 0.05])*oom;
ax3.YLim = (nanmean(round(eta(:, 1) / oom, 1)) + [-0.05 0.05])*oom;

xtickformat(ax2, '%.1f')
ytickformat(ax1, '%.2f')

ax1.YLabel.Position(1) = -0.8;

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
saveas(fh, '../figures/Figure_S8/FigS8a.png')


% Plot the adsorption rates
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% Draw grey areas
x = [1 N(end)*1.1 1];
y = [1 N(end)*1.1 N(end)*1.1];
patch(x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')

x = [1  N(end)*1.1        N(end)*1.1];
y = [1 (N(end)*1.1)^(1 / 3) 1];
patch(x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')

% Draw reference lines
plot(ax, [1 N(end)*1.1],  [1 N(end)*1.1] / N(ref),           'k:',  'LineWidth', 2, 'DisplayName', 'Well-mixed')
plot(ax, [1 N(end)*1.1], ([1 N(end)*1.1] / N(ref)).^(0.500), 'k',   'LineWidth', 1, 'DisplayName', 'Gaussian Chain')
plot(ax, [1 N(end)*1.1], ([1 N(end)*1.1] / N(ref)).^(1 / 3), 'k--', 'LineWidth', 2, 'DisplayName', 'Sphere')


% Plot the measured adsorption rates
for t = 1:numel(theta) + 2
    r = min(sum(~isnan(f_eta(:, :, t)), 2));
    if t <= numel(theta)
        name = sprintf('<\\theta> = %.1f \\pi (n = %d)', theta(t) / 2, r);
        symbol = 's';
    elseif t == numel(theta) + 1
        name = sprintf('Well-Mixed (n = %d)', r);
        symbol = '^';
    elseif t == numel(theta) + 2
        name = sprintf('Colony (n = %d)', r);
        symbol = 'o';
    end
    e = errorbar(ax, N(1:end), r_eta(t, :), d_r_eta(t, :), symbol, 'LineWidth', 1.5, 'DisplayName', name, 'MarkerFaceColor', 'auto');
end

xlabel(ax, 'N')
ylabel(ax, sprintf('\\eta / \\eta_{%d}', N(ref)))
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');

xlim(ax, [1, N(end) * 1.1])
ylim(ax, [1 / N(ref), N(end) / N(ref) * 1.1])

ax.XTick = N(1:2:end);

ax.LineWidth = 1.5;
ax.FontSize = 16;

set(ax,'Layer','top');
pause(0.1); fh.Position = [10 50 480 420]; pause(0.1);
saveas(fh, '../figures/Figure_5/Fig5a.png')


% Plot the adsorption scaling
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

ax1.Position = [0.15  0.17 0.1    0.8];
ax2.Position = [0.315 0.17 0.4875 0.8];
ax3.Position = [0.86  0.17 0.1    0.8];

% Draw grey areas
x = [0 1 1 0];
y = [1 1 2 2];
patch(ax1, x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')
patch(ax3, x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')

y = [1 1 0 0] * 1/3;
patch(ax1, x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')
patch(ax3, x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')

x = [-1 1 1 -1];
y = [1 1 2 2];
patch(ax2, x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')

y = [1 1 0 0] * 1/3;
patch(ax2, x, y, [0.9 0.9 0.9], 'EdgeColor', 'None')

% Draw reference lines
plot(ax1, [0 1], [1 1],         ':k',  'LineWidth', 2, 'DisplayName', 'Well Mixed')
plot(ax1, [0 1], [1 1] * 1 / 3, '--k', 'LineWidth', 2, 'DisplayName', 'Sphere')
plot(ax1, [0 1], [1 1] * 0.500, 'k',   'LineWidth', 1, 'DisplayName', 'Gaussian Chain')

plot(ax2, [0 1] - [0.05 -0.05], [1 1],         ':k',  'LineWidth', 2, 'DisplayName', 'Well Mixed')
plot(ax2, [0 1] - [0.05 -0.05], [1 1] * 1 / 3, '--k', 'LineWidth', 2, 'DisplayName', 'Sphere')
plot(ax2, [0 1] - [0.05 -0.05], [1 1] * 0.500, 'k',   'LineWidth', 1, 'DisplayName', 'Gaussian Chain')

plot(ax3, [0 1], [1 1],         ':k',  'LineWidth', 2, 'DisplayName', 'Well Mixed')
plot(ax3, [0 1], [1 1] * 1 / 3, '--k', 'LineWidth', 2, 'DisplayName', 'Sphere')
plot(ax3, [0 1], [1 1] * 0.500, 'k',   'LineWidth', 1, 'DisplayName', 'Gaussian Chain')

% Measure and plot the adsorption exponents
exponent = nan(size(theta));
d_exponent = nan(size(theta));
cc = lines(numel(theta) + 2);
for t = 1:(numel(theta) + 2)

    I = isnan(eta(t, :));

    if all(I)
        continue
    end

    f = fit(N(~I)', eta(t, ~I)', 'power1', 'weight', (eta(t, ~I)' ./ d_eta(t, ~I)').^2, 'StartPoint', [1 1]);
    exponent(t) = f.b;

    c = confint(f, 0.68);
    d_exponent(t) = diff(c(:, 2)) / 2;

    if t <= numel(theta)
        errorbar(ax2, theta(t) / 2, exponent(t), d_exponent(t), 's', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'auto')
    elseif t == numel(theta) + 1
        errorbar(ax1, 0.5,          exponent(t), d_exponent(t), '^', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'auto')
    elseif t == numel(theta) + 2
        errorbar(ax3, 0.5,          exponent(t), d_exponent(t), 'o', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'auto')
    end

end

ylabel(ax1, '\gamma')
xlabel(ax2, '\Theta / 2\pi')

ax1.XTickLabel = [];
ax3.XTickLabel = [];

text(ax1, 0.5, 0.225, {'Well', 'Mixed'}, 'HorizontalAlignment', 'Center', 'FontSize', 11);
text(ax3, 0.5, 0.225, {'Spherical', 'Colony'}, 'HorizontalAlignment', 'Center', 'FontSize', 11);

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

pause(0.1);
ax1.YLim = [0.3 1.05];
ax2.YLim = [0.3 1.05];
ax3.YLim = [0.3 1.05];

xtickformat(ax2, '%.1f')
ytickformat(ax1, '%.1f')

set(ax1,'Layer','top');
set(ax2,'Layer','top');
set(ax3,'Layer','top');
pause(0.1); fh.Position = [10 50 480 420]; pause(0.1);
saveas(fh, '../figures/Figure_5/Fig5b.png')


% Plot the number of hits
fh = figure;
fh.Resize = 'off';
ax1 = axes;
ax2 = axes;
ax3 = axes;

ax1.Box = 'on';
ax2.Box = 'on';
ax3.Box = 'on';

hits = hits / 1e5;
imagesc(ax1, 0.5,       1:numel(N), squeeze(nanmedian(hits(:, :, end-1), 2)));
imagesc(ax2, theta / 2, 1:numel(N), squeeze(nanmedian(hits(:, :, 1:numel(theta)), 2)));
imagesc(ax3, 0.5,       1:numel(N), squeeze(nanmedian(hits(:, :, end), 2)));
chits = squeeze(nanmedian(hits(:, :, :), 2));
ax1.CLim = [0.055 0.105];
ax2.CLim = [0.055 0.105];
ax3.CLim = [0.055 0.105];

c = colorbar(ax3);
c.LineWidth = 1;
c.Ticks = 0.06:0.01:1;
c.TickLabels = cellfun(@(s)sprintf('%.2f', str2double(s)), c.TickLabels, 'UniformOutput', false);

ylabel(ax1, 'N')
xlabel(ax2, '\Theta / 2\pi')

ax1.XTickLabel = [];
ax3.XTickLabel = [];

oom = mode(10.^floor(log10(eta(:, 1))));
text(ax1, 0.5, 0, {'Well', 'Mixed'}, 'HorizontalAlignment', 'Center', 'FontSize', 11);
text(ax3, 0.5, 0, {'Spherical', 'Colony'}, 'HorizontalAlignment', 'Center', 'FontSize', 11);

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
ax2.XLim = [(theta(1)/2) 0.5] + [-0.05 0.05];
ax3.XLim = [0 1];

ax1.XTick = [0 1];
ax2.XTick = 0:0.1:0.5;
ax3.XTick = [0 1];

ax1.YTick = 1:2:numel(N);
ax2.YTick = 1:2:numel(N);
ax3.YTick = 1:2:numel(N);

ax1.YLim = [1 numel(N)] + [-0.5 0.5];
ax2.YLim = [1 numel(N)] + [-0.5 0.5];
ax3.YLim = [1 numel(N)] + [-0.5 0.5];

ax1.YDir = 'Normal';
ax2.YDir = 'Normal';
ax3.YDir = 'Normal';

xtickformat(ax2, '%.1f')
ytickformat(ax1, '%.1f')

ax1.YTickLabel = cellfun(@(N)sprintf('10^%d', N), num2cell(log10(N(ax1.YTick))), 'UniformOutput', false);

ax1.Position = [0.15  0.15 0.09  0.8];
ax2.Position = [0.285 0.15 0.439 0.8];
ax3.Position = [0.77  0.15 0.09  0.8];
c.Position([1 3]) = [0.88 0.025];

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
saveas(fh, '../figures/Figure_S8/FigS8b.png')


% Store the exponent
save('AdsorptionExponent.mat', 'exponent', 'd_exponent')