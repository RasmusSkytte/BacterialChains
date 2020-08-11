% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_S7', 'dir')
    mkdir ../figures/Figure_S7
end

% Time-steps to check
dT = 4:0.5:7.5;
repeats = 3;

% Make arrays to store the fitted values in
eta = nan(length(dT), repeats);
d_eta = nan(length(dT), repeats);

fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for t = 1:numel(dT)
    for r = 1:repeats

        % Load data from the run
        data = importdata(sprintf('%s/TimeStepEstimation/dT_1e-%.2f/repeat_%d/ColonySize.txt', path, dT(t), r-1));

        % Get time and phage numbers
        T = data(:, 1);
        P = data(:, 2);

        % Open the log to get the meta data
        flog = importdata(sprintf('%s/TimeStepEstimation/dT_1e-%.2f/repeat_%d/log.txt', path, dT(t), r-1));
        V = flog.data(5).^3;

        % Perform fit
        I = ceil(size(data, 1) / 2):size(data, 1);
        f = fit(T(I) - T(I(1)), P(I) / P(I(1)), 'exp(-a*x)', 'StartPoint', 1e5 / V);

        % Store the fitted values
        eta(t, r) = f.a * V;
        d_eta(t, r) = diff(confint(f, 0.68)) / 2;

    end
end

% Plot adsorption rates
errorbar(10.^(-dT), mean(eta(:, :), 2), std(eta(:, :), [], 2) / sqrt(repeats), 'kx', 'LineWidth', 1.5, 'MarkerFaceColor', 'Auto', 'DisplayName', sprintf('n = %d', repeats))

% Label the figure
xlabel('{\Delta}T (h)')
ylabel('\eta ({\mu}m^3/h)')
ax.XScale = 'log';
ax.LineWidth = 1.5;
ax.FontSize = 16;

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
print(fh, '../figures/Figure_S7/FigS7.tif', '-dtiff', '-r900')