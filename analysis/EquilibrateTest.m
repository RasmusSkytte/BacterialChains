% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_S9', 'dir')
    mkdir ../figures/Figure_S9
end

T_Eq = [0 1 2 3 4 5 6];
N    = [1 32];

% Time-steps to check
repeats = 3;

% Make arrays to store the fitted values in
eta   = nan(numel(T_Eq), repeats);
d_eta = nan(numel(T_Eq), repeats);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');ls
ls

% Loop over runs
for n = 1:numel(N)

    fh = figure;
    fh.Resize = 'off';
    hold on; box on;
    ax = gca;

    for t = 1:numel(T_Eq)
        for r = 1:repeats

            % Load data from the run
            data = importdata(sprintf('%s/EquilibriumTest/T_%d/N_%d/repeat_%d/ColonySize.txt', path, T_Eq(t), N(n), r-1));

            % Get time and phage numbers
            T = data(:, 1);
            P = data(:, 2);

            % Open the log to get the meta data
            flog = importdata(sprintf('%s/EquilibriumTest/T_%d/N_%d/repeat_%d/log.txt', path, T_Eq(t), N(n), r-1));
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
    errorbar(T_Eq, mean(eta, 2), std(eta, [], 2) / sqrt(repeats), 'kx', 'LineWidth', 1.5, 'MarkerFaceColor', 'Auto', 'DisplayName', sprintf('n = %d', repeats))

    % Label the figure
    xlabel('T_{Eq} (h)')
    ylabel('\eta ({\mu}m^3/h)')
    ax.LineWidth = 1.5;
    ax.FontSize = 16;
    ax.Position(1) = 0.2;
    ax.XTick = 0:6;

    ytickformat(ax, '%.2f');
    ax.YLabel.Position(1) = 1.2*ax.YLabel.Position(1);

    pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
    if n == 1
        print(fh, '../figures/Figure_S9/FigS9a.tif', '-dtiff', '-r900')
    elseif n == 2
        print(fh, '../figures/Figure_S9/FigS9b.tif', '-dtiff', '-r900')
    end

end