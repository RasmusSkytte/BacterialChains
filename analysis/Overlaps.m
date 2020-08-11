% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_S1', 'dir')
    mkdir ../figures/Figure_S1
end

% Define the simulations
N = [3 10 32 100 316 1000];
theta = 0.2:0.2:1;

repeats = 100;

% Make arrays to store the fitted values in
max_overlaps  = nan(numel(N), repeats, numel(theta) + 2);
median_overlaps = nan(numel(N), repeats, numel(theta) + 2);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for i = 1:numel(N)
    for t = 1:(numel(theta)+2)
        for j = 1:repeats

            % Define the path for the data
            if t == numel(theta) + 1
                dpath = sprintf('%s/WellMixed/N_%d', path, N(i));
            elseif t == numel(theta) + 2
                dpath = sprintf('%s/SphericalColony/N_%d', path, N(i));
            else
                dpath = sprintf('%s/Chain/N_%d/theta_%.3f_pi', path, N(i), theta(t));
            end

            % Check if data exists
            if exist(sprintf('%s/repeat_%d/Overlaps.txt', dpath, j-1), 'file')
                overlap = importdata(sprintf('%s/repeat_%d/Overlaps.txt', dpath, j-1));
                flog = importdata(sprintf('%s/repeat_%d/log.txt', dpath, j-1), '=');
            else
                continue
            end

            % Store the largest and median overlap
            max_overlaps(i, j, t)  = max(overlap);
            median_overlaps(i, j, t) = median(overlap);
        end
    end
end

% Store the radius of the cells
R = flog.data(10);

% Test for overlaps
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

overlap = squeeze(nanmean(max_overlaps, 2))' / R;
d_overlap = squeeze(nanstd(max_overlaps, [], 2))' ./ (squeeze(sum(~isnan(max_overlaps), 2))' * R);

for t = 1:numel(theta) + 2
    if t <= numel(theta)
        name = sprintf('\\Theta / 2 = %.1f \\pi', theta(t) / 2);
        symbol = 's';
    elseif t == numel(theta) + 1
        name = 'Well-Mixed';
        symbol = '^';
    elseif t == numel(theta) + 2
        name = 'Spherical Colony';
        symbol = 'o';
    end
    e = errorbar(ax, N(1:end), overlap(t, :), d_overlap(t, :), symbol, 'LineWidth', 1.5, 'DisplayName', name, 'MarkerFaceColor', 'auto');
end
, '-dtiff', '-r900'
xlabel(ax, 'N')
ylabel(ax, 'max(overlap) / R')

set(ax, 'XScale', 'log');

ax.YTick = 0:0.01:0.02;

ax.YLim = [0 0.02];
ax.XLim = [1 1000];

ax.LineWidth = 1.5;
ax.FontSize = 16;

l = legend('Location', 'NorthEastOutside');
ax.Position = [0.12 0.17 0.55 0.8];

pause(0.1); fh.Position = [10 50 720 420]; l.Position(1) = 0.7; pause(0.1);
print(fh, '../figures/Figure_S1/FigS1.tif', '-dtiff', '-r900')