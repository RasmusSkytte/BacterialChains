% Clean up
close all;
clearvars;

D = 1e4;
xint = [1e-2 1e4];

if ~exist('../figures/Figure_5', 'dir')
    mkdir ../figures/Figure_5
end

if ~exist('../figures/Figure_S4', 'dir')
    mkdir ../figures/Figure_S4
end

if ~exist('../figures/Figure_S5', 'dir')
    mkdir ../figures/Figure_S5
end

if ~exist('../figures/Figure_S10', 'dir')
    mkdir ../figures/Figure_S10
end

% Define the simulations
N = [3 10 32 100 316 1000];
theta = 0.2:0.2:1;

repeats = 100;

% Define colors
cc = lines(numel(theta));

% Make arrays to store the persistence length
P     = nan(numel(N), repeats, numel(theta));
d_P   = nan(numel(N), repeats, numel(theta));

Rg    = nan(numel(N), repeats, numel(theta));


cc = lines(numel(theta) + 2);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for i = 1:numel(N)
    for t = 1:numel(theta)
        for j = 1:repeats

            % Define the path for the data
            dpath = sprintf('%s/Chain/N_%d/theta_%.3f_pi', path, N(i), theta(t));

            % Check if data exists
            if exist(sprintf('%s/repeat_%d/PersistenceLength.txt', dpath, j-1), 'file')
                data = importdata(sprintf('%s/repeat_%d/PersistenceLength.txt', dpath, j-1));

                % Store the fitted values
                P(i, j, t)   = data(1, 1);
                d_P(i, j, t) = data(1, 2);
            end

            % Check if data exists
            if exist(sprintf('%s/repeat_%d/GyrationRadius.txt', dpath, j-1), 'file')
                Rg(i, j, t) = importdata(sprintf('%s/repeat_%d/GyrationRadius.txt', dpath, j-1));
            end

            % Show example of persistence length fits
            if (N(i) == 32 && j == 1 && any(round(theta(t), 1) == [0.2 0.6]))

                % Load cell coordinates
                data = importdata(sprintf('%s/repeat_%d/CellData.txt', dpath, j-1));

                % Convert to chain coordiantes
                xx = reshape(data(:, 2:7)', 3, 2*size(data,1))';

                % Compute arc length along the chain
                arclen = [0; cumsum(sqrt(sum((xx(1:(end-1), :) - xx(2:end, :)).^2, 2)))];

                distance = [];
                R2 = [];
                for m = 1:(numel(arclen)-1)
                    R2 = [R2; sum((xx((m+1):end, :) - xx(m, :)).^2, 2)];
                    distance = [distance; arclen((m+1):end, :) - arclen(m)];
                end

                fh = figure;
                fh.Resize = 'off';
                hold on; box on;
                ax = gca;
                ax.Position = [0.20 0.17 0.77 0.8];

                l = linspace(min(distance), max(distance));

                scatter(distance, R2, '.', 'MarkerEdgeColor', cc(t, :))
                plot(l, 2 * P(i, j, t) * l .* (1 - P(i, j, t) ./ l .* (1 - exp(- l / P(i, j, t)))), 'k', 'LineWidth', 2)


                ax.Position(4) = 0.76;

                xlabel('l ({\mu}m)')
                ylabel('R^2 ({\mu}m^2)')
                ax.LineWidth = 1.5;
                ax.FontSize = 16;

                ytickformat(ax, '%.1f')

                pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
                if round(theta(t), 1) == 0.2
                    saveas(fh, '../figures/Figure_S4/FigS4a.png')
                elseif round(theta(t), 1) == 0.6
                    saveas(fh, '../figures/Figure_S4/FigS4b.png')
                end

            end
        end
    end
end

% Take weighted mean over persistence length measurements
%     w = 1./d_P.^2;
%     P = squeeze(sum(w.*P, 2))' ./ squeeze(sum(w, 2))';
%     d_P = sqrt(1./squeeze(sum(w, 2)))';

d_P = squeeze(std(P, [], 2))' / sqrt(repeats);
P   = squeeze(mean(P, 2))';

% Determine SEM and mean of gyration radius
d_Rg = squeeze(std(Rg, [], 2))' ./ sqrt(repeats);
Rg   = squeeze(mean(Rg, 2))';


% Determine differnces in relative error
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

scatter(d_Rg(:) ./ Rg(:), d_P(:) ./ P(:), 'filled')

xlabel(ax, 'Relative error in R_g')
ylabel(ax, 'Relative error in P')

ax.LineWidth = 1.5;
ax.FontSize = 16;


% Measure in cell numbers
ell = 1.5+0.9;
P   =   P ./ ell;
d_P = d_P ./ ell;

% Compute adjusted length
NN = repmat(N, numel(theta), 1);
N_adj = NN ./ P;

% Keep track of those more longer than P, and those longer than 10 P
I1  = N_adj > 1;
I10 = N_adj > 10;

% Load measured adsorption rates
load(['AdsorptionRate.mat']);

% Extract the relevant adsorption rates
eta   =   eta(1:numel(theta), 1+(1:numel(N)));
d_eta = d_eta(1:numel(theta), 1+(1:numel(N)));

% Determine differnces in relative error
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

scatter(d_eta(:) ./ eta(:), d_P(:) ./ P(:), 'filled')

xlabel(ax, 'Relative error in \eta')
ylabel(ax, 'Relative error in P')

ax.LineWidth = 1.5;
ax.FontSize = 16;


% Plot R_g as a function of N
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

for i = 1:numel(theta)
    errorbar(ax, NN(i, :), Rg(i, :), d_Rg(i, :), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
end
xlim([1 N(end)])

ax.Position(4) = 0.76;

ax.XScale = 'log';
ax.YScale = 'log';

xlabel(ax, 'N')
ylabel(ax, 'R_g ({\mu}m)')
ax.LineWidth = 1.5;
ax.FontSize = 16;

pause(0.1); fh.Position = [10 50 480 420]; pause(0.1);
saveas(fh, '../figures/Figure_S5/FigS5a.png')


% Plot R_g' as a function of N'
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% Compute x and y errors for error bars
% x = N / P
x = N_adj;

% dx comes only from P
dx = (d_P ./ P) .* x;

% reference values is R_g / P
y  = Rg ./ P;

% dy comes from P and Rg
dy = sqrt((d_Rg ./ Rg).^2 + (d_P ./ P).^2) .* y;

for i = 1:numel(theta)
    errorbar(ax, x(i, :), y(i, :), dy(i, :), dy(i, :), dx(i, :), dx(i, :), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
end
xlim([1 N(end)])

ax.Position(4) = 0.76;

ax.XScale = 'log';
ax.YScale = 'log';

xlabel(ax, 'N^\prime')
ylabel(ax, 'R_g^\prime ({\mu}m)')
ax.LineWidth = 1.5;
ax.FontSize = 16;

pause(0.1); fh.Position = [10 50 480 420]; pause(0.1);
saveas(fh, '../figures/Figure_S5/FigS5b.png')



% Plot R_g / P vs fitted relation
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% Compute x and y errors for error bars
% x = N / P
x = N_adj;

% dx comes only from P
dx = (d_P ./ P) .* x;

% reference values is R_g / P
y  = Rg ./ P;

% dy comes from P and Rg
dy = sqrt((d_Rg ./ Rg).^2 + (d_P ./ P).^2) .* y;

% Use our fit on the long chains
xx  =  x(I10(:));
dxx = dx(I10(:));

yy  =  y(I10(:));
dyy = dy(I10(:));

% Define fit function
fun   = @(x, p) p(1) .* x.^p(2);
% and the error on the fit function
d_fun = @(x, dx, p) abs( (p(1) * p(2) .* x.^(p(2)-1)) .* dx);

% Extract the best parameters set p
fprintf('-------------------------------\n')
fprintf('--- R_g / P vs N / P (long) ---\n')
[p_Rg, ~] = fit_ratio(yy, dyy, @(p)fun(xx, p), @(p)d_fun(xx, dxx, p), [0.5 0.5]);

% Compute the ratio
ratio = y ./ fun(x, p_Rg);

% Compute the error on the ratio
d_ratio = sqrt((dy ./ y).^2 + (d_fun(x, dx, p_Rg) ./ fun(x, p_Rg)).^2) .* ratio;

% Plot reference line
plot(xint, [1 1], 'k', 'LineWidth', 1);

for i = 1:numel(theta)
    errorbar(ax, x(i,  I10(i, :)), ratio(i,  I10(i, :)), d_ratio(i,  I10(i, :)), d_ratio(i,  I10(i, :)), dx(i,  I10(i, :)), dx(i,  I10(i, :)), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
    errorbar(ax, x(i, ~I10(i, :)), ratio(i, ~I10(i, :)), d_ratio(i, ~I10(i, :)), d_ratio(i, ~I10(i, :)), dx(i, ~I10(i, :)), dx(i, ~I10(i, :)), '^', 'MarkerFaceColor', 'None',   'Color', cc(i, :), 'LineWidth', 1.5);
end
ax.Position(4) = 0.76;

ax.XScale = 'log';

xlabel(ax, 'N^\prime')
ylabel(ax, 'R_g^\prime / k(N^\prime)^\zeta')
ax.LineWidth = 1.5;
ax.FontSize = 16;

pause(0.1); fh.Position = [10 50 480 420]; pause(0.1);
saveas(fh, '../figures/Figure_S5/FigS5c.png')




% Plot R_g / P vs fitted relation (small N_adj)
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% Compute x and y errors for error bars
% x = N / P
x = N_adj;

% dx comes only from P
dx = (d_P ./ P) .* x;

% reference values is R_g / P
y  = Rg ./ P;

% dy comes from P and Rg
dy = sqrt((d_Rg ./ Rg).^2 + (d_P ./ P).^2) .* y;

% Use our fit on the short chains
xx  =  x(~I10(:));
dxx = dx(~I10(:));

yy  =  y(~I10(:));
dyy = dy(~I10(:));

% Define fit function
fun   = @(x, p) p(1) .* x.^p(2);
% and the error on the fit function
d_fun = @(x, dx, p) abs( (p(1) * p(2) .* x.^(p(2)-1)) .* dx);

% Extract the best parameters set p
fprintf('-------------------------------\n');
fprintf('--- R_g / P vs N / P (short) --\n')
[p, ~] = fit_ratio(yy, dyy, @(p)fun(xx, p), @(p)d_fun(xx, dxx, p), [0.5 0.5]);

% Compute the ratio
ratio = y ./ fun(x, p);

% Compute the error on the ratio
d_ratio = sqrt((dy ./ y).^2 + (d_fun(x, dx, p) ./ fun(x, p)).^2) .* ratio;

% Plot reference line
plot(xint, [1 1], 'k', 'LineWidth', 1);

for i = 1:numel(theta)
    errorbar(ax, x(i, ~I10(i, :)), ratio(i, ~I10(i, :)), d_ratio(i, ~I10(i, :)), d_ratio(i, ~I10(i, :)), dx(i, ~I10(i, :)), dx(i, ~I10(i, :)), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
    errorbar(ax, x(i,  I10(i, :)), ratio(i,  I10(i, :)), d_ratio(i,  I10(i, :)), d_ratio(i,  I10(i, :)), dx(i,  I10(i, :)), dx(i,  I10(i, :)), '^', 'MarkerFaceColor', 'None',   'Color', cc(i, :), 'LineWidth', 1.5);
end
ax.Position(4) = 0.76;

ax.XScale = 'log';

xlabel(ax, 'N^\prime')
ylabel(ax, 'R_g^\prime / k(N^\prime)^\zeta')
ax.LineWidth = 1.5;
ax.FontSize = 16;



% Plot eta / P vs N / P
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% Compute x and y errors for error bars
% x = N / P
x = N_adj;

% dx comes only from P
dx = (d_P ./ P) .* x;

% reference values is eta / P
y  = eta ./ P;

% dy comes from P and eta
dy = sqrt((d_eta ./ eta).^2 + (d_P ./ P).^2) .* y;

% Use our fit on the long chains
xx  =  x(I10(:));
dxx = dx(I10(:));

yy  =  y(I10(:));
dyy = dy(I10(:));

% Define fit function
fun   = @(x, p) p(1) .* x.^p(2);
% and the error on the fit function
d_fun = @(x, dx, p) abs( (p(1) * p(2) .* x.^(p(2)-1)) .* dx);

% Extract the best parameters set p
fprintf('-------------------------------\n');
fprintf('--- eta / P vs N / P (long) ---\n')
[p3D, ~] = fit_ratio(yy, dyy, @(p)fun(xx, p), @(p)d_fun(xx, dxx, p), [1e5 0.5]);

% Use our fit on the short chains
xx  =  x(~I10(:));
dxx = dx(~I10(:));

yy  =  y(~I10(:));
dyy = dy(~I10(:));

% Define fit function
fun   = @(x, p) p(1) .* x.^p(2);
% and the error on the fit function
d_fun = @(x, dx, p) abs( (p(1) * p(2) .* x.^(p(2)-1)) .* dx);

% Extract the best parameters set p
fprintf('-------------------------------\n');
fprintf('--- eta / P vs N / P (short) --\n')
[p2D, ~] = fit_ratio(yy, dyy, @(p)fun(xx, p), @(p)d_fun(xx, dxx, p), [5e4 0.75]);

% Plot reference line
plot(ax, xint, fun(xint, p3D), 'k',   'LineWidth', 1.5, 'DisplayName', 'Gaussian Chain');
plot(ax, xint, fun(xint, p2D), 'k--', 'LineWidth', 1.5, 'DisplayName', 'N < 3 P');

cc = lines(numel(theta));
for i = 1:numel(theta)
    errorbar(ax, x(i,  I10(i, :)), y(i,  I10(i, :)), dy(i,  I10(i, :)), dy(i,  I10(i, :)), dx(i,  I10(i, :)), dx(i,  I10(i, :)), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
    errorbar(ax, x(i, ~I10(i, :)), y(i, ~I10(i, :)), dy(i, ~I10(i, :)), dy(i, ~I10(i, :)), dx(i, ~I10(i, :)), dx(i, ~I10(i, :)), '^', 'MarkerFaceColor', 'None',   'Color', cc(i, :), 'LineWidth', 1.5);
end

ax.XScale = 'log';
ax.YScale = 'log';

xlabel(ax, 'N^\prime')
ylabel(ax, '\eta^\prime({\mu}m^3/h)')
ax.LineWidth = 1.5;
ax.FontSize = 16;


% Plot von Smoluchowski prediction based on R_g
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% x = N / P
x = N_adj;

% dx comes only from P
dx = (d_P ./ P) .* x;

% y values is eta / (4 * pi * D * Rg)
y  = eta ./ (4 * pi * D * Rg);

% dy comes from eta and Rg
dy = sqrt((d_eta ./ eta).^2 + (d_P ./ P).^2) .* y;

% Plot reference line
plot(xint, [1 1], 'k', 'LineWidth', 1);

for i = 1:numel(theta)
    errorbar(ax, x(i, :), y(i, :), dy(i, :), dy(i, :), dx(i, :), dx(i, :), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
end

ax.Position(4) = 0.76;

ax.XScale = 'log';
ax.YScale = 'lin';

xlabel(ax, 'N^\prime')
ylabel(ax, '\eta / 4 \pi D R_g')
ax.LineWidth = 1.5;
ax.FontSize = 16;

ytickformat(ax, '%.1f')

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
saveas(fh, '../figures/Figure_S10/FigS10.png')




% Plot our fitted prediction for eta
fh = figure;
fh.Resize = 'off';
hold on; box on;
ax = gca;
ax.Position = [0.17 0.17 0.79 0.8];

% Compute x and y errors for error bars
% x = P
x  = P;
dx = d_P;

% reference values is eta
y  = eta;

% dy comes only from eta
dy = d_eta;

% Use our fit on the long chains
xx  =  x(I10(:));
dxx = dx(I10(:));
n   = NN(I10(:));

yy  =  y(I10(:));
dyy = dy(I10(:));

% Define fit function
fun   = @(x, n, p) 4 * pi * D * x .* (n ./ x).^p;

% and the error on the fit function
d_fun = @(x, n, dx, p) 4 * pi * D * (1-p) * (n ./ x).^p .* dx;

% Extract the best parameters set p
fprintf('-------------------------------\n');
fprintf('-- P(N / P)^g vs N / P (long) -\n')
[p, ~] = fit_ratio(yy, dyy, @(p)fun(xx, n, p), @(p)d_fun(xx, n, dxx, p), 0.5);

% Compute the ratio
ratio = y ./ fun(x, NN, p);

% Change x to N_adj for plotting
x = N_adj;

% dx comes only from P
dx = (d_P ./ P) .* x;

% Compute the error on the ratio
d_ratio = sqrt((dy ./ y).^2 + (d_fun(x, NN, dx, p) ./ fun(x, NN, p)).^2) .* ratio;

% Plot reference line
plot(xint, [1 1], 'k', 'LineWidth', 1);

for i = 1:numel(theta)
    errorbar(ax, x(i,  I10(i, :)), ratio(i,  I10(i, :)), d_ratio(i,  I10(i, :)), d_ratio(i,  I10(i, :)), dx(i,  I10(i, :)), dx(i,  I10(i, :)), 's', 'MarkerFaceColor', cc(i, :), 'Color', cc(i, :), 'LineWidth', 1.5);
    errorbar(ax, x(i, ~I10(i, :)), ratio(i, ~I10(i, :)), d_ratio(i, ~I10(i, :)), d_ratio(i, ~I10(i, :)), dx(i, ~I10(i, :)), dx(i, ~I10(i, :)), '^', 'MarkerFaceColor', 'None',   'Color', cc(i, :), 'LineWidth', 1.5);
end

ax.XScale = 'log';

xlabel(ax, 'N^\prime')
ylabel(ax, '\eta  / 4 \pi D P (N^\prime)^{\gamma_P}')
ax.LineWidth = 1.5;
ax.FontSize = 16;

pause(0.1); fh.Position = [10 50 480 420]; pause(0.1);
saveas(fh, '../figures/Figure_5/Fig5c.png')


% Store the radius of gyration
save('RadiusOfGyration.mat', 'Rg', 'd_Rg')