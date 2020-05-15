% Clean up
close all;
clearvars;

if ~exist('../figures/Figure_7', 'dir')
    mkdir ../figures/Figure_7
end

if ~exist('../figures/Figure_8', 'dir')
    mkdir ../figures/Figure_8
end

% Define the simulations
N = [3 10 32 100 316 1000];
theta = 0.2:0.2:1;

samples     = 5;
beta        = 1e2;
nBootstraps = 1e3;

% Set the seed
rng(0)

% Make arrays to store the test
m_Burstsize    = zeros(numel(N), numel(theta) + 2);
s_Burstsize    = zeros(numel(N), numel(theta) + 2);
m_Escape       = zeros(numel(N), numel(theta) + 2);
s_Escape       = zeros(numel(N), numel(theta) + 2);
m_MultipleHits = zeros(numel(N), numel(theta) + 2);
s_MultipleHits = zeros(numel(N), numel(theta) + 2);

Adsorption = nan(numel(theta) + 2, samples, 10001);
weights = nan(numel(theta) + 2, 5);

% Prepare path
path = strrep(pwd, 'analysis', 'cpp/data');

% Loop over runs
for n = 1:numel(N)

    for t = 1:(numel(theta) + 2)

        % Define the path for the data
        if t == numel(theta) + 1
            dpath = sprintf('%s/WellMixed/N_%d/repeat_0', path, N(n));
        elseif t == numel(theta) + 2
            dpath = sprintf('%s/SphericalColony/N_%d/repeat_0', path, N(n));
        else
            dpath = sprintf('%s/Chain/N_%d/theta_%.3f_pi/repeat_0', path, N(n), theta(t));
        end

        % Check if reference data exists
        if exist(sprintf('%s/PhageLocation.txt', dpath), 'file')
            rdata = importdata(sprintf('%s/PhageLocation.txt', dpath));
        else
            continue
        end

        % Compute the ID of lysed cells
        I = floor((0:(min(N(n), samples)-1)) * (N(n) - 1) / (min(N(n), samples) - 1)) + 1;

        % Get hit probabilities from reference data.
        w = histcounts(rdata(:, end), 0.5:(N(n) + 0.5));
        w = w(I)/sum(w(I));


        for j = 1:min(N(n), samples)

            % Check if data exists
            if exist(sprintf('%s/lysis_%d/PhageLocation.txt', dpath, j-1), 'file')
                data = importdata(sprintf('%s/lysis_%d/PhageLocation.txt', dpath, j-1));
            else
                continue
            end

            % Check if data exists
            if exist(sprintf('%s/lysis_%d/ColonySize.txt', dpath, j-1), 'file')
                sdata = importdata(sprintf('%s/lysis_%d/ColonySize.txt', dpath, j-1));
            else
                continue
            end

            % Store the adsorption fraction
            if N(n) == 10
                Adsorption(t, j, :) = sdata(:, 2);
                weights(t, :) = w;
            end

            % Perform bootstrap loop
            Burstsizes   = nan(1, nBootstraps);
            Escapes      = nan(1, nBootstraps);
            MultipleHits = nan(1, nBootstraps);

            for i = 1:nBootstraps

                % Choose beta random phages (without replacement)
                r = randi(size(data, 1), beta, 1);

                % Plot examples
                if ( (i == 1) && (t == 2 || t == numel(theta) || t == numel(theta) + 2 ) && (N(n) == 10) && (j == 3))

                    % Load full phage data
                    if exist(sprintf('%s/lysis_%d/PhageData.txt', dpath, j-1), 'file')
                        pdata = importdata(sprintf('%s/lysis_%d/PhageData.txt', dpath, j-1));
                    else
                        sprintf('%s/lysis_%d/PhageData.txt', dpath, j-1)
                        error('!');
                    end

                    % Load log
                    if exist(sprintf('%s/lysis_%d/log.txt', dpath, j-1), 'file')
                        log = importdata(sprintf('%s/lysis_%d/log.txt', dpath, j-1));
                        R = log.data(10);
                    else
                        error('!');
                    end

                    % Load full cell data
                    if exist(sprintf('%s/CellData.txt', dpath), 'file')
                        cdata = importdata(sprintf('%s/CellData.txt', dpath));
                    else
                        error('!');
                    end

                    % Load log
                    if exist(sprintf('%s/log.txt', dpath), 'file')
                        clog = importdata(sprintf('%s/log.txt', dpath));
                    else
                        error('!');
                    end

                    % Store the ID of the cells
                    ID = [1:(I(j)-1) nan I(j):(N(n)-1)];

                    fh = figure();
                    fh.Resize = 'off';
                    ax = axes;
                    ax.NextPlot = 'add';

                    % Plot the bacteria as reference
                    for l = 1:N(n)
                        x1 = cdata(l, 2) - clog.data(5) / 2;
                        y1 = cdata(l, 3) - clog.data(5) / 2;
                        z1 = cdata(l, 4) - clog.data(5) / 2;
                        x2 = cdata(l, 5) - clog.data(5) / 2;
                        y2 = cdata(l, 6) - clog.data(5) / 2;
                        z2 = cdata(l, 7) - clog.data(5) / 2;
                        L  = norm([x2 - x1, y2 - y1, z2 - z1]);

                        % Generate cell shape
                        x = linspace(-L/2 - R, L/2 + R, 50);
                        d = ones(size(x)) * R;
                        d(abs(x) > L/2) = sqrt(abs(R^2 - (abs(x(abs(x) > L/2)) - L/2).^2));
                        [X,Y,Z] = cylinder(d, 50);

                        % Rescale the z axis to match cell length
                        Z = (L + 2 * R) * (Z - 0.5);

                        % Rotate cell shape to match angle of cell and center it at right
                        % coordinates
                        z0 = [x2 - x1, y2 - y1, z2 - z1] / L;

                        a = cross([0 0 1], z0);
                        a = a / norm(a);
                        x = a(1); y = a(2); z = a(3);
                        c = [0 0 1] * z0';              % cos(theta)
                        d = 1-c;                        % 1 - cos(theta)
                        s = sqrt(1-c^2);                % sin(theta)

                        Rot =  [x^2*d+c     x*y*d-z*s   x*z*d+y*s;
                                x*y*d+z*s   y^2*d+c     y*z*d-x*s;
                                x*z*d-y*s   y*z*d+x*s   z^2*d+c  ];

                        for x = 1:size(X,1)
                            for y = 1:size(X,2)
                                v = Rot*[X(x,y) Y(x,y) Z(x,y)]';

                                X(x,y) = v(1) + (x1 + x2) / 2;
                                Y(x,y) = v(2) + (y1 + y2) / 2;
                                Z(x,y) = v(3) + (z1 + z2) / 2;

                                if imag(X(x,y)) > 0
                                    disp('!')
                                end
                                if imag(Y(x,y)) > 0
                                    disp('!')
                                end
                                if imag(Z(x,y)) > 0
                                    disp('!')
                                end
                            end
                        end

                        % Plot the cell
                        if l == I(j)
                            surf(ax, X, Y, Z, ones(size(X)), 'FaceColor', [0.60 0.60 0.60], 'EdgeColor', 'None', 'FaceAlpha', 0.5);
                        else
                            surf(ax, X, Y, Z, ones(size(X)), 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'None', 'FaceAlpha', 0.5);
                        end

                    end

                    x0 = log.data(5) / 2;
                    scatter3(ax, pdata(r(data(r) > 0), 2) - x0, pdata(r(data(r) > 0), 3) - x0, pdata(r(data(r) > 0), 4) - x0, '.r', 'SizeData', 50);
                    axis equal;
                    ax.Visible = 'off';
                    if t <= numel(theta) && theta(t) == 0.2
                        view([40 90])
                        ax.Position = [-0.35 -0.25 1.5 1.5];
                    elseif t <= numel(theta) && theta(t) == 0.4
                        view([140 0])
                        ax.Position = [-0.1 -0.15 1.3 1.3];
                    elseif t == numel(theta)
                        fh.Position(4) = 280;
                        view([1 1 1])
                        ax.Position = [0 0 1.1 1.1];
                    elseif t == numel(theta) + 2
                        fh.Position(4) = 140;
                        view([12 60])
                        ax.Position = [-0.125 -0.15 1.25 1.25];
                    end

                    pause(0.1); fh.Position(1:3) = [10 50 280]; pause(0.1);
                    if t == 2
                        saveas(fh, '../figures/Figure_7/Fig7b.png')
                    elseif t == numel(theta)
                        saveas(fh, '../figures/Figure_7/Fig7c.png')
                    elseif t == numel(theta) + 2
                        saveas(fh, '../figures/Figure_7/Fig7d.png')
                    end
                end

                % Get beta random phage fates
                hits = data(r);

                % Get the distribution of hits
                hitDist = histcounts(hits, unique([-1.5 -0.5 unique(hits)'+0.5]));
                nInfected = numel(hitDist) - 1;

                % Count the effective burst size
                Burstsizes(i)   = (hitDist(1) + nInfected) / beta;

                % Count fraction of hits that is multiple hits
                MultipleHits(i) = sum(hitDist(2:end) > 1) / nInfected;

                % Count the fraction of phages that escapes
                Escapes(i)   = hitDist(1) / beta;
            end

            % Store the average and std
            m_Burstsize(n, t)    = m_Burstsize(n, t)    + w(j) * mean(Burstsizes);
            s_Burstsize(n, t)    = s_Burstsize(n, t)    + w(j) * std(Burstsizes);

            m_Escape(n, t)       = m_Escape(n, t)       + w(j) * mean(Burstsizes);
            s_Escape(n, t)       = s_Escape(n, t)       + w(j) * std(Burstsizes);

            m_MultipleHits(n, t) = m_MultipleHits(n, t) + w(j) * nanmean(MultipleHits);
            s_MultipleHits(n, t) = s_MultipleHits(n, t) + w(j) * nanstd(MultipleHits);

        end
    end
end

% Plot the adsorption for N = 100
fh = figure;
fh.Resize = 'off';
ax = axes;
ax.Box = 'on';
ax.NextPlot = 'add';

cc = lines(numel(theta) + 2);

T = sdata(:, 1) - sdata(1, 1);
%     Adsorption = Adsorption(:, 1:3, :);
for t = 1:(numel(theta)+2)
    wm = sum(  squeeze(Adsorption(t, :, :))                                           .* repmat(weights(t, :)', 1, size(Adsorption, 3)));
    ws = sum( (squeeze(Adsorption(t, :, :)) - repmat(wm, numel(weights(t, :)), 1)).^2 .* repmat(weights(t, :)', 1, size(Adsorption, 3))) * numel(weights(t, :)) / (numel(weights(t, :)) - 1);
    ws = sqrt(ws) / sqrt(samples);

    if t <= numel(theta)
        h = errorbar(ax, T * 3600, wm, ws, 's', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'Auto');
    elseif t == numel(theta) + 1
        h = errorbar(ax, T * 3600, wm, ws, '^', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'Auto');
    elseif t == numel(theta) + 2
        h = errorbar(ax, T * 3600, wm, ws, 'o', 'LineWidth', 1.5, 'Color', cc(t, :), 'MarkerFaceColor', 'Auto');
    end

    axes(ax);
    hs = shadedErrorBar(h.XData, h.YData, [h.YPositiveDelta; h.YNegativeDelta], 'lineProps', {'Color', h.Color, 'LineWidth', 1.5});
    delete(h);

end

ylabel(ax, 'Free phage particles')
xlabel(ax, 'Time (s)')

ax.LineWidth = 1.5;
ax.FontSize = 16;

ax.XLim = [T(1) 10];
ax.YLim = [1000 wm(1)];

ax.YScale = 'log';

xtickformat(ax, '%.1f')

annotation('textarrow', [0.6 0.55], [0.88 0.92], 'String', 'Well-mixed', 'FontSize', 14);
annotation('textarrow', [0.7 0.6],  [0.52 0.67], 'String', 'Chains', 'FontSize', 14);
annotation('textarrow', [0.4 0.45], [0.23 0.28], 'String', 'Spherical Colony', 'FontSize', 14);

pause(0.1); fh.Position = [10 50 560 420]; pause(0.1);
saveas(fh, '../figures/Figure_8/Fig8b.png')

% Plot the effective burstsize
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

ax1.Position = [0.16  0.15 0.1    0.775];
ax2.Position = [0.325 0.15 0.4875 0.775];
ax3.Position = [0.87  0.15 0.1    0.775];

cc = lines(numel(theta) + 2); cc(6, :) = cc(3, :); cc(4, :) = cc(2, :); cc(2, :) = cc(1, :);

for n = [2 3 4] %1:numel(N)

    for t = 1:(numel(theta) + 2)

        if t <= numel(theta)
            errorbar(ax2, theta(t) / 2, m_Burstsize(n, t), s_Burstsize(n, t) / sqrt(min(N(n), samples)), 's', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        elseif t == numel(theta) + 1
            errorbar(ax1, 0.5,          m_Burstsize(n, t), s_Burstsize(n, t) / sqrt(min(N(n), samples)), '^', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        elseif t == numel(theta) + 2
            errorbar(ax3, 0.5,          m_Burstsize(n, t), s_Burstsize(n, t) / sqrt(min(N(n), samples)), 'o', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        end

    end

end

ylabel(ax1, '\beta_{eff} / \beta')
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
saveas(fh, '../figures/Figure_7/Fig7a.png')


% Plot the multiple hits
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

ax1.Position = [0.16  0.15 0.1    0.775];
ax2.Position = [0.325 0.15 0.4875 0.775];
ax3.Position = [0.87  0.15 0.1    0.775];

cc = lines(numel(theta) + 2); cc(6, :) = cc(3, :); cc(4, :) = cc(2, :); cc(2, :) = cc(1, :);
for n = [2 3 4] %1:numel(N)

    for t = 1:(numel(theta) + 2)

        if t <= numel(theta)
            errorbar(ax2, theta(t) / 2, m_MultipleHits(n, t), s_MultipleHits(n, t) / sqrt(min(N(n), samples)), 's', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        elseif t == numel(theta) + 1
            errorbar(ax1, 0.5,          m_MultipleHits(n, t), s_MultipleHits(n, t) / sqrt(min(N(n), samples)), '^', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        elseif t == numel(theta) + 2
            errorbar(ax3, 0.5,          m_MultipleHits(n, t), s_MultipleHits(n, t) / sqrt(min(N(n), samples)), 'o', 'LineWidth', 1.5, 'Color', cc(n, :), 'MarkerFaceColor', 'Auto')
        end

    end

end

ylabel(ax1, '# Super-infected / # Infected')
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
saveas(fh, '../figures/Figure_8/Fig8a.png')

% Store the effective burst size
save('EffectiveBurstSize.mat', 'm_Escape', 's_Escape')