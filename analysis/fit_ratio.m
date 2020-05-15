function [p, dp] = fit_ratio(y, dy, f, df, p0)

    % Minimize the function abs(y / f(x, p) - 1)
    
    % Define squared difference
    squared_difference = @(p) (y ./ f(p) - 1).^2;
    
    % Define the error of the ratio y / f(x, p)
    squared_errors = @(p)  ( (dy ./ y).^2 + (df(p) ./ f(p)).^2 ) .* (y ./ f(p)).^2;
    
    % Define chi^2
    chi2 = @(p) sum(squared_difference(p) ./ squared_errors(p));
    
    % Perform minimization
    options = optimoptions('fmincon', 'Display', 'off');
    p = fmincon(chi2, p0, [], [], [], [], zeros(1, numel(p0)), inf(1, numel(p0)), [], options);
    
    % Estimate the error on p
    min_chi2 = chi2(p);
    
    % Allocate array for p errors
    dp = nan(2, numel(p));
    
    % Locate where chi^2 = min(chi^2) + 1
    for i = 1:numel(p)
       dp(1, i) =  fmincon(@(delta)(chi2(p + delta*circshift(1:numel(p) == 1, i-1)) - (min_chi2 + 1)), 0, [], [], [], [], 0,   Inf, [], options);
       dp(2, i) = -fmincon(@(delta)(chi2(p + delta*circshift(1:numel(p) == 1, i-1)) - (min_chi2 + 1)), 0, [], [], [], [], -p(i), 0, [], options);
    end
    
    % Take the average deviation as the error
    dp = mean(dp);
    
    % Report fit result
    fprintf('-------------------------------\n')
    fprintf('---------- Ratio fit-----------\n')
    fprintf('chi^2 = %.3g, ndf = %d\n', chi2(p), numel(y)-numel(p))
    fprintf('chi^2 / ndf = %.3g\n', chi2(p) / (numel(y)-numel(p)))
    fprintf('--\n')
    for i = 1:numel(p)
        fprintf('p_%d: (%.5g +- %.5g)\n', i, p(i), dp(i))
    end
end