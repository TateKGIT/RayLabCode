function results = bobaFRET(perMolHist, centers, startGuess, nBoot, seed)
%BOBAFRET  Bootstrap-based analysis (BOBA) of a molecule-weighted FRET histogram.
%
%   results = bobaFRET(perMolHist, centers, startGuess)
%   results = bobaFRET(perMolHist, centers, startGuess, nBoot, seed)
%
%   Inputs
%     perMolHist  nMolecules x nBins. Row i is molecule i's histogram, normalized
%                 exactly as in SM_MakeHistogram_new (counts ./ number of frames).
%     centers     1 x nBins bin centers (histogramData(1,:) in the GUI).
%     startGuess  [mu1 fwhm1 mu2 fwhm2 ...]  starting values for the Gaussian fit.
%                 Same convention as the GUI: the "sigma" box is really the FWHM.
%     nBoot       number of bootstrap replicates (default 1000).
%     seed        optional integer for rng(), for reproducible results.
%
%   Procedure
%     1. Fit the sum of K Gaussians to the mean histogram of all molecules.
%     2. Repeat nBoot times: draw nMolecules molecules WITH replacement, average
%        their histograms (molecular weighting), refit starting from step 1's
%        solution, and record each state's center, FWHM and population.
%     3. Report the mean, standard deviation and 95% percentile interval over
%        the replicates (sigma_BOBA).
%
%   Population of state k = area of Gaussian k / total area of all Gaussians.
%   Amplitudes are solved by non-negative linear least squares; only centers and
%   widths are searched by fminsearch (as in the GUI, but with lsqnonneg instead
%   of a penalty for negative heights). States are sorted by center in every fit,
%   so state k means the same state in every replicate.
%
%   The per-bin bootstrap standard error of the histogram (results.hist.se) is what
%   the GUI draws as error bars on the histogram bars.
%
%   Output: struct with fields  full, boot, summary, hist, nMolecules  (see code).

if nargin < 4 || isempty(nBoot), nBoot = 1000; end
if nargin >= 5 && ~isempty(seed), rng(seed); end

centers    = centers(:)';
startGuess = startGuess(:)';
if mod(numel(startGuess), 2) ~= 0
    error('startGuess must be [mu1 fwhm1 mu2 fwhm2 ...] (an even number of values).');
end
if size(perMolHist, 2) ~= numel(centers)
    error('perMolHist has %d columns but centers has %d elements.', size(perMolHist, 2), numel(centers));
end
nMol = size(perMolHist, 1);
K    = numel(startGuess) / 2;
if nMol < 2
    error('Need at least 2 molecules to bootstrap.');
end

%% 1. Fit to the full data set
yAll = mean(perMolHist, 1);
[full, okFull] = fitGaussians(centers, yAll, startGuess);
if ~okFull
    warning('bobaFRET:fullFit', 'Fit to the full data set did not converge; check startGuess.');
end
SS_res  = sum((yAll - full.yhat).^2);
SS_tot  = sum((yAll - mean(yAll)).^2);
full.R2 = 1 - SS_res / SS_tot;

%% 2. Bootstrap over molecules
pop  = nan(nBoot, K);
mu   = nan(nBoot, K);
fwhm = nan(nBoot, K);
yBoot = zeros(nBoot, numel(centers));         % bootstrap histograms, for per-bin error bars
for b = 1:nBoot
    idx = randi(nMol, nMol, 1);               % resample molecules with replacement
    y   = mean(perMolHist(idx, :), 1);        % bootstrap cumulated histogram
    yBoot(b, :) = y;
    [fit, ok] = fitGaussians(centers, y, full.lambda);
    if ok
        pop(b, :)  = fit.pop;
        mu(b, :)   = fit.mu;
        fwhm(b, :) = fit.fwhm;
    end
end
valid   = all(~isnan(pop), 2);
nValid  = sum(valid);
if nValid < 0.9 * nBoot
    warning('bobaFRET:failures', '%d of %d replicates did not converge; results may be unreliable.', nBoot - nValid, nBoot);
end
if nValid < 2
    error('Too few converged replicates to summarize.');
end

%% 3. Summarize
summary.popMean  = mean(pop(valid, :), 1);
summary.popStd   = std(pop(valid, :), 0, 1);          % sigma_BOBA of the population
summary.popCI95  = [pctl(pop(valid, :), 2.5); pctl(pop(valid, :), 97.5)];
summary.muMean   = mean(mu(valid, :), 1);
summary.muStd    = std(mu(valid, :), 0, 1);
summary.fwhmMean = mean(fwhm(valid, :), 1);
summary.fwhmStd  = std(fwhm(valid, :), 0, 1);

results.full        = full;
results.boot.pop    = pop;
results.boot.mu     = mu;
results.boot.fwhm   = fwhm;
results.boot.nBoot  = nBoot;
results.boot.nValid = nValid;
results.summary     = summary;
results.nMolecules  = nMol;

% Per-bin uncertainty of the histogram itself (all replicates; independent of the fits)
results.hist.centers = centers;
results.hist.mean    = yAll;
results.hist.se      = std(yBoot, 0, 1);                     % bootstrap standard error of each bin
results.hist.ci95    = [pctl(yBoot, 2.5); pctl(yBoot, 97.5)]; % 95% percentile interval of each bin

%% Print a short report
fprintf('\nBOBA FRET: %d molecules, %d/%d replicates converged. Full-data fit R^2 = %.4f\n', nMol, nValid, nBoot, full.R2);
fprintf('%-6s %-18s %-18s %-24s %-s\n', 'State', 'Center (mean+/-sd)', 'FWHM (mean+/-sd)', 'Population (mean+/-sd)', '95% interval');
for k = 1:K
    fprintf('%-6d %6.3f +/- %-9.3f %6.3f +/- %-9.3f %6.3f +/- %-15.3f [%.3f, %.3f]\n', k, ...
        summary.muMean(k), summary.muStd(k), summary.fwhmMean(k), summary.fwhmStd(k), ...
        summary.popMean(k), summary.popStd(k), summary.popCI95(1, k), summary.popCI95(2, k));
end
fprintf('\n');
end


%% ======================== local functions ========================

function [fit, ok] = fitGaussians(t, y, lambda0)
% Fit a sum of Gaussians to histogram y (centers t). lambda = [mu1 w1 mu2 w2 ...],
% w = FWHM. Returns states sorted by center.
opts = optimset('TolX', 1e-4, 'TolFun', 1e-4, 'MaxIter', 5000, ...
                'MaxFunEvals', 5000 * numel(lambda0), 'Display', 'off');
[lambda, ~, flag] = fminsearch(@(l) fitError(l, t, y), lambda0, opts);
[err, c] = fitError(lambda, t, y);
ok = (flag > 0) && (err < 1e2);               % err = 1e3 marks an out-of-bounds solution

mu = lambda(1:2:end);
w  = abs(lambda(2:2:end));
[mu, order] = sort(mu);
w  = w(order);
c  = c(order)';
area = c .* 0.60056120439323 .* w .* sqrt(pi);  % area of c*exp(-((x-mu)/(0.60056*w))^2)

fit.lambda = reshape([mu; w], 1, []);          % sorted, ready to use as next start
fit.mu     = mu(:)';
fit.fwhm   = w(:)';
fit.sigma  = w(:)' / (2 * sqrt(2 * log(2)));   % conventional standard deviation
fit.amp    = c;
fit.pop    = area / sum(area);
fit.yhat   = zeros(size(t));
for k = 1:numel(mu)
    fit.yhat = fit.yhat + c(k) * gaussianFWHM(t, mu(k), w(k));
end
if sum(area) <= 0 || any(isnan(fit.pop))
    ok = false;
end
end

function [err, c] = fitError(lambda, t, y)
% Residual norm with amplitudes solved by non-negative least squares.
K  = numel(lambda) / 2;
mu = lambda(1:2:end);
w  = abs(lambda(2:2:end));
dt = t(2) - t(1);
% Keep centers inside the histogram and widths between one bin and the full span.
if any(mu < t(1)) || any(mu > t(end)) || any(w < dt) || any(w > (t(end) - t(1)))
    err = 1e3;
    c   = zeros(K, 1);
    return
end
A = zeros(numel(t), K);
for k = 1:K
    A(:, k) = gaussianFWHM(t, mu(k), w(k))';
end
c   = lsqnonneg(A, y(:));
err = norm(A * c - y(:));
end

function g = gaussianFWHM(x, peakPosition, width)
% Same model as gaussian() in SM_MakeHistogram_new: width is the FWHM.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
end

function q = pctl(X, p)
% Percentile p (0-100) of each column, linear interpolation (no toolbox needed).
X   = sort(X, 1);
n   = size(X, 1);
pos = 1 + (n - 1) * p / 100;
lo  = floor(pos);
hi  = ceil(pos);
q   = X(lo, :) + (pos - lo) .* (X(hi, :) - X(lo, :));
end
