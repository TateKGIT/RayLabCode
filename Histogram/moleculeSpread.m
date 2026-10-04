function results = moleculeSpread(folder, nBoot, seed, range, showPlot)
%MOLECULESPREAD  Spread of per-molecule FRET across molecules, with bootstrap uncertainty.
%
%   results = moleculeSpread()                 % asks for the folder
%   results = moleculeSpread(folder, nBoot, seed, range, showPlot)
%
%   folder  folder containing the *_forvBFRET.dat files (one molecule per file),
%           OR a cell array of full file paths; omit or [] to pick a folder in a
%           dialog.
%   nBoot   bootstrap replicates (default 1000)
%   seed    optional integer for rng(), for reproducible results
%   range   [lo hi] FRET window for the in-range mean (default [-0.2 1.2],
%           the GUI's histogram window)
%   showPlot  true (default) to draw histograms of the per-molecule values
%
%   For each molecule, FRET = col2 ./ (col2 + col1), as in SM_MakeHistogram_new,
%   and two summary values are computed:
%     median : median over ALL frames (robust to the A+D ~ 0 outliers)
%     mean   : mean over frames inside 'range' (biased low if many frames are
%              cut off at the upper edge, so read it with fracInRange)
%   The spread is std() across molecules of each value. The bootstrap resamples
%   MOLECULES with replacement to give a standard error and a 95% percentile
%   interval for that spread (and a standard error for the across-molecule mean).
%
%   The spread includes measurement noise in each molecule's value as well as real
%   molecule-to-molecule differences; this function does not separate the two.

if nargin < 1 || isempty(folder)
    folder = uigetdir(pwd, 'Select the folder containing the *_forvBFRET.dat files');
    if isequal(folder, 0), results = []; return; end
end
if nargin < 2 || isempty(nBoot), nBoot = 1000; end
if nargin >= 3 && ~isempty(seed), rng(seed); end
if nargin < 4 || isempty(range), range = [-0.2 1.2]; end
if nargin < 5 || isempty(showPlot), showPlot = true; end

if iscell(folder)
    paths = folder(:);                                   % list of full file paths
else
    files = dir(fullfile(folder, '*_forvBFRET.dat'));
    paths = fullfile(folder, {files.name})';
end
nMol = numel(paths);
if nMol < 2, error('Need at least 2 *_forvBFRET.dat files.'); end
names = cell(nMol, 1);

medE = nan(nMol, 1);  meanE = nan(nMol, 1);  fracIn = nan(nMol, 1);  nFrames = nan(nMol, 1);
for i = 1:nMol
    [~, nm, ext] = fileparts(paths{i});
    names{i} = [nm ext];
    data = importdata(paths{i});
    E    = data(:,2) ./ (data(:,2) + data(:,1));
    inR  = E >= range(1) & E <= range(2);          % NaN/Inf compare false, so they drop out
    medE(i)    = median(E, 'omitnan');
    meanE(i)   = mean(E(inR));                     % NaN if no frame is in range
    fracIn(i)  = mean(inR);
    nFrames(i) = size(data, 1);
end

results.names     = names;
results.nFrames   = nFrames;
results.medianE   = medE;
results.meanE     = meanE;
results.fracInRange = fracIn;
results.median    = spreadBoot(medE,  nBoot);
results.mean      = spreadBoot(meanE, nBoot);

fprintf('\nAcross-molecule spread of FRET: %d molecules, %d bootstrap replicates\n', nMol, nBoot);
fprintf('%-28s %8s %10s %10s %18s\n', 'Per-molecule value', 'mean', 'SD', 'SE of SD', '95% interval of SD');
printRow('median (all frames)', results.median);
printRow(sprintf('mean (FRET in [%g, %g])', range(1), range(2)), results.mean);
fprintf('In-range fraction per molecule: mean %.2f, min %.2f\n\n', mean(fracIn), min(fracIn));

if showPlot
    figure;
    subplot(1, 2, 1); histogram(medE);  xlabel('Per-molecule median FRET'); ylabel('Molecules');
    subplot(1, 2, 2); histogram(meanE); xlabel('Per-molecule in-range mean FRET'); ylabel('Molecules');
end
end


%% ======================== local functions ========================

function s = spreadBoot(x, nBoot)
% Bootstrap over molecules for the SD (and mean) of the values in x.
x = x(~isnan(x));
n = numel(x);
s.n    = n;
s.mean = mean(x);
s.sd   = std(x);                                   % SD across molecules (N-1)
sdB = zeros(nBoot, 1);  mB = zeros(nBoot, 1);
for b = 1:nBoot
    xb     = x(randi(n, n, 1));                    % resample molecules with replacement
    sdB(b) = std(xb);
    mB(b)  = mean(xb);
end
s.sdSE    = std(sdB);                              % bootstrap standard error of the SD
s.sdCI95  = pctl(sdB, [2.5 97.5]);
s.meanSE  = std(mB);                               % bootstrap standard error of the mean
s.sdBoot  = sdB;
end

function printRow(label, s)
fprintf('%-28s %8.3f %10.3f %10.3f %9.3f - %-8.3f\n', label, s.mean, s.sd, s.sdSE, s.sdCI95(1), s.sdCI95(2));
end

function q = pctl(x, p)
% Percentiles (0-100) of vector x, linear interpolation (no toolbox needed).
x   = sort(x(:));
n   = numel(x);
pos = 1 + (n - 1) * p / 100;
lo  = floor(pos);
hi  = ceil(pos);
q   = x(lo)' + (pos - lo) .* (x(hi)' - x(lo)');
end
