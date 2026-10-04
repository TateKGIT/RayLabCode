% bobaFRET_run.m  -- build per-molecule histograms like SM_MakeHistogram_new, then run bobaFRET.
% Edit the settings below, then run the script.

folder    = 'C:\data\';        % folder holding the *_forvBFRET.dat files (one molecule per file)
binSize   = 0.02;              % same as the GUI's bin size box
startMu   = [0.9 0.45];        % same as the GUI's "Centers" box
startFWHM = [0.2 0.3];         % same as the GUI's "sigma" box (it is the FWHM)
nBoot     = 1000;
seed      = 1;                 % fixed seed = reproducible; change or remove for a fresh draw

% --- histogram settings copied from CalculateHistogram (keep them identical) ---
edges   = -0.2 + binSize/2 : binSize : 1.2 + binSize/2;
centers = edges(1) + binSize/2 : binSize : edges(end) - binSize/2;

files = dir(fullfile(folder, '*_forvBFRET.dat'));
if isempty(files), error('No *_forvBFRET.dat files found in %s', folder); end

perMolHist = zeros(numel(files), numel(centers));
for i = 1:numel(files)
    data = importdata(fullfile(folder, files(i).name));
    FRET = data(:,2) ./ (data(:,2) + data(:,1));
    perMolHist(i,:) = histcounts(FRET, edges) ./ size(data,1);   % normalized per molecule
end

% Sanity check: this should reproduce the GUI's saved _Histdata.dat (column 2)
guiHist = mean(perMolHist, 1);
fprintf('%d molecules; mean in-range fraction per molecule = %.3f (min %.3f)\n', ...
    numel(files), mean(sum(perMolHist,2)), min(sum(perMolHist,2)));

startGuess = reshape([startMu(:) startFWHM(:)]', 1, []);          % [mu1 w1 mu2 w2 ...]
results = bobaFRET(perMolHist, centers, startGuess, nBoot, seed);

save(fullfile(folder, 'boba_results.mat'), 'results', 'perMolHist', 'centers');

% Quick look at the bootstrap distributions of the populations
K = size(results.boot.pop, 2);
figure;
for k = 1:K
    subplot(1, K, k);
    histogram(results.boot.pop(:,k));
    xlabel(sprintf('Population of state %d', k));
    ylabel('Bootstrap replicates');
end

% Histogram with bootstrap error bars (+/- 1 standard error of each bin)
figure;
bar(centers, results.hist.mean, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [0.65 0.65 0.65]);
hold on;
errorbar(centers, results.hist.mean, results.hist.se, 'k', 'LineStyle', 'none', 'CapSize', 4);
hold off;
xlabel('FRET Efficiency'); ylabel('Normalized Count');
