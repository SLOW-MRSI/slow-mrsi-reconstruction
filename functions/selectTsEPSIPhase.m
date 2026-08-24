function [Pha_sel, Pha_s, X_ppm_met] = selectTsEPSIPhase(xfData_par_Gau, cfg, runParallel, progress)
%SELECTTSEPSIPHASE Choose voxel-wise phase from the edited metabolite band.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

X_ppm_met = linspace(-cfg.spectralBandwidthHz/2,+cfg.spectralBandwidthHz/2,size(xfData_par_Gau,4)) ...
    / cfg.scannerFrequencyMHz + (cfg.waterPpm - cfg.metabolitePpmOffset);
Pha_s = linspaceGW(1,-1,cfg.phaseSearchSamples);

[n1Max,n2Max,n3Max,~] = size(xfData_par_Gau);
Pha_sel = zeros(n1Max,n2Max,n3Max,'like',1);
ppmMask = (X_ppm_met > cfg.phaseSearchPpm(1)) & (X_ppm_met < cfg.phaseSearchPpm(2));
phaseFactors = exp(1i*pi*Pha_s(:).');
phaseApodizationHz = cfg.phaseApodizationHz;

updateTsEPSIProgress(progress, 0.87, 'Phase selection', 0, 'Apodizing spectra');
xfData_par_Apo = apodizeSpectraTimeDim(xfData_par_Gau, phaseApodizationHz);
clear xfData_par_Gau

if runParallel
    pool = gcp('nocreate');
    runParallel = ~isempty(pool);
end

if runParallel
    doneCount = 0;
    futures = parallel.FevalFuture.empty;
    nextRow = 1;
    nQueued = min(n1Max, pool.NumWorkers);

    updateTsEPSIProgress(progress, 0.89, 'Phase selection', 0, 'Starting parallel row selection');
    for nq = 1:nQueued
        futures(end+1) = submitPhaseRow(nextRow); %#ok<AGROW>
        nextRow = nextRow + 1;
    end

    while doneCount < n1Max
        [completedIdx, n1, rowPhase] = fetchNext(futures);
        Pha_sel(n1,:,:) = rowPhase;

        doneCount = doneCount + 1;
        updateTsEPSIProgress(progress, 0.89 + 0.08*doneCount/n1Max, 'Phase selection', doneCount/n1Max, ...
            sprintf('Rows complete %d/%d', doneCount, n1Max));

        if nextRow <= n1Max
            futures(completedIdx) = submitPhaseRow(nextRow);
            nextRow = nextRow + 1;
        else
            futures(completedIdx) = [];
        end
    end
else
    for n1 = 1:n1Max
        updateTsEPSIProgress(progress, 0.89 + 0.08*(n1-1)/n1Max, 'Phase selection', (n1-1)/n1Max, ...
            sprintf('Row %d/%d', n1, n1Max));
        for n2 = 1:n2Max
            for n3 = 1:n3Max
                Pha_sel(n1,n2,n3) = selectOneVoxelPhase(xfData_par_Apo, ppmMask, phaseFactors, n1, n2, n3);
            end
        end
    end
end
updateTsEPSIProgress(progress, 0.97, 'Phase selection', 1, 'Voxel-wise phase selection complete');

    function future = submitPhaseRow(rowIndex)
        future = parfeval(pool, @selectOnePhaseRow, 2, rowIndex, ...
            xfData_par_Apo(rowIndex,:,:,:), ppmMask, phaseFactors, n2Max, n3Max);
    end
end

function spectraApo = apodizeSpectraTimeDim(spectra, apodizationHz)
nPts = size(spectra, 4);
t = linspace(0, 390/1000/1000*nPts, nPts);
weight = reshape(exp(-apodizationHz.*t), [1, 1, 1, nPts]);
if isa(spectra, 'single')
    weight = single(weight);
end

timeData = ifft(ifftshift(spectra, 4), [], 4) .* weight;
spectraApo = fftshift(fft(timeData, [], 4), 4);
end

function [rowIndex, rowPhase] = selectOnePhaseRow(rowIndex, rowData, ppmMask, phaseFactors, n2Max, n3Max)
rowPhase = zeros(1, n2Max, n3Max);
for n2 = 1:n2Max
    for n3 = 1:n3Max
        rowPhase(1,n2,n3) = selectOneVoxelPhaseFromRow(rowData, ppmMask, phaseFactors, n2, n3);
    end
end
end

function idx = selectOneVoxelPhase(xfData_par_Apo, ppmMask, phaseFactors, n1, n2, n3)
sig = squeeze(xfData_par_Apo(n1,n2,n3,:));
sigMask = sig(ppmMask);
temp_s = min(real(sigMask(:) * phaseFactors), [], 1);
[~,idx] = max(temp_s);
end

function idx = selectOneVoxelPhaseFromRow(rowData, ppmMask, phaseFactors, n2, n3)
sig = squeeze(rowData(1,n2,n3,:));
sigMask = sig(ppmMask);
temp_s = min(real(sigMask(:) * phaseFactors), [], 1);
[~,idx] = max(temp_s);
end
