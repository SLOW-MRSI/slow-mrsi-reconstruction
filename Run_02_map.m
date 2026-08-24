% SLOW-MRSI processing script
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

% -------------------------------------------------------------------------
% NOTE ON THE MAPS PRODUCED BY THIS SCRIPT
%
% The maps generated here (GABA, Glx, Cr, Cho, water) are DRAFT maps.
% Each voxel value comes from simple peak integration: the signal magnitude
% is summed over a fixed ppm window, with a small baseline offset removed.
%
% What these maps are for:
%   - selecting and checking the location of regions of interest,
%   - a fast overall quality check of the dataset (SNR, artefacts, lipid
%     contamination, frequency and phase behaviour) before further analysis.
%
% For quantification, use the filtered spectra saved by this script
% (xfData_Filt.mat) and process them with proper spectral fitting.
% -------------------------------------------------------------------------

clc
clear

addpath('functions')

%% User settings

processedDataFolder = fullfile(pwd, 'processedData');
if ~exist(processedDataFolder, 'dir')
    mkdir(processedDataFolder)
end

inputFile = fullfile(processedDataFolder, 'xfData_raw.mat');

% Manual overrides for parameters loaded from Run_main.
cfgManual = struct();
cfgManual.metabolitePpmOffset = 1.65;
zeroOrderPha = 0.0;

% Spatial smoothing after frequency/phase correction.
medianKernel = [1 1 1];
mapDisplayPercentile = 50;

% Display and interactive spectrum viewer settings.
showBrainMaskPreview = 'on';  % Options: 'on' or 'off'
adjustBrainMask = 'on';       % Options: 'on' or 'off'
brainMaskThresholdFraction = 0.06;
brainMaskSkullMarginVoxels = 3;

mapFigurePosition = [1, 1];
mapFigureSize = [11, 8];
spectrumFigureSize = [12, 5];

% Lipid/water signal ratio from the strongest per-voxel spectral peaks.
lipidWaterTopFraction = 0.01;

outputFilteredFile = fullfile(processedDataFolder, 'xfData_Filt.mat');

% Spatial smoothing parameters for the single loaded dataset.
gaussianSigma    = [1.4, 1.4, 1]*0.50;
gaussianSigmaDif = [1, 1, 1]*0.65*1.0;
phaseFilterSigma = [1, 1, 1]*0.75;

fprintf('Gaussian kernel settings:\n \n');
reportGaussianKernel('main spatial filter', gaussianSigma);
reportGaussianKernel('difference spatial filter', gaussianSigmaDif);
reportGaussianKernel('phase-reference filter', phaseFilterSigma);
fprintf('moving median filter: kernel = [%d %d %d], voxels = %d\n \n', ...
    medianKernel(1), medianKernel(2), medianKernel(3), prod(medianKernel));


%% Load data

data = load(inputFile);

if isfield(data, 'xfData')
    xfData = data.xfData;
    xfData_ful = xfData.ful;
    xfData_par = xfData.par;
    xfData_wat = xfData.wat;
else
    xfData_ful = data.xfData_ful;
    xfData_par = data.xfData_par;
    xfData_wat = data.xfData_wat;
end

xfData_dif = xfData_par - xfData_ful;

%% Lipid-to-water ratio from top per-voxel spectral peaks

lipidVoxelPeak = max(abs(xfData_ful), [], 4);
waterVoxelPeak = max(abs(xfData_wat), [], 4);
lipidPeakValues = lipidVoxelPeak(:);
waterPeakValues = waterVoxelPeak(:);
lipidPeakValues = lipidPeakValues(isfinite(lipidPeakValues));
waterPeakValues = waterPeakValues(isfinite(waterPeakValues));

lipidTopCount = max(1, ceil(numel(lipidPeakValues) * lipidWaterTopFraction));
waterTopCount = max(1, ceil(numel(waterPeakValues) * lipidWaterTopFraction));

lipidPeakValues = sort(lipidPeakValues, 'descend');
waterPeakValues = sort(waterPeakValues, 'descend');

lipidTopMedian = median(lipidPeakValues(1:lipidTopCount));
waterTopMedian = median(waterPeakValues(1:waterTopCount));
lipidToWaterRatio = lipidTopMedian / waterTopMedian;

lipidWaterRatio = struct();
lipidWaterRatio.topFraction = lipidWaterTopFraction;
lipidWaterRatio.method = 'top fraction of per-voxel spectral peak magnitude';
lipidWaterRatio.lipidTopMedian = lipidTopMedian;
lipidWaterRatio.waterTopMedian = waterTopMedian;
lipidWaterRatio.lipidToWaterRatio = lipidToWaterRatio;
lipidWaterRatio.lipidTopCount = lipidTopCount;
lipidWaterRatio.waterTopCount = waterTopCount;
lipidWaterRatio.lipidVoxelPeak = lipidVoxelPeak;
lipidWaterRatio.waterVoxelPeak = waterVoxelPeak;
save(fullfile(processedDataFolder, 'lipidWaterRatio.mat'), 'lipidWaterRatio')

if isfield(data, 'cfg')
    cfg = data.cfg;
else
    warning('Run_map:MissingCfg', ...
        'No cfg saved in %s. Using legacy Run_map defaults for reconstruction parameters.', inputFile);
    cfg = struct();
end

manualFields = fieldnames(cfgManual);
for idxField = 1:numel(manualFields)
    cfg.(manualFields{idxField}) = cfgManual.(manualFields{idxField});
end
scannerFrequencyManuallySet = isfield(cfgManual, 'scannerFrequencyMHz') && ~isempty(cfgManual.scannerFrequencyMHz);
spectralBandwidthManuallySet = isfield(cfgManual, 'spectralBandwidthHz') && ~isempty(cfgManual.spectralBandwidthHz);

if ~scannerFrequencyManuallySet && isfield(data, 'twix_head')
    scannerFrequencyMHz = getScannerFrequencyMHzFromTwixHeader(data.twix_head);
    if ~isempty(scannerFrequencyMHz)
        cfg.scannerFrequencyMHz = scannerFrequencyMHz;
        fprintf('Scanner frequency: %.6f MHz \n', cfg.scannerFrequencyMHz);
    end
end

if ~isfield(cfg, 'scannerFrequencyMHz') || isempty(cfg.scannerFrequencyMHz)
    cfg.scannerFrequencyMHz = 123.2;
end
if ~isfield(cfg, 'spectralBandwidthHz')
    if cfg.scannerFrequencyMHz > 297
        cfg.spectralBandwidthHz = 1/500/2*1e6;
    else
        cfg.spectralBandwidthHz = 1/800/2*1e6;
    end
end
fprintf('Spectral bandwidth: %.6f Hz\n', cfg.spectralBandwidthHz);
if ~isfield(cfg, 'waterPpm') || isempty(cfg.waterPpm)
    cfg.waterPpm = 4.65;
end
if ~isfield(cfg, 'metabolitePpmOffset') || isempty(cfg.metabolitePpmOffset)
    cfg.metabolitePpmOffset = cfgManual.metabolitePpmOffset;
end
if ~isfield(cfg, 'phaseSearchSamples') || isempty(cfg.phaseSearchSamples)
    cfg.phaseSearchSamples = 200;
end
if ~isfield(cfg, 'phaseApodizationHz') || isempty(cfg.phaseApodizationHz)
    cfg.phaseApodizationHz = 10;
end
if ~isfield(cfg, 'phaseSearchPpm') || isempty(cfg.phaseSearchPpm)
    cfg.phaseSearchPpm = [2.6, 3.5];
end
if cfg.scannerFrequencyMHz > 297
    outputLinewidthHz = struct('ful', 15, 'par', 15, 'dif', 25, 'wat', 10);
else
    outputLinewidthHz = struct('ful', 10, 'par', 10, 'dif', 15, 'wat', 10);
end

EPSIDur = 1e6 / cfg.spectralBandwidthHz;
BW = cfg.spectralBandwidthHz;
B0 = cfg.scannerFrequencyMHz;

X_ppm_wat = linspace(-BW/2, +BW/2, size(xfData_ful,4))/B0 + cfg.waterPpm;
X_ppm = X_ppm_wat - cfg.metabolitePpmOffset;
%% Whole-volume spectral-quality metric

% Apodize the raw difference spectrum in the time domain with the same
% linewidth setting used for the filtered difference-spectrum output.
spectralQualityQC = struct();
spectralQualityQC.excludedHighestFraction = 0.01;

ppmRMS = (X_ppm > 3.2) & (X_ppm < 3.6);
xfData_dif_Apo = applySpectralApodization4D(xfData_dif, EPSIDur, outputLinewidthHz.dif);
spectralQualityRms = sqrt(mean(abs(xfData_dif_Apo(:,:,:,ppmRMS)).^2, 4));
spectralQualityWaterAbsSum = abs(sum(xfData_wat, 4));
spectralQualityMatrix = spectralQualityRms .* spectralQualityWaterAbsSum;

spectralQualityValues = sort(spectralQualityMatrix(isfinite(spectralQualityMatrix)));
if isempty(spectralQualityValues)
    spectralQualityKeepCount = 0;
    spectralQualityMean = NaN;
else
    spectralQualityKeepCount = max(1, floor((1-spectralQualityQC.excludedHighestFraction) * numel(spectralQualityValues)));
    spectralQualityMean = mean(spectralQualityValues(end/2:spectralQualityKeepCount)) * 1e4;
end


spectralQualityQC.method = 'rms(apodized(xfData_par - xfData_ful), 4) .* abs(sum(xfData_wat, 4))';
spectralQualityQC.apodizationLinewidthHz = outputLinewidthHz.dif;
spectralQualityQC.voxelwiseMatrix = spectralQualityMatrix;
spectralQualityQC.meanAllVoxels = spectralQualityMean;
save(fullfile(processedDataFolder, 'spectralQualityQC.mat'), 'spectralQualityQC')

fprintf('\nLipid/water ratio from top %.1f%% per-voxel spectral peaks: lipid median %.6g, water median %.6g, ratio %.6g\n', ...
    100*lipidWaterTopFraction, lipidTopMedian, waterTopMedian, lipidToWaterRatio);
fprintf('Spectral-quality QC mean: %.6g\n', spectralQualityMean);

phaseSearchSamples = cfg.phaseSearchSamples;
phaseApodizationHz = cfg.phaseApodizationHz;
phaseSearchPpm = cfg.phaseSearchPpm;


%% Phase-reference smoothing

xfData_par_phaseFiltered = applyComplexSpatialFilter(xfData_par, 'gaussian', phaseFilterSigma, medianKernel);

%% Brain mask

ppmMaskWater = (X_ppm_wat > 2.6) & (X_ppm_wat < 6.3);
Map_wat = sum(xfData_wat(:,:,:,ppmMaskWater), 4);
Map_wat_abs = abs(Map_wat);
Map_wat_s = imgaussfilt3(Map_wat_abs, 0.3);

if strcmpi(adjustBrainMask, 'on')
    [brainMask, brainMaskInfo, mask_main, mask_wat] = adjustBrainMaskInteractively( ...
        Map_wat_s, brainMaskThresholdFraction, brainMaskSkullMarginVoxels, showBrainMaskPreview);
end
if ~exist('brainMask', 'var')
    [brainMask, mask_main, mask_wat, thresholdValue, brainMaskSkullMarginVoxels] = buildBrainMaskFromWaterMap( ...
        Map_wat_s, brainMaskThresholdFraction, brainMaskSkullMarginVoxels);
    brainMaskInfo = struct('thresholdFraction', brainMaskThresholdFraction, ...
        'skullMarginVoxels', brainMaskSkullMarginVoxels, 'thresholdValue', thresholdValue);
    plotBrainMaskPreview(brainMask, showBrainMaskPreview, Map_wat_s, brainMaskInfo);
end
fprintf('\nBrain mask: threshold fraction %.4g, skull margin %d voxels\n', ...
    brainMaskInfo.thresholdFraction, brainMaskInfo.skullMarginVoxels);

%% Frequency and phase correction

Pha_s = linspaceGW(1, -1, phaseSearchSamples);
phaseFactors = exp(1i*pi*Pha_s(:).');
mid = round(0.5 + size(xfData_ful,4)/2);

[n1Max, n2Max, n3Max, ~] = size(xfData_ful);
Pha_sel = zeros(n1Max, n2Max, n3Max, 'like', 1);
ppmMeta = (X_ppm > phaseSearchPpm(1)) & (X_ppm < phaseSearchPpm(2));
ppmLip = (X_ppm < 1.85);
dPPM   = X_ppm(2) - X_ppm(1);

xfData_ful_Corr = zeros(size(xfData_ful), 'like', xfData_ful);
xfData_par_Corr = zeros(size(xfData_par), 'like', xfData_par);
xfData_wat_Corr = zeros(size(xfData_wat), 'like', xfData_wat);
xfData_wat_pkCentered = zeros(size(xfData_wat), 'like', xfData_wat);

clear temp temp_ful temp_par temp_wat
parfor n1 = 1:n1Max
    for n2 = 1:n2Max
        for n3 = 1:n3Max
            sigWat = squeeze(xfData_wat(n1,n2,n3,:));
            sigPar = squeeze(xfData_par_phaseFiltered(n1,n2,n3,:));

            [~, idxWater] = max(abs(sigWat));
            shiftAmt = mid - double(idxWater);
            sigWatShifted = circshift(sigWat, shiftAmt);
            sigParShifted = circshift(sigPar, shiftAmt);
            sigParApoShifted = squeeze(apodizeGW_v2(sigParShifted, EPSIDur, phaseApodizationHz, 0));
            xfData_wat_pkCentered(n1,n2,n3,:) = sigWatShifted;

            sigMeta = sigParApoShifted(ppmMeta);
            temp_s = min(real(sigMeta(:) * phaseFactors), [], 1);
            [~, idxPhase] = max(temp_s);
            Pha_sel(n1,n2,n3) = idxPhase;

            if  abs(shiftAmt)*dPPM  > 0.4
                shiftAmt = 0;
            end

            temp_ful = circshift(squeeze(xfData_ful(n1,n2,n3,:)), shiftAmt);
            temp_par = circshift(squeeze(xfData_par(n1,n2,n3,:)), shiftAmt);
            temp_wat = sigWatShifted;

            if brainMask(n1,n2,n3)
                phaseFactor = phaseFactors(idxPhase);
            else
                phaseFactor = 0;
            end

            if brainMask(n1,n2,n3) == 0
                xfData_ful_Corr(n1,n2,n3,:) = temp_ful .* phaseFactor*0;
                xfData_par_Corr(n1,n2,n3,:) = temp_par .* phaseFactor*0;
            else
                xfData_ful_Corr(n1,n2,n3,:) = temp_ful .* phaseFactor .* exp(2i*pi*zeroOrderPha);
                xfData_par_Corr(n1,n2,n3,:) = temp_par .* phaseFactor .* exp(2i*pi*zeroOrderPha);

            end
            xfData_wat_Corr(n1,n2,n3,:) = temp_wat .* phaseFactor .* exp(2i*pi*zeroOrderPha);
        end
    end
end

%% Spatial smoothing

fprintf('4D spatial filter: full/par/water = gaussian sigma [%g %g %g]; difference = median [%d %d %d] then gaussian sigma [%g %g %g]\n', ...
    gaussianSigma(1), gaussianSigma(2), gaussianSigma(3), ...
    medianKernel(1), medianKernel(2), medianKernel(3), ...
    gaussianSigmaDif(1), gaussianSigmaDif(2), gaussianSigmaDif(3));

xfData_ful_Filt = applyComplexSpatialFilter(xfData_ful_Corr, 'gaussian', gaussianSigma, medianKernel);
xfData_par_Filt = applyComplexSpatialFilter(xfData_par_Corr, 'gaussian', gaussianSigma, medianKernel);
xfData_wat_Filt = applyComplexSpatialFilter(xfData_wat_Corr, 'gaussian', gaussianSigma, medianKernel);

xfData_dif_Raw = xfData_par_Corr - xfData_ful_Corr;
xfData_dif_Med = applyComplexSpatialFilter(xfData_dif_Raw, 'median', gaussianSigmaDif, medianKernel);
xfData_dif_Filt = applyComplexSpatialFilter(xfData_dif_Med, 'gaussian', gaussianSigmaDif, medianKernel);
clear xfData_dif_Raw xfData_dif_Med

fprintf('Preparing apodized filtered output: full/par/water = %g Hz, difference = %g Hz\n', ...
    outputLinewidthHz.ful, outputLinewidthHz.dif);
xfData_Filt = struct();
xfData_Filt.ful = applySpectralApodization4D(xfData_ful_Filt, EPSIDur, outputLinewidthHz.ful);
xfData_Filt.par = applySpectralApodization4D(xfData_par_Filt, EPSIDur, outputLinewidthHz.par);
xfData_Filt.dif = applySpectralApodization4D(xfData_dif_Filt, EPSIDur, outputLinewidthHz.dif);
xfData_Filt.wat = applySpectralApodization4D(xfData_wat_Filt, EPSIDur, outputLinewidthHz.wat);
xfData_Filt.linewidthHz = outputLinewidthHz;
xfData_Filt.EPSIDur = EPSIDur;
xfData_Filt.BW = BW;
xfData_Filt.B0 = B0;
xfData_Filt.X_ppm = X_ppm;
xfData_Filt.X_ppm_wat = X_ppm_wat;
save(outputFilteredFile, 'xfData_Filt', '-v7.3')
fprintf('Saved apodized filtered spectra to %s\n', outputFilteredFile);


%% Generate maps

ppmGABA = (X_ppm > 2.85) & (X_ppm < 3.16);
ppmGlx = (X_ppm > 3.67) & (X_ppm < 3.88);
ppmCho = (X_ppm > 3.11) & (X_ppm < 3.35);
ppmCr = (X_ppm > 2.85) & (X_ppm < 3.11);

if islogical(ppmGABA), ppmGABA = ppmGABA(:).'; end
if islogical(ppmGlx), ppmGlx = ppmGlx(:).'; end
if islogical(ppmCho), ppmCho = ppmCho(:).'; end
if islogical(ppmCr), ppmCr = ppmCr(:).'; end
if islogical(ppmLip), ppmLip = ppmLip(:).'; end

Map_GABA = zeros(n1Max, n2Max, n3Max, 'like', real(xfData_dif_Filt(1)));
Map_Glx = zeros(n1Max, n2Max, n3Max, 'like', real(xfData_dif_Filt(1)));
Map_Cr = zeros(n1Max, n2Max, n3Max, 'like', real(xfData_dif_Filt(1)));
Map_Cho = zeros(n1Max, n2Max, n3Max, 'like', real(xfData_dif_Filt(1)));
Map_Lip = zeros(n1Max, n2Max, n3Max, 'like', real(xfData_dif_Filt(1)));

parfor n1 = 1:n1Max
    for n2 = 1:n2Max
        for n3 = 1:n3Max
            sigDif = squeeze(xfData_dif_Filt(n1,n2,n3,:));
            sigDifApo = apodizeGW_v2(sigDif,  EPSIDur, 20,0);

            sigSum = squeeze(xfData_ful_Filt(n1,n2,n3,:) + xfData_par_Filt(n1,n2,n3,:));
            sigSumApo = 0.5*apodizeGW_v2(sigSum, EPSIDur,10, 0);

            Map_GABA(n1,n2,n3) = sum(abs(sigDifApo(ppmGABA))) - 0.5*sum(ppmGABA*min(abs(sigDifApo(ppmGABA))));
            Map_Glx(n1,n2,n3) = sum(abs(sigDifApo(ppmGlx)))   - 0.5*sum(ppmGlx*min(abs(sigDifApo(ppmGlx))));
            Map_Cr(n1,n2,n3) = sum(abs(sigSumApo(ppmCr)))     - 0.5*sum(ppmCr*min(abs(sigSumApo(ppmCr))));
            Map_Cho(n1,n2,n3) = sum(abs(sigSumApo(ppmCho)))   - 0.5*sum(ppmCho*min(abs(sigSumApo(ppmCho))));
            Map_Lip(n1,n2,n3) = sum(abs(sigSumApo(ppmLip)));
        end
    end
end

%% Display maps

Mask = Map_Lip.*0 + 1;
Mask(Map_Lip > 0.20) = 0; % not used here

Map_GABA_filt = Map_GABA.*Mask;
Map_Glx_filt = Map_Glx.*Mask;
Map_Cr_filt = Map_Cr.*Mask;
Map_Cho_filt = Map_Cho.*Mask;

Nz = size(Map_GABA,3);
mapFigureNumber = 5;
figMap = figure(mapFigureNumber); clf(figMap)
connectMapSpectrumViewer('setupFigure', figMap, spectrumFigureSize);
t = tiledlayout(5, Nz, 'Padding', 'none', 'TileSpacing', 'tight');

Map_wat_display = abs(sum(xfData_wat(:,:,:,:), 4));
Sc1 = [
    percentileUpperLimit(Map_GABA_filt, mapDisplayPercentile), ...
    percentileUpperLimit(Map_Glx_filt, mapDisplayPercentile), ...
    percentileUpperLimit(Map_Cr_filt, mapDisplayPercentile), ...
    percentileUpperLimit(Map_Cho_filt, mapDisplayPercentile), ...
    percentileUpperLimit(Map_wat_display(Map_GABA_filt>0), mapDisplayPercentile)];
fprintf('Display scale upper limits from %.1fth percentile: GABA %.4g, Glx %.4g, Cr %.4g, Cho %.4g, water %.4g\n', ...
    mapDisplayPercentile, Sc1(1), Sc1(2), Sc1(3), Sc1(4), Sc1(5));
fprintf('Click any map voxel in figure %d to show spectra. Use keyboard arrow keys to move the selected voxel.\n', ...
    mapFigureNumber);

Sc2 = [2.8, 2.8 ,2,2,2]*1.3;
for nz = 1:Nz
    ax = nexttile(t, nz);
    img = imagesc(ax, Map_GABA_filt(:,:,nz), [0, Sc1(1)*Sc2(1)]);
    connectMapSpectrumViewer('enableClick', img, 'GABA', nz);
    title(ax, 'GABA');
    axis(ax, 'off'); axis(ax, 'image'); axis(ax, 'tight');

    ax = nexttile(t, Nz + nz);
    img = imagesc(ax, Map_Glx_filt(:,:,nz), [0, Sc1(2)*Sc2(2)]);
    connectMapSpectrumViewer('enableClick', img, 'Glx', nz);
    title(ax, 'Glx');
    axis(ax, 'off'); axis(ax, 'image'); axis(ax, 'tight');

    ax = nexttile(t, 2*Nz + nz);
    img = imagesc(ax, Map_Cr_filt(:,:,nz), [0, Sc1(3)*Sc2(3)]);
    connectMapSpectrumViewer('enableClick', img, 'Cr', nz);
    title(ax, 'Cr');
    axis(ax, 'off'); axis(ax, 'image'); axis(ax, 'tight');

    ax = nexttile(t, 3*Nz + nz);
    img = imagesc(ax, Map_Cho_filt(:,:,nz), [0, Sc1(4)*Sc2(4)]);
    connectMapSpectrumViewer('enableClick', img, 'Cho', nz);
    title(ax, 'Cho');
    axis(ax, 'off'); axis(ax, 'image'); axis(ax, 'tight');

    ax = nexttile(t, 4*Nz + nz);
    img = imagesc(ax, Map_wat_display(:,:,nz), [0, Sc1(5)*Sc2(5)]);
    connectMapSpectrumViewer('enableClick', img, 'water', nz);
    title(ax, 'water');
    axis(ax, 'off'); axis(ax, 'image'); axis(ax, 'tight');

    colormap('turbo')
end

set(figMap, 'color', 'w')
set(figMap, 'Units', 'inches', 'Position', [mapFigurePosition, mapFigureSize]);
print(fullfile(processedDataFolder, 'Figure_draftMap'), '-dpng', '-r300');
