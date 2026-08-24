% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

clc
clear

addpath('functions')
tic
%% User settings

% Set these paths for each computer/user.
cfg = struct();
dataPath = [pwd filesep];

rawDataPath = pwd;    % adjuest the path for your dataset
cfg.measID  = 19;       % adjuest the ID for your dataset

processedDataFolder = fullfile(pwd, 'processedData');
if ~exist(processedDataFolder, 'dir')
    mkdir(processedDataFolder)
end

cfg.maxAverages = 100;
cfg.useParallel = true;
cfg.showFigures = true;
cfg.showProgress = true;
cfg.saveFile = fullfile(processedDataFolder, 'xfData_raw.mat');

% Leave empty to use the sequence default: Nt/1.25.
cfg.metaboliteTimePoints = [];

% Reconstruction constants. Most users should not need to change these.
cfg.phaseSearchPpm = [2.6, 3.35];
cfg.phaseSearchSamples = 200;
cfg.phaseApodizationHz = 10;
cfg.spectralBandwidthHz = [];
cfg.scannerFrequencyMHz = [];
cfg.addPpm   = 0.0;          %  cfg.addPpm   = 0.15 for phantom
cfg.waterPpm = 4.65 + cfg.addPpm;
cfg.metabolitePpmOffset = 1.65-cfg.addPpm;
cfg.gaussianSigma = [1.5*0.8, 1.5*0.8, 1.5*0.8*0.08];

% Preview voxel grid for the diagnostic spectrum plots.
cfg.previewX = 8 + (-1:1)*2;
cfg.previewY = 8 + (-1:1)*2;
cfg.previewZ = 1;

%% Load one scan: metabolite and water data are separated from this scan

runParallel = setupTsEPSIParallel(cfg.useParallel);
progress = startTsEPSIProgress(cfg.showProgress, 'SLOW-MRSI reconstruction');
progressCleanup = onCleanup(@() closeTsEPSIProgress(progress));

twixFileNameCsi = ['fid_csi2d_', num2str(cfg.measID), '_twix.mat'];
opt = struct('reload', 1, 'measID', cfg.measID);

updateTsEPSIProgress(progress, 0.02, 'Loading raw data', 0, ...
    sprintf('Reading measurement %d', cfg.measID));
tic;
[rawData_met, CsiInfo, twix_obj_met] = loadEPSI_v1(twixFileNameCsi, dataPath, rawDataPath, opt);
toc;
updateTsEPSIProgress(progress, 0.08, 'Loading raw data', 1, 'Raw data loaded');

twix_head = twix_obj_met.hdr;
[scannerFrequencyMHz] = getScannerFrequencyMHzFromTwixHeader(twix_head);
if ~isempty(scannerFrequencyMHz)
    cfg.scannerFrequencyMHz = scannerFrequencyMHz;
    fprintf('Scanner frequency: %.6f MHz\n', cfg.scannerFrequencyMHz);
else
    cfg.scannerFrequencyMHz = 297.2;
    warning('Run_main:MissingScannerFrequency', ...
        'Could not find scanner frequency in twix header. Using fallback %.6f MHz.', cfg.scannerFrequencyMHz);
end
if cfg.scannerFrequencyMHz > 297
    if CsiInfo.Ny == 48
        cfg.spectralBandwidthHz = 1/500/2*1e6;
    elseif CsiInfo.Ny == 60
        cfg.spectralBandwidthHz = 1/540/2*1e6;
    elseif CsiInfo.Ny == 15
        cfg.spectralBandwidthHz = 1/260/2*1e6;
    end
else
    cfg.spectralBandwidthHz = 1/800/2*1e6;
end
fprintf('Spectral bandwidth: %.6f Hz\n', cfg.spectralBandwidthHz);

%% Reconstruction

[ktData_met, ktData_wat, reconInfo] = assembleTsEPSIData(rawData_met, CsiInfo, twix_obj_met, cfg, progress, runParallel);
clear rawData_met twix_obj_met

[ktData_shifted, ~] = alignTsEPSIEchoes(ktData_met, ktData_wat, reconInfo, runParallel, progress);
clear ktData_met ktData_wat

ktData_corrected = correctTsEPSIKx(ktData_shifted, reconInfo, progress);
clear ktData_shifted

xfData = reconstructTsEPSISpectra(ktData_corrected, runParallel, progress);
clear ktData_corrected

[xfDataGau, gaussianInfo] = applyTsEPSIGaussianFilter(xfData.ful, xfData.par, cfg.gaussianSigma, runParallel, progress);

[Pha_sel, Pha_s, X_ppm_met] = selectTsEPSIPhase(xfDataGau.par, cfg, runParallel, progress);

%% Diagnostic plots

if cfg.showFigures
    updateTsEPSIProgress(progress, 0.98, 'Diagnostic plots', 0, 'Generating preview figures');
    plotTsEPSIPreviewSpectra(xfDataGau, Pha_sel, Pha_s, X_ppm_met, cfg);
    plotTsEPSIMaps(xfData, reconInfo.Nz);
end

%% Save data

xfData_ful = xfData.ful;
xfData_par = xfData.par;
xfData_wat = xfData.wat;
updateTsEPSIProgress(progress, 0.99, 'Saving data', 0, cfg.saveFile);
save(cfg.saveFile, 'xfData_ful', 'xfData_par', 'xfData_wat', 'cfg', 'reconInfo', 'twix_head')
updateTsEPSIProgress(progress, 1.00, 'Complete', 1, 'Reconstruction finished');
clear progressCleanup

toc
