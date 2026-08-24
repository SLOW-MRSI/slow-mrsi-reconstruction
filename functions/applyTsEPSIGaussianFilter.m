function [xfDataGau, gaussianInfo] = applyTsEPSIGaussianFilter(xfData_ful, xfData_par, sigmaVec, runParallel, progress)
%APPLYTSEPSIGAUSSIANFILTER Spatial Gaussian smoothing for spectra.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

sigmaX = sigmaVec(1);
sigmaY = sigmaVec(2);
sigmaZ = sigmaVec(3);

filterSizeX = 2*ceil(2*sigmaX)+1;
filterSizeY = 2*ceil(2*sigmaY)+1;
filterSizeZ = 2*ceil(2*sigmaZ)+1;

x = -floor(filterSizeX/2):floor(filterSizeX/2);
y = -floor(filterSizeY/2):floor(filterSizeY/2);
z = -floor(filterSizeZ/2):floor(filterSizeZ/2);
[xx, yy, zz] = ndgrid(x, y, z);

gaussianKernel = exp(-((xx.^2)/(2*sigmaX^2) + (yy.^2)/(2*sigmaY^2) + (zz.^2)/(2*sigmaZ^2)));
gaussianKernel = gaussianKernel / max(gaussianKernel(:));

xfData_ful_Gau = zeros(size(xfData_ful), 'like', xfData_ful);
xfData_par_Gau = zeros(size(xfData_par), 'like', xfData_par);

if runParallel
    pool = gcp('nocreate');
    runParallel = ~isempty(pool);
end

doneCount = 0;
nTotal = size(xfData_ful,4) + size(xfData_par,4);

if runParallel
    updateTsEPSIProgress(progress, 0.75, 'Gaussian filtering', 0, 'Starting parallel filtering');
    xfData_ful_Gau = filterVolumeBatch(xfData_ful, xfData_ful_Gau, 'SLOW-ful spectra');
    xfData_par_Gau = filterVolumeBatch(xfData_par, xfData_par_Gau, 'SLOW-par spectra');
else
    doneCount = 0;
    for volIdx = 1:size(xfData_ful,4)
        xfData_ful_Gau(:,:,:,volIdx) = filterComplexVolume(xfData_ful(:,:,:,volIdx), sigmaVec);
        doneCount = doneCount + 1;
        updateTsEPSIProgress(progress, 0.75 + 0.12*doneCount/nTotal, 'Gaussian filtering', doneCount/nTotal, ...
            sprintf('Full spectra %d/%d', volIdx, size(xfData_ful,4)));
    end

    for volIdx = 1:size(xfData_par,4)
        xfData_par_Gau(:,:,:,volIdx) = filterComplexVolume(xfData_par(:,:,:,volIdx), sigmaVec);
        doneCount = doneCount + 1;
        updateTsEPSIProgress(progress, 0.75 + 0.12*doneCount/nTotal, 'Gaussian filtering', doneCount/nTotal, ...
            sprintf('Edited spectra %d/%d', volIdx, size(xfData_par,4)));
    end
end

updateTsEPSIProgress(progress, 0.87, 'Gaussian filtering', 1, 'Computing filtered difference');
xfData_dif_Gau = xfData_par_Gau - xfData_ful_Gau;
updateTsEPSIProgress(progress, 0.87, 'Gaussian filtering', 1, 'Gaussian filtering complete');

xfDataGau.ful = xfData_ful_Gau;
xfDataGau.par = xfData_par_Gau;
xfDataGau.dif = xfData_dif_Gau;

gaussianInfo.kernel = gaussianKernel;
gaussianInfo.effectiveVolume = sum(gaussianKernel(:));

    function outputData = filterVolumeBatch(inputData, outputData, labelText)
        nVol = size(inputData,4);
        futures = parallel.FevalFuture.empty;
        nextVol = 1;
        nQueued = min(nVol, pool.NumWorkers);

        for nq = 1:nQueued
            futures(end+1) = parfeval(pool, @filterComplexVolumeWithIndex, 2, nextVol, inputData(:,:,:,nextVol), sigmaVec); %#ok<AGROW>
            nextVol = nextVol + 1;
        end

        completedInBatch = 0;
        while completedInBatch < nVol
            [completedIdx, completedVolIdx, filteredVol] = fetchNext(futures);
            outputData(:,:,:,completedVolIdx) = filteredVol;

            doneCount = doneCount + 1;
            completedInBatch = completedInBatch + 1;
            updateTsEPSIProgress(progress, 0.75 + 0.12*doneCount/nTotal, 'Gaussian filtering', doneCount/nTotal, ...
                sprintf('%s %d/%d', labelText, completedInBatch, nVol));

            if nextVol <= nVol
                futures(completedIdx) = parfeval(pool, @filterComplexVolumeWithIndex, 2, nextVol, inputData(:,:,:,nextVol), sigmaVec);
                nextVol = nextVol + 1;
            else
                futures(completedIdx) = [];
            end
        end
    end

end

function [idx, out] = filterComplexVolumeWithIndex(idx, vol, sigmaVec)
out = filterComplexVolume(vol, sigmaVec);
end

function out = filterComplexVolume(vol, sigmaVec)
realPart = imgaussfilt3(real(vol), sigmaVec);
imagPart = imgaussfilt3(imag(vol), sigmaVec);
out = realPart + 1i*imagPart;
end
