function xfData = reconstructTsEPSISpectra(ktData_corrected, runParallel, progress)
%RECONSTRUCTTSEPSISPECTRA Transform to image space, coil-combine, then FFT.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

updateTsEPSIProgress(progress, 0.60, 'Spatial FFT', 0/3, 'Water data');
ktData_wat = ktData_corrected.wat_e1 + ktData_corrected.wat_e2;
xtData_wat = mfft3d(ktData_wat);
clear ktData_wat

updateTsEPSIProgress(progress, 0.62, 'Spatial FFT', 1/3, 'SLOW-ful spectra');
ktData_ful = ktData_corrected.met_s1_e1 + ktData_corrected.met_s1_e2;
xtData_ful = mfft3d(ktData_ful);
clear ktData_ful

updateTsEPSIProgress(progress, 0.64, 'Spatial FFT', 2/3, 'SLOW-par spectra');
ktData_par = ktData_corrected.met_s2_e1 + ktData_corrected.met_s2_e2;
xtData_par = mfft3d(ktData_par);
clear ktData_par

xtData_comb_wat = zeros([size(xtData_wat,1),size(xtData_wat,2),size(xtData_wat,3),size(xtData_wat,4)], 'like', xtData_wat);
xtData_comb_ful = zeros(size(xtData_comb_wat), 'like', xtData_wat);
xtData_comb_par = zeros(size(xtData_comb_wat), 'like', xtData_wat);

Nz = size(xtData_wat,3);
if runParallel
    pool = gcp('nocreate');
    runParallel = ~isempty(pool);
end

if runParallel
    doneCount = 0;
    futures = parallel.FevalFuture.empty;
    nextSlice = 1;
    nQueued = min(Nz, pool.NumWorkers);

    updateTsEPSIProgress(progress, 0.66, 'Coil combination', 0, 'Starting parallel slice combination');
    for nq = 1:nQueued
        futures(end+1) = submitCoilCombine(nextSlice); %#ok<AGROW>
        nextSlice = nextSlice + 1;
    end

    while doneCount < Nz
        [completedIdx, n1, watSlice, fulSlice, parSlice] = fetchNext(futures);
        xtData_comb_wat(:,:,n1,:) = watSlice;
        xtData_comb_ful(:,:,n1,:) = fulSlice;
        xtData_comb_par(:,:,n1,:) = parSlice;

        doneCount = doneCount + 1;
        updateTsEPSIProgress(progress, 0.66 + 0.07*doneCount/Nz, 'Coil combination', doneCount/Nz, ...
            sprintf('Slices complete %d/%d', doneCount, Nz));

        if nextSlice <= Nz
            futures(completedIdx) = submitCoilCombine(nextSlice);
            nextSlice = nextSlice + 1;
        else
            futures(completedIdx) = [];
        end
    end
else
    for n1 = 1:Nz
        updateTsEPSIProgress(progress, 0.66 + 0.07*(n1-1)/Nz, 'Coil combination', (n1-1)/Nz, ...
            sprintf('Slice %d/%d', n1, Nz));
        [xtData_comb_wat(:,:,n1,:),xtData_comb_ful(:,:,n1,:),xtData_comb_par(:,:,n1,:)] = CoilCombineShort(xtData_wat(:,:,n1,:,:),xtData_ful(:,:,n1,:,:),xtData_par(:,:,n1,:,:),'SNR','linear');
    end
end

updateTsEPSIProgress(progress, 0.74, 'Spectral FFT', 0, 'Transforming combined spectra');
xfData.ful = flip(fftshift(fft(xtData_comb_ful,[],4),4),2); % flip RL direction to match the default MRI displayed direction
xfData.par = flip(fftshift(fft(xtData_comb_par,[],4),4),2);
xfData.wat = flip(fftshift(fft(xtData_comb_wat,[],4),4),2);
updateTsEPSIProgress(progress, 0.75, 'Spectral FFT', 1, 'Spectral reconstruction complete');

    function future = submitCoilCombine(sliceIndex)
        future = parfeval(pool, @coilCombineSliceWithIndex, 4, sliceIndex, ...
            xtData_wat(:,:,sliceIndex,:,:), xtData_ful(:,:,sliceIndex,:,:), xtData_par(:,:,sliceIndex,:,:));
    end
end

function [sliceIndex, watSlice, fulSlice, parSlice] = coilCombineSliceWithIndex(sliceIndex, watSliceInput, fulSliceInput, parSliceInput)
[watSlice, fulSlice, parSlice] = CoilCombineShort(watSliceInput, fulSliceInput, parSliceInput, 'SNR', 'linear');
end
