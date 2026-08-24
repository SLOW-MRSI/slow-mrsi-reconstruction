function [ktData_shifted, echoInfo] = alignTsEPSIEchoes(ktData_met, ktData_wat, reconInfo, runParallel, progress)
%ALIGNTSEPSIECHOES Align odd/even EPSI echoes using the embedded water signal.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

Nx = reconInfo.Nx;
Nz = reconInfo.Nz;
Nc = reconInfo.Nc;

Ny = size(ktData_wat.e1,2);
Ny_s = round((Ny+1)/2) + (-1:1);
Nz_s = round((Nz+1)/2);

wat_e1 = ktData_wat.e1;
wat_e2 = ktData_wat.e2;
met_s1_e1 = ktData_met.s1_e1;
met_s1_e2 = ktData_met.s1_e2;
met_s2_e1 = ktData_met.s2_e1;
met_s2_e2 = ktData_met.s2_e2;

ktData_wat_e1_sum = smoothdata(sum(wat_e1(:,Ny_s,Nz_s,:,:),2), 1, 'gaussian', 15);
ktData_wat_e2_sum = smoothdata(sum(wat_e2(:,Ny_s,Nz_s,:,:),2), 1, 'gaussian', 15);

targetIdx = Nx+1;

wat_e1_sh = zeros(size(wat_e1), 'like', wat_e1);
wat_e2_sh = zeros(size(wat_e2), 'like', wat_e2);
met_s1_e1_sh = zeros(size(met_s1_e1), 'like', met_s1_e1);
met_s1_e2_sh = zeros(size(met_s1_e2), 'like', met_s1_e2);
met_s2_e1_sh = zeros(size(met_s2_e1), 'like', met_s2_e1);
met_s2_e2_sh = zeros(size(met_s2_e2), 'like', met_s2_e2);
maxIdx1_all = zeros(1, Nc);
maxIdx2_all = zeros(1, Nc);
shift1_all = zeros(1, Nc);
shift2_all = zeros(1, Nc);

if runParallel
    pool = gcp('nocreate');
    runParallel = ~isempty(pool);
end

if runParallel
    updateTsEPSIProgress(progress, 0.30, 'Echo alignment', 0, 'Detecting echo shifts');
    for nc = 1:Nc
        [maxIdx1, maxIdx2, sh1, sh2] = getEchoShift(ktData_wat_e1_sum, ktData_wat_e2_sum, targetIdx, nc);
        maxIdx1_all(nc) = maxIdx1;
        maxIdx2_all(nc) = maxIdx2;
        shift1_all(nc) = sh1;
        shift2_all(nc) = sh2;
        updateTsEPSIProgress(progress, 0.30 + 0.02*nc/Nc, 'Echo alignment', 0.15*nc/Nc, ...
            sprintf('Shift detection %d/%d', nc, Nc));
    end

    updateTsEPSIProgress(progress, 0.32, 'Echo alignment', 0.15, 'Starting parallel echo shifts');
    futures = parallel.FevalFuture.empty;
    nextCoil = 1;
    nQueued = min(Nc, pool.NumWorkers);
    for nq = 1:nQueued
        futures(end+1) = submitEchoShift(nextCoil); %#ok<AGROW>
        nextCoil = nextCoil + 1;
    end

    doneCount = 0;
    while doneCount < Nc
        [completedIdx, nc, wat_e1_c, wat_e2_c, met_s1_e1_c, met_s1_e2_c, met_s2_e1_c, met_s2_e2_c] = fetchNext(futures);

        wat_e1_sh(:,:,:,:,nc) = wat_e1_c;
        wat_e2_sh(:,:,:,:,nc) = wat_e2_c;
        met_s1_e1_sh(:,:,:,:,nc) = met_s1_e1_c;
        met_s1_e2_sh(:,:,:,:,nc) = met_s1_e2_c;
        met_s2_e1_sh(:,:,:,:,nc) = met_s2_e1_c;
        met_s2_e2_sh(:,:,:,:,nc) = met_s2_e2_c;

        doneCount = doneCount + 1;
        updateTsEPSIProgress(progress, 0.32 + 0.10*doneCount/Nc, 'Echo alignment', 0.15 + 0.85*doneCount/Nc, ...
            sprintf('Echo-shifted coils %d/%d', doneCount, Nc));

        if nextCoil <= Nc
            futures(completedIdx) = submitEchoShift(nextCoil);
            nextCoil = nextCoil + 1;
        else
            futures(completedIdx) = [];
        end
    end
else
    for nc = 1:Nc
        updateTsEPSIProgress(progress, 0.30 + 0.12*(nc-1)/Nc, 'Echo alignment', (nc-1)/Nc, ...
            sprintf('Coil %d/%d', nc, Nc));
        [maxIdx1, maxIdx2, sh1, sh2] = getEchoShift(ktData_wat_e1_sum, ktData_wat_e2_sum, targetIdx, nc);

        wat_e1_sh(:,:,:,:,nc) = circshift(wat_e1(:,:,:,:,nc), sh1, 1);
        wat_e2_sh(:,:,:,:,nc) = circshift(wat_e2(:,:,:,:,nc), sh2, 1);
        met_s1_e1_sh(:,:,:,:,nc) = circshift(met_s1_e1(:,:,:,:,nc), sh1, 1);
        met_s1_e2_sh(:,:,:,:,nc) = circshift(met_s1_e2(:,:,:,:,nc), sh2, 1);
        met_s2_e1_sh(:,:,:,:,nc) = circshift(met_s2_e1(:,:,:,:,nc), sh1, 1);
        met_s2_e2_sh(:,:,:,:,nc) = circshift(met_s2_e2(:,:,:,:,nc), sh2, 1);

        maxIdx1_all(nc) = maxIdx1;
        maxIdx2_all(nc) = maxIdx2;
    end
end
updateTsEPSIProgress(progress, 0.42, 'Echo alignment', 1, 'Echo shifts applied');

ktData_shifted.wat_e1 = wat_e1_sh;
ktData_shifted.wat_e2 = wat_e2_sh;
ktData_shifted.met_s1_e1 = met_s1_e1_sh;
ktData_shifted.met_s1_e2 = met_s1_e2_sh;
ktData_shifted.met_s2_e1 = met_s2_e1_sh;
ktData_shifted.met_s2_e2 = met_s2_e2_sh;

echoInfo.maxIdx1 = maxIdx1_all;
echoInfo.maxIdx2 = maxIdx2_all;

    function future = submitEchoShift(coilIndex)
        future = parfeval(pool, @shiftOneEchoCoil, 7, coilIndex, shift1_all(coilIndex), shift2_all(coilIndex), ...
            wat_e1(:,:,:,:,coilIndex), wat_e2(:,:,:,:,coilIndex), ...
            met_s1_e1(:,:,:,:,coilIndex), met_s1_e2(:,:,:,:,coilIndex), ...
            met_s2_e1(:,:,:,:,coilIndex), met_s2_e2(:,:,:,:,coilIndex));
    end
end

function [maxIdx1, maxIdx2, sh1, sh2] = getEchoShift(ktData_wat_e1_sum, ktData_wat_e2_sum, targetIdx, nc)
temp1 = squeeze(ktData_wat_e1_sum(:,1,1,1,nc));
temp2 = squeeze(ktData_wat_e2_sum(:,1,1,1,nc));
[~, maxIdx1] = max(abs(temp1));
[~, maxIdx2] = max(abs(temp2));
sh1 = targetIdx - maxIdx1;
sh2 = targetIdx - maxIdx2;
end

function [coilIndex, wat_e1_c, wat_e2_c, met_s1_e1_c, met_s1_e2_c, met_s2_e1_c, met_s2_e2_c] = ...
    shiftOneEchoCoil(coilIndex, sh1, sh2, wat_e1_c, wat_e2_c, met_s1_e1_c, met_s1_e2_c, met_s2_e1_c, met_s2_e2_c)
wat_e1_c = circshift(wat_e1_c, sh1, 1);
wat_e2_c = circshift(wat_e2_c, sh2, 1);
met_s1_e1_c = circshift(met_s1_e1_c, sh1, 1);
met_s1_e2_c = circshift(met_s1_e2_c, sh2, 1);
met_s2_e1_c = circshift(met_s2_e1_c, sh1, 1);
met_s2_e2_c = circshift(met_s2_e2_c, sh2, 1);
end
