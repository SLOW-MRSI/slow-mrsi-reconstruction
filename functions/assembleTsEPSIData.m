function [ktData_met, ktData_wat, reconInfo] = assembleTsEPSIData(rawData_met, CsiInfo, twix_obj_met, cfg, progress, runParallel)
%ASSEMBLETSEPSIDATA Sort raw SLOW ts-EPSI data into metabolite and water arrays.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

Nx = CsiInfo.Ny;
Ny = CsiInfo.Nx;
Nz = CsiInfo.Nz;
Nt = twix_obj_met.hdr.MeasYaps.sSpecPara.lVectorSize*1.25;
Nc = CsiInfo.Nc;
Nkx = Nx*2;

if isempty(cfg.metaboliteTimePoints)
    Ntime = round(Nt/1.25);
else
    Ntime = cfg.metaboliteTimePoints;
end

Nz_s = twix_obj_met.image.Par(1:2:end);
Ny_s = twix_obj_met.image.Seg(1:2:end);
Na_s = twix_obj_met.image.Ave(1:2:end);

Nz_s = reshape(Nz_s,Nt,[]);
Nz_s = reshape(Nz_s(:,1:2:end),1,[]);
Ny_s = reshape(Ny_s,Nt,[]);
Ny_s = reshape(Ny_s(:,1:2:end),1,[]);
Na_s = reshape(Na_s,Nt,[]);
Na_s = reshape(Na_s(:,1:2:end),1,[]);

[routeMatrix, averageCount] = makeSortingRoute(Ny_s, Nz_s, Na_s, Ny, Nz, Nt, cfg.maxAverages, progress);
denominator = reshape(averageCount, [1, Ny, Nz, Nt]) + 1e-50;

waterTimeIdx = Ntime+1:Nt;
nWaterTime = min(numel(waterTimeIdx), Ntime);
waterTimeIdx = waterTimeIdx(1:nWaterTime);

ktData_met_s1_e1 = zeros(Nkx,Ny,Nz,Ntime,Nc, 'like', rawData_met);
ktData_met_s1_e2 = zeros(Nkx,Ny,Nz,Ntime,Nc, 'like', rawData_met);
ktData_met_s2_e1 = zeros(Nkx,Ny,Nz,Ntime,Nc, 'like', rawData_met);
ktData_met_s2_e2 = zeros(Nkx,Ny,Nz,Ntime,Nc, 'like', rawData_met);
ktData_wat_e1 = zeros(Nkx,Ny,Nz,Ntime,Nc, 'like', rawData_met);
ktData_wat_e2 = zeros(Nkx,Ny,Nz,Ntime,Nc, 'like', rawData_met);

if runParallel
    pool = gcp('nocreate');
    runParallel = ~isempty(pool);
end

if runParallel
    doneCount = 0;
    futures = parallel.FevalFuture.empty;
    nextCoil = 1;
    nQueued = min(Nc, pool.NumWorkers);

    updateTsEPSIProgress(progress, 0.12, 'Sorting raw data', 0, 'Starting parallel coil sorting');
    for nq = 1:nQueued
        futures(end+1) = submitSortCoil(nextCoil); %#ok<AGROW>
        nextCoil = nextCoil + 1;
    end

    while doneCount < Nc
        [completedIdx, nc, met_s1_e1, met_s1_e2, met_s2_e1, met_s2_e2, wat_e1, wat_e2] = fetchNext(futures);
        ktData_met_s1_e1(:,:,:,:,nc) = met_s1_e1;
        ktData_met_s1_e2(:,:,:,:,nc) = met_s1_e2;
        ktData_met_s2_e1(:,:,:,:,nc) = met_s2_e1;
        ktData_met_s2_e2(:,:,:,:,nc) = met_s2_e2;
        ktData_wat_e1(:,:,:,:,nc) = wat_e1;
        ktData_wat_e2(:,:,:,:,nc) = wat_e2;

        doneCount = doneCount + 1;
        updateTsEPSIProgress(progress, 0.12 + 0.18*doneCount/Nc, 'Sorting raw data', doneCount/Nc, ...
            sprintf('Coils complete %d/%d', doneCount, Nc));

        if nextCoil <= Nc
            futures(completedIdx) = submitSortCoil(nextCoil);
            nextCoil = nextCoil + 1;
        else
            futures(completedIdx) = [];
        end
    end
else
    for nc = 1:Nc
        updateTsEPSIProgress(progress, 0.12 + 0.18*(nc-1)/Nc, 'Sorting raw data', (nc-1)/Nc, ...
            sprintf('Coil %d/%d', nc, Nc));

        [ktData_met_s1_e1(:,:,:,:,nc), ktData_met_s1_e2(:,:,:,:,nc), ...
            ktData_met_s2_e1(:,:,:,:,nc), ktData_met_s2_e2(:,:,:,:,nc), ...
            ktData_wat_e1(:,:,:,:,nc), ktData_wat_e2(:,:,:,:,nc)] = ...
            sortOneCoil(rawData_met(:,nc,:), routeMatrix, denominator, waterTimeIdx, Ntime, Nkx, Ny, Nz, Nt);
    end
end

ktData_met.s1_e1 = ktData_met_s1_e1;
ktData_met.s1_e2 = ktData_met_s1_e2;
ktData_met.s2_e1 = ktData_met_s2_e1;
ktData_met.s2_e2 = ktData_met_s2_e2;
ktData_wat.e1 = ktData_wat_e1;
ktData_wat.e2 = ktData_wat_e2;
updateTsEPSIProgress(progress, 0.30, 'Sorting raw data', 1, 'Metabolite and embedded water data prepared');

reconInfo = struct();
reconInfo.Nx = Nx;
reconInfo.Nz = Nz;
reconInfo.Nt = Nt;
reconInfo.Nc = Nc;
reconInfo.Ntime = Ntime;
reconInfo.nWaterTime = nWaterTime;

    function future = submitSortCoil(coilIndex)
        future = parfeval(pool, @sortOneCoilWithIndex, 7, coilIndex, rawData_met(:,coilIndex,:), ...
            routeMatrix, denominator, waterTimeIdx, Ntime, Nkx, Ny, Nz, Nt);
    end
end

function [routeMatrix, averageCount] = makeSortingRoute(Ny_s, Nz_s, Na_s, Ny, Nz, Nt, maxAverages, progress)
updateTsEPSIProgress(progress, 0.08, 'Preparing sorting route', 0.05, 'Indexing acquired samples');
numSamples = numel(Ny_s);

updateTsEPSIProgress(progress, 0.09, 'Preparing sorting route', 0.30, 'Computing time indices');
timeIndex = repmat(1:Nt, 1, numSamples/Nt);

updateTsEPSIProgress(progress, 0.10, 'Preparing sorting route', 0.50, 'Finding valid averages');
valid = Na_s <= maxAverages;

updateTsEPSIProgress(progress, 0.11, 'Preparing sorting route', 0.70, 'Computing target positions');
sampleIndex = sub2ind([Ny, Nz, Nt], Ny_s(valid), Nz_s(valid), timeIndex(valid));

updateTsEPSIProgress(progress, 0.115, 'Preparing sorting route', 0.85, 'Building sparse sorting matrix');
routeMatrix = sparse(find(valid), sampleIndex, 1, numSamples, Ny*Nz*Nt);

updateTsEPSIProgress(progress, 0.118, 'Preparing sorting route', 0.95, 'Counting averages per voxel');
averageCount = full(sum(routeMatrix, 1));

updateTsEPSIProgress(progress, 0.12, 'Preparing sorting route', 1, 'Sorting route ready');
end

function [coilIndex, met_s1_e1, met_s1_e2, met_s2_e1, met_s2_e2, wat_e1, wat_e2] = sortOneCoilWithIndex(coilIndex, rawData_coil, routeMatrix, denominator, waterTimeIdx, Ntime, Nkx, Ny, Nz, Nt)
[met_s1_e1, met_s1_e2, met_s2_e1, met_s2_e2, wat_e1, wat_e2] = ...
    sortOneCoil(rawData_coil, routeMatrix, denominator, waterTimeIdx, Ntime, Nkx, Ny, Nz, Nt);
end

function [met_s1_e1, met_s1_e2, met_s2_e1, met_s2_e2, wat_e1, wat_e2] = sortOneCoil(rawData_coil, routeMatrix, denominator, waterTimeIdx, Ntime, Nkx, Ny, Nz, Nt)
rawData_coil = squeeze(rawData_coil);
rawData_1 = rawData_coil(:,1:2:end);
rawData_2 = flip(rawData_coil(:,2:2:end),1);

temp = reshape(rawData_1,Nkx,Nt,[]);
rawData_e1_s1 = reshape(temp(:,:,1:2:end),Nkx,[]);
rawData_e1_s2 = reshape(temp(:,:,2:2:end),Nkx,[]);

temp2 = reshape(rawData_2,Nkx,Nt,[]);
rawData_e2_s1 = reshape(temp2(:,:,1:2:end),Nkx,[]);
rawData_e2_s2 = reshape(temp2(:,:,2:2:end),Nkx,[]);

% MATLAB does not support single*sparse mtimes, so route in double and
% cast back immediately. The sparse route avoids the old slow scatter loop.
acc_e1_s1 = cast(reshape(double(rawData_e1_s1) * routeMatrix, Nkx, Ny, Nz, Nt), 'like', rawData_coil);
acc_e1_s2 = cast(reshape(double(rawData_e1_s2) * routeMatrix, Nkx, Ny, Nz, Nt), 'like', rawData_coil);
acc_e2_s1 = cast(reshape(double(rawData_e2_s1) * routeMatrix, Nkx, Ny, Nz, Nt), 'like', rawData_coil);
acc_e2_s2 = cast(reshape(double(rawData_e2_s2) * routeMatrix, Nkx, Ny, Nz, Nt), 'like', rawData_coil);

met_s1_e1 = acc_e1_s1(:,:,:,1:Ntime) ./ denominator(:,:,:,1:Ntime);
met_s1_e2 = acc_e2_s1(:,:,:,1:Ntime) ./ denominator(:,:,:,1:Ntime);
met_s2_e1 = acc_e1_s2(:,:,:,1:Ntime) ./ denominator(:,:,:,1:Ntime);
met_s2_e2 = acc_e2_s2(:,:,:,1:Ntime) ./ denominator(:,:,:,1:Ntime);

wat_e1 = zeros(Nkx,Ny,Nz,Ntime, 'like', rawData_coil);
wat_e2 = zeros(Nkx,Ny,Nz,Ntime, 'like', rawData_coil);
wat_e1(:,:,:,1:numel(waterTimeIdx)) = acc_e1_s1(:,:,:,waterTimeIdx) ./ denominator(:,:,:,waterTimeIdx);
wat_e2(:,:,:,1:numel(waterTimeIdx)) = acc_e2_s1(:,:,:,waterTimeIdx) ./ denominator(:,:,:,waterTimeIdx);
end
