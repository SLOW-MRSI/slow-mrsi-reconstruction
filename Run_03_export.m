clc
clear

%% User settings

% Select spectra source for SID conversion:
%   'raw'      -> load xfData_raw.mat variables xfData_ful/xfData_par/xfData_wat
%   'filtered' -> load xfData_Filt.mat struct xfData_Filt.ful/par/dif/wat
xfDataSource = 'filtered';

% SpectrIm uses ImagePositionPatient as the outer grid corner rather than
% the DICOM voxel-center convention. Values are [x y z] in half-voxel
% units; positive values move the first SID voxel along increasing index.
spectrImHalfVoxelShiftXYZ = [-1 0 0];

processedDataFolder = fullfile(pwd, 'processedData');
if ~exist(processedDataFolder, 'dir')
    mkdir(processedDataFolder)
end
spectrImOutputFolder = fullfile(pwd, 'data4spectrIm');
if ~exist(spectrImOutputFolder, 'dir')
    mkdir(spectrImOutputFolder)
end
xfDataRawFile = fullfile(processedDataFolder, 'xfData_raw.mat');
xfDataFilteredFile = fullfile(processedDataFolder, 'xfData_Filt.mat');
sidSeries = {'ful', 'par', 'dif', 'wat'};

%% Load spectral data and Twix geometry

xfDataSid = loadXfDataForSid(xfDataSource, xfDataRawFile, xfDataFilteredFile);

twixData = load(xfDataRawFile, 'twix_head');
assertStructFields(twixData, {'twix_head'}, xfDataRawFile)
twix_head = twixData.twix_head;

% Set true when the dataset's x/y labels are switched.
% false: rows=readout/base matrix, columns=phase matrix.
% true:  rows=phase matrix, columns=readout/base matrix.
switchXY = true;

spatialInfo = extractTwixSpatialInfo(twix_head);
if switchXY
    nominalXMatrix = spatialInfo.summary.kSpace.phoenix.phaseEncodingLines;
else
    nominalXMatrix = spatialInfo.summary.kSpace.phoenix.baseResolution;
end
xfDataSid = applyXOversamplingOption(xfDataSid, sidSeries, ...
    'remove', nominalXMatrix);
printXfDataSidSummary(xfDataSid, sidSeries)

exportSpatialSize = size(xfDataSid.ful);
dummyDicomOutputFolder = fullfile(spectrImOutputFolder, ...
    'MRSIDummyDicom');
spatialInfo.dummyDicom = buildDummyDicomGeometryFromTwix( ...
    spatialInfo, switchXY, exportSpatialSize(1:3), 'remove', ...
    spectrImHalfVoxelShiftXYZ);
writeZeroDicomSeriesFromTwix(spatialInfo.dummyDicom, dummyDicomOutputFolder)
sidOutputFolder = fullfile(spectrImOutputFolder, sprintf('xfData_%s', ...
    xfDataSid.sourceMode));
writeXfDataSidSeries(xfDataSid, sidSeries, sidOutputFolder)

printDummyDicomGeometry(spatialInfo.dummyDicom, dummyDicomOutputFolder)
fprintf('\nSID files are in:\n%s\n', sidOutputFolder)

fprintf('\nFound %d spatial header candidates in spatialInfo.headerCandidates.\n', ...
    height(spatialInfo.headerCandidates))
fprintf('Use disp(spatialInfo.headerCandidates) to inspect the full list.\n')
fprintf('Use spatialInfo.dummyDicom.sliceTable for the Twix-derived DICOM geometry.\n')


resonanceFreq = floor(xfDataSid.B0*1e6/1e4)*1e4;
fprintf('\n--------------------------------------User guide for spectrIm--------------------------------------\n')
fprintf('Please enter this "Resonance frequency" = %.0f when loading the MRSI data using spectrIm.\n', ...
    resonanceFreq)
fprintf('Please enter this "Sampling frequency" = %d when loading the MRSI data using spectrIm.\n', ...
    xfDataSid.BW)

%% functions
function xfDataSid = loadXfDataForSid(sourceMode, rawFile, filteredFile)
% Load raw or filtered spectra into a consistent ful/par/dif/wat struct.
    switch lower(sourceMode)
        case {'raw', 'xfdata_raw'}
            inputFile = rawFile;
            data = load(inputFile);
            requiredVars = {'xfData_ful', 'xfData_par', 'xfData_wat'};
            assertStructFields(data, requiredVars, inputFile)

            xfDataSid = struct();
            xfDataSid.sourceMode = 'raw';
            xfDataSid.sourceFile = inputFile;
            xfDataSid.ful = data.xfData_ful;
            xfDataSid.par = data.xfData_par;
            xfDataSid.dif = data.xfData_par - data.xfData_ful;
            xfDataSid.wat = data.xfData_wat;

        case {'filtered', 'filt', 'xfdata_filt'}
            inputFile = filteredFile;
            data = load(inputFile, 'xfData_Filt');
            assertStructFields(data, {'xfData_Filt'}, inputFile)
            assertStructFields(data.xfData_Filt, {'ful', 'par', 'dif', 'wat'}, ...
                'xfData_Filt')

            xfDataSid = data.xfData_Filt;
            xfDataSid.sourceMode = 'filtered';
            xfDataSid.sourceFile = inputFile;

        otherwise
            error('Unknown xfDataSource "%s". Use "raw" or "filtered".', sourceMode);
    end

    validateXfDataSidSeries(xfDataSid)
end

function validateXfDataSidSeries(xfDataSid)
    seriesNames = {'ful', 'par', 'dif', 'wat'};
    referenceSize = size(xfDataSid.ful);
    for k = 1:numel(seriesNames)
        seriesName = seriesNames{k};
        if ndims(xfDataSid.(seriesName)) ~= 4
            error('xfDataSid.%s must be a 4-D array [Nx Ny Nz Npts].', seriesName);
        end
        if ~isequal(size(xfDataSid.(seriesName)), referenceSize)
            error('xfDataSid.%s size [%s] does not match ful size [%s].', ...
                seriesName, num2str(size(xfDataSid.(seriesName))), ...
                num2str(referenceSize));
        end
    end
end

function xfDataSid = applyXOversamplingOption(xfDataSid, seriesNames, ...
    option, nominalXMatrix)
% Keep the acquired x grid or center-crop it to the nominal Twix grid.
    option = validatestring(option, {'keep', 'remove'});
    nominalXMatrix = round(nominalXMatrix);
    if ~isfinite(nominalXMatrix) || nominalXMatrix < 1
        error('A valid Twix nominal x matrix is required for x oversampling export.')
    end

    inputNx = size(xfDataSid.ful, 1);
    switch option
        case 'keep'
            if inputNx <= nominalXMatrix
                error(['xOversampling is ''keep'', but the loaded data have ' ...
                    '%d x voxels (Twix base matrix %d). Oversampling cannot be ' ...
                    'restored after it was removed; export data reconstructed at 96 voxels.'], ...
                    inputNx, nominalXMatrix)
            end
        case 'remove'
            if inputNx < nominalXMatrix
                error('Loaded x dimension %d is smaller than Twix nominal x matrix %d.', ...
                    inputNx, nominalXMatrix)
            elseif inputNx > nominalXMatrix
                firstX = floor((inputNx - nominalXMatrix) / 2) + 1;
                lastX = firstX + nominalXMatrix - 1;
                for k = 1:numel(seriesNames)
                    name = seriesNames{k};
                    xfDataSid.(name) = xfDataSid.(name)(firstX:lastX,:,:,:);
                end
                fprintf('Removed x oversampling by center crop: %d -> %d voxels.\n', ...
                    inputNx, nominalXMatrix)
            end
    end
    xfDataSid.xOversampling = option;
end

function printXfDataSidSummary(xfDataSid, seriesNames)
    fprintf('SID spectra source: %s (%s)\n', xfDataSid.sourceMode, ...
        xfDataSid.sourceFile)
    for k = 1:numel(seriesNames)
        seriesName = seriesNames{k};
        dataSize = size(xfDataSid.(seriesName));
        fprintf('  %-3s: %s [%s]\n', seriesName, class(xfDataSid.(seriesName)), ...
            num2str(dataSize))
    end
end

function assertStructFields(value, fieldNames, sourceName)
    for k = 1:numel(fieldNames)
        if ~isstruct(value) || ~isfield(value, fieldNames{k})
            error('%s is missing required field/variable "%s".', sourceName, ...
                fieldNames{k});
        end
    end
end

function writeXfDataSidSeries(xfDataSid, seriesNames, outputFolder)
% Write selected 4-D complex spectra as MIDAS/SpectrIm .sid float32 files.
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder)
    end

    for k = 1:numel(seriesNames)
        seriesName = seriesNames{k};
        outputFile = fullfile(outputFolder, sprintf('Matlab_%s.sid', seriesName));
        fprintf('Writing %s [%s] -> %s\n', seriesName, ...
            num2str(size(xfDataSid.(seriesName))), outputFile)
        writeSidComplexFloat32(xfDataSid.(seriesName), outputFile)
    end
end

function writeSidComplexFloat32(data, outputFile)
% SID order follows Process_3: z -> x -> y -> spectralPoint, with real/imag
% interleaved as float32.
    if ndims(data) ~= 4
        error('SID input must be 4-D: x, y, z, spectralPoint.')
    end

    matrixSize = size(data);
    nX = matrixSize(1);
    nY = matrixSize(2);
    nZ = matrixSize(3);
    nSpectral = matrixSize(4);

    fileID = fopen(outputFile, 'w');
    if fileID < 0
        error('Could not open %s for writing.', outputFile)
    end
    cleanupObj = onCleanup(@() fclose(fileID));

    interleaved = zeros(2*nSpectral, 1, 'single');
    for z = 1:nZ
        for xIndex = 1:nX
            for yIndex = 1:nY
                spectrum = squeeze(data(xIndex, yIndex, z, :));
                interleaved(1:2:end) = flip(single(real(spectrum))); % note it is flip in frequency domain to fit spectrIm
                interleaved(2:2:end) = -flip(single(imag(spectrum)));
                fwrite(fileID, interleaved, 'float32');
            end
        end
    end

    expectedBytes = 2*nSpectral*nX*nY*nZ*4;
    fileInfo = dir(outputFile);
    if isempty(fileInfo) || fileInfo.bytes ~= expectedBytes
        error('Unexpected file size for %s. Expected %d bytes.', ...
            outputFile, expectedBytes)
    end
end

function spatialInfo = extractTwixSpatialInfo(twixHead)
% Extract FOV, position, normal, rotation, and related spatial header fields.
    spatialInfo = struct();
    spatialInfo.summary = struct();

    spatialInfo.summary.config = getFlatTwixGeometry(twixHead, {'Config'});
    spatialInfo.summary.meas = getFlatTwixGeometry(twixHead, {'Meas'});
    spatialInfo.summary.dicom = getFlatTwixGeometry(twixHead, {'Dicom'});

    spatialInfo.summary.phoenixSlice = getNestedGeometry( ...
        twixHead, {'Phoenix', 'sSliceArray', 'asSlice', 1});
    spatialInfo.summary.phoenixVoi = getNestedGeometry( ...
        twixHead, {'Phoenix', 'sSpecPara', 'sVoI'});
    spatialInfo.summary.measYapsSlice = getNestedGeometry( ...
        twixHead, {'MeasYaps', 'sSliceArray', 'asSlice', 1});
    spatialInfo.summary.measYapsVoi = getNestedGeometry( ...
        twixHead, {'MeasYaps', 'sSpecPara', 'sVoI'});
    spatialInfo.summary.kSpace = getKSpaceGeometry(twixHead);

    spatialInfo.headerCandidates = findSpatialHeaderCandidates(twixHead);
end

function printDummyDicomGeometry(dummyDicom, outputFolder)
    fprintf('\nDummy DICOM series\n')
    fprintf('  output folder: %s\n', outputFolder)
    fprintf('  x/y switched:          %d\n', dummyDicom.switchXY)
    fprintf('  SpectrIm x/y/z shift:   [%s %s %s] half voxels\n', ...
        valueText(dummyDicom.spectrImHalfVoxelShiftXYZ(1)), ...
        valueText(dummyDicom.spectrImHalfVoxelShiftXYZ(2)), ...
        valueText(dummyDicom.spectrImHalfVoxelShiftXYZ(3)))
    fprintf('  matrix rows x columns: %d x %d\n', ...
        dummyDicom.rows, dummyDicom.columns)
    fprintf('  number of slices:      %d\n', dummyDicom.nSlices)
    fprintf('  pixel spacing row/col: [%s %s] mm\n', ...
        valueText(dummyDicom.pixelSpacing_RowColumn_mm(1)), ...
        valueText(dummyDicom.pixelSpacing_RowColumn_mm(2)))
    fprintf('  FOV row/column:        [%s %s] mm\n', ...
        valueText(dummyDicom.fov_RowColumn_mm(1)), ...
        valueText(dummyDicom.fov_RowColumn_mm(2)))
    fprintf('  slice thickness:       %s mm\n', ...
        valueText(dummyDicom.sliceThickness_mm))
    fprintf('  spacing between slices: %s mm\n', ...
        valueText(dummyDicom.spacingBetweenSlices_mm))
    fprintf('  Twix volume center Sag/Cor/Tra: ')
    printNumericVector(dummyDicom.volumeCenter_SagCorTra_mm, 'mm')
    fprintf('  DICOM orientation row:  ')
    printNumericVector(dummyDicom.rowDirection_XYZ, '')
    fprintf('  DICOM orientation col:  ')
    printNumericVector(dummyDicom.columnDirection_XYZ, '')
    fprintf('  DICOM normal:           ')
    printNumericVector(dummyDicom.normal_XYZ, '')
    fprintf('  tilt from transverse:   %s deg\n', ...
        valueText(dummyDicom.tiltFromTransverse_deg))

    fprintf('\nTwix-calculated dummy DICOM slice geometry\n')
    disp(dummyDicom.sliceTable(:, {'SliceIndex', 'InstanceNumber', ...
        'SliceLocation_mm', 'ImagePosition_X_mm', 'ImagePosition_Y_mm', ...
        'ImagePosition_Z_mm', 'Rows', 'Columns'}))
end

function dummyDicom = buildDummyDicomGeometryFromTwix(spatialInfo, switchXY, ...
    dataSpatialSize, xOversampling, spectrImHalfVoxelShiftXYZ)
% Calculate DICOM spatial metadata from Twix header values only.
    if nargin < 2
        switchXY = false;
    end
    if nargin < 5 || isempty(spectrImHalfVoxelShiftXYZ)
        spectrImHalfVoxelShiftXYZ = [0 0 0];
    end
    validateattributes(spectrImHalfVoxelShiftXYZ, {'numeric'}, ...
        {'real', 'finite', 'vector', 'numel', 3})
    spectrImHalfVoxelShiftXYZ = double(reshape(spectrImHalfVoxelShiftXYZ, 1, []));

    voi = spatialInfo.summary.phoenixVoi;
    slicePrescription = spatialInfo.summary.phoenixSlice;
    kSpace = spatialInfo.summary.kSpace.phoenix;

    nominalReadoutMatrix = round(kSpace.baseResolution);
    nominalPhaseMatrix = round(kSpace.phaseEncodingLines);
    if switchXY
        phaseMatrix = dataSpatialSize(1);
        readoutMatrix = dataSpatialSize(2);
    else
        readoutMatrix = dataSpatialSize(1);
        phaseMatrix = dataSpatialSize(2);
    end
    nSlices = dataSpatialSize(3);

    nominalXMatrix = nominalReadoutMatrix;
    if switchXY
        nominalXMatrix = nominalPhaseMatrix;
    end
    if strcmpi(xOversampling, 'remove') && dataSpatialSize(1) ~= nominalXMatrix
        error('Removed-oversampling export must have %d x voxels, but data have %d.', ...
            nominalXMatrix, dataSpatialSize(1))
    end

    % Oversampling adds x samples outside the nominal FOV while preserving
    % voxel size. With switchXY=true, x is the phase-labelled direction.
    readoutFactor = readoutMatrix / nominalReadoutMatrix;
    phaseFactor = phaseMatrix / nominalPhaseMatrix;
    readFov = voi.readFov_mm * readoutFactor;
    phaseFov = voi.phaseFov_mm * phaseFactor;

    if switchXY
        rows = phaseMatrix;
        columns = readoutMatrix;
        rowFov = phaseFov;
        columnFov = readFov;
    else
        rows = readoutMatrix;
        columns = phaseMatrix;
        rowFov = readFov;
        columnFov = phaseFov;
    end

    rowSpacing = rowFov / rows;
    columnSpacing = columnFov / columns;

    volumeCenter = voi.position_SagCorTra_mm;
    volumeCenter(~isfinite(volumeCenter)) = 0;

    normal = voi.normal_SagCorTra;
    normal(~isfinite(normal)) = 0;
    if norm(normal) == 0
        normal = [0 0 1];
    else
        normal = normal ./ norm(normal);
    end

    [rowDirection, columnDirection] = patientAxesFromNormal(normal);

    slabThickness = slicePrescription.sliceThickness_mm;
    if ~isfinite(slabThickness)
        slabThickness = voi.sliceThickness_mm;
    end
    sliceSpacing = slabThickness / nSlices;
    sliceOffsets = ((1:nSlices) - (nSlices + 1) / 2) .* sliceSpacing;

    records = repmat(emptyDummyDicomSliceRecord(), nSlices, 1);
    for s = 1:nSlices
        sliceCenter = volumeCenter + normal .* sliceOffsets(s);
        % SpectrIm interprets ImagePositionPatient as the outer grid corner.
        % The optional [x y z] adjustment is in half-voxel units. x is the
        % first SID dimension, y is the second, and z follows slice index.
        imagePosition = sliceCenter ...
            - rowDirection .* (columnFov / 2) ...
            - columnDirection .* (rowFov / 2) ...
            + columnDirection .* (spectrImHalfVoxelShiftXYZ(1) * rowSpacing / 2) ...
            + rowDirection .* (spectrImHalfVoxelShiftXYZ(2) * columnSpacing / 2) ...
            + normal .* (spectrImHalfVoxelShiftXYZ(3) * sliceSpacing / 2);

        records(s).SliceIndex = s;
        records(s).InstanceNumber = s;
        records(s).SliceLocation_mm = dot(sliceCenter, normal);
        records(s).ImagePosition_X_mm = imagePosition(1);
        records(s).ImagePosition_Y_mm = imagePosition(2);
        records(s).ImagePosition_Z_mm = imagePosition(3);
        records(s).Center_Sag_mm = sliceCenter(1);
        records(s).Center_Cor_mm = sliceCenter(2);
        records(s).Center_Tra_mm = sliceCenter(3);
        records(s).Rows = rows;
        records(s).Columns = columns;
    end

    dummyDicom = struct();
    dummyDicom.switchXY = switchXY;
    dummyDicom.spectrImHalfVoxelShiftXYZ = spectrImHalfVoxelShiftXYZ;
    dummyDicom.nominalReadoutMatrix = nominalReadoutMatrix;
    dummyDicom.nominalPhaseMatrix = nominalPhaseMatrix;
    dummyDicom.readoutOversamplingFactor = readoutFactor;
    dummyDicom.phaseOversamplingFactor = phaseFactor;
    dummyDicom.readoutMatrix = readoutMatrix;
    dummyDicom.phaseMatrix = phaseMatrix;
    dummyDicom.rows = rows;
    dummyDicom.columns = columns;
    dummyDicom.nSlices = nSlices;
    dummyDicom.fov_RowColumn_mm = [rowFov columnFov];
    dummyDicom.pixelSpacing_RowColumn_mm = [rowSpacing columnSpacing];
    dummyDicom.sliceThickness_mm = sliceSpacing;
    dummyDicom.spacingBetweenSlices_mm = sliceSpacing;
    dummyDicom.volumeCenter_SagCorTra_mm = volumeCenter;
    dummyDicom.rowDirection_XYZ = rowDirection;
    dummyDicom.columnDirection_XYZ = columnDirection;
    dummyDicom.normal_XYZ = normal;
    dummyDicom.tiltFromTransverse_deg = acosd(max(min(abs(normal(3)), 1), -1));
    dummyDicom.sliceTable = struct2table(records);
end

function [rowDirection, columnDirection] = patientAxesFromNormal(normal)
    normal = normal ./ norm(normal);

    readoutReference = [1 0 0];
    rowDirection = readoutReference - dot(readoutReference, normal) .* normal;

    if norm(rowDirection) < 1e-6
        readoutReference = [0 1 0];
        rowDirection = readoutReference - dot(readoutReference, normal) .* normal;
    end
    rowDirection = rowDirection ./ norm(rowDirection);

    columnDirection = cross(normal, rowDirection);
    columnDirection = columnDirection ./ norm(columnDirection);
end

function writeZeroDicomSeriesFromTwix(dummyDicom, outputFolder)
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder)
    end

    studyUID = dicomuid;
    seriesUID = dicomuid;
    frameUID = dicomuid;
    zeroImage = zeros(dummyDicom.rows, dummyDicom.columns, 'uint16');
    for s = 1:dummyDicom.nSlices
        row = dummyDicom.sliceTable(s, :);
        metadata = makeTwixDummyDicomMetadata(dummyDicom, row, studyUID, ...
            seriesUID, frameUID);
        outputFile = fullfile(outputFolder, sprintf('dummy_MRSI_slice_%03d.dcm', s));
        dicomwrite(zeroImage, outputFile, metadata, 'CreateMode', 'Create')
    end
end

function metadata = makeTwixDummyDicomMetadata(dummyDicom, row, studyUID, ...
    seriesUID, frameUID)
        metadata = struct();
        metadata.PatientName = 'TwixDummy';
        metadata.PatientID = 'TwixDummy';
        metadata.PatientBirthDate = '19000101';
        metadata.PatientSex = 'O';
        metadata.Modality = 'MR';
        metadata.Manufacturer = 'SIEMENS';
        metadata.StudyID = '1';
        metadata.SeriesNumber = 1;
        metadata.AcquisitionNumber = 1;
        metadata.ImagesInAcquisition = dummyDicom.nSlices;
        metadata.StudyDate = '20260706';
        metadata.StudyTime = '000000';
        metadata.SeriesDate = metadata.StudyDate;
        metadata.SeriesTime = metadata.StudyTime;
        metadata.AcquisitionDate = metadata.StudyDate;
        metadata.AcquisitionTime = metadata.StudyTime;
        metadata.ContentDate = metadata.StudyDate;
        metadata.ContentTime = metadata.StudyTime;
        metadata.SeriesDescription = 'Dummy MRSI spatial reference';
        metadata.ProtocolName = 'TwixDummyMRSI';
        metadata.ImageType = 'DERIVED\SECONDARY\DUMMY';
        metadata.StudyInstanceUID = studyUID;
        metadata.SeriesInstanceUID = seriesUID;
        metadata.FrameOfReferenceUID = frameUID;
        metadata.SOPInstanceUID = dicomuid;
        metadata.InstanceNumber = row.InstanceNumber;
        metadata.Rows = dummyDicom.rows;
        metadata.Columns = dummyDicom.columns;
        metadata.PixelSpacing = flip(dummyDicom.pixelSpacing_RowColumn_mm(:));
        metadata.SliceThickness = dummyDicom.sliceThickness_mm;
        metadata.SpacingBetweenSlices = dummyDicom.spacingBetweenSlices_mm;
        metadata.ImageOrientationPatient = [dummyDicom.rowDirection_XYZ(:); ...
            dummyDicom.columnDirection_XYZ(:)];
        metadata.ImagePositionPatient = [row.ImagePosition_X_mm; ...
            row.ImagePosition_Y_mm; row.ImagePosition_Z_mm];
        metadata.SliceLocation = row.SliceLocation_mm;
        metadata.SamplesPerPixel = 1;
        metadata.PhotometricInterpretation = 'MONOCHROME2';
        metadata.BitsAllocated = 16;
        metadata.BitsStored = 16;
        metadata.HighBit = 15;
        metadata.PixelRepresentation = 0;
end

function record = emptyDummyDicomSliceRecord()
    record = struct();
    record.SliceIndex = NaN;
    record.InstanceNumber = NaN;
    record.SliceLocation_mm = NaN;
    record.ImagePosition_X_mm = NaN;
    record.ImagePosition_Y_mm = NaN;
    record.ImagePosition_Z_mm = NaN;
    record.Center_Sag_mm = NaN;
    record.Center_Cor_mm = NaN;
    record.Center_Tra_mm = NaN;
    record.Rows = NaN;
    record.Columns = NaN;
end

function printNumericVector(values, unitText)
    fprintf('[%s %s %s]', valueText(values(1)), valueText(values(2)), ...
        valueText(values(3)))
    if ~isempty(unitText)
        fprintf(' %s', unitText)
    end
    fprintf('\n')
end

function printSagCorCenterWarning(volumeCenter_SagCorTra_mm)
    sagCorTolerance_mm = 1e-6;
    sagCorCenter_mm = volumeCenter_SagCorTra_mm(1:2);
    if any(isfinite(sagCorCenter_mm) & abs(sagCorCenter_mm) > sagCorTolerance_mm)
        fprintf(['\n[WARNING]: Twix volume center Sag/Cor is [%s %s] mm, ' ...
            'not [0 0] mm. \nCorrect acquisition of this sequence requires ' ...
            'the FOV location in Sag and Cor to be 0 mm.\n'], ...
            valueText(sagCorCenter_mm(1)), valueText(sagCorCenter_mm(2)))
    end
end

function text = valueText(value)
    if isnumeric(value) && isscalar(value) && isnan(value)
        text = 'n/a';
    elseif isnumeric(value) || islogical(value)
        text = num2str(value, '%.12g');
    else
        text = char(string(value));
    end
end

function geom = getFlatTwixGeometry(twixHead, branchPath)
    branch = getByPath(twixHead, branchPath);
    geom = emptyGeometryStruct();
    geom.source = strjoin(branchPath, '.');

    if isempty(branch) || ~isstruct(branch)
        return
    end

    geom.readFov_mm = firstFieldValue(branch, {'RoFOV', 'ReadFoV', ...
        'VoiReadoutFOV', 'VoI_RoFOV', 'Adj_RoFOV', 'dReadoutFOV', ...
        'dReadFOE', 'dTotalFOVUser_mm', 'readout_fov'});
    geom.phaseFov_mm = firstFieldValue(branch, {'PeFOV', 'PhaseFoV', ...
        'VoiPhaseFOV', 'VoI_PeFOV', 'Adj_PeFOV', 'dPhaseFOV', ...
        'dPhaseFOE'});
    geom.sliceThickness_mm = firstFieldValue(branch, {'VoI_SliceThickness', ...
        'VoiThickness', 'Adj_SliceThickness', 'SliceThickness', ...
        'dThickness', 'thickness', 'sl_thick'});
    geom.inPlaneRotation_rad = firstFieldValue(branch, ...
        {'VoI_InPlaneRotAngle', 'VoiInPlaneRot', 'Adj_InPlaneRotAngle', ...
        'dInPlaneRot'});

    geom.position_SagCorTra_mm = [ ...
        firstFieldValue(branch, {'VoI_Position_Sag', 'VoiPositionSag', ...
            'Adj_Position_Sag', 'Position_Sag', 'dSag'}), ...
        firstFieldValue(branch, {'VoI_Position_Cor', 'VoiPositionCor', ...
            'Adj_Position_Cor', 'Position_Cor', 'dCor'}), ...
        firstFieldValue(branch, {'VoI_Position_Tra', 'VoiPositionTra', ...
            'Adj_Position_Tra', 'Position_Tra', 'dTra'})];

    geom.normal_SagCorTra = [ ...
        firstFieldValue(branch, {'VoI_Normal_Sag', 'VoiNormalSag', ...
            'Adj_Normal_Sag', 'Normal_Sag'}), ...
        firstFieldValue(branch, {'VoI_Normal_Cor', 'VoiNormalCor', ...
            'Adj_Normal_Cor', 'Normal_Cor'}), ...
        firstFieldValue(branch, {'VoI_Normal_Tra', 'VoiNormalTra', ...
            'Adj_Normal_Tra', 'Normal_Tra'})];

    geom.origin_XYZ_mm = [ ...
        firstFieldValue(branch, {'SBCSOriginPositionX', 'lSBCSOriginPositionX'}), ...
        firstFieldValue(branch, {'SBCSOriginPositionY', 'lSBCSOriginPositionY'}), ...
        firstFieldValue(branch, {'SBCSOriginPositionZ', 'lSBCSOriginPositionZ'})];

    geom.tablePosition_SagCorTra_mm = [ ...
        firstFieldValue(branch, {'GlobalTablePosSag', 'lGlobalTablePosSag'}), ...
        firstFieldValue(branch, {'GlobalTablePosCor', 'lGlobalTablePosCor'}), ...
        firstFieldValue(branch, {'GlobalTablePosTra', 'lGlobalTablePosTra'})];
end

function geom = getNestedGeometry(twixHead, blockPath)
    block = getByPath(twixHead, blockPath);
    geom = emptyGeometryStruct();
    geom.source = pathToText(blockPath);

    if isempty(block) || ~isstruct(block)
        return
    end

    geom.readFov_mm = firstFieldValue(block, {'dReadoutFOV', 'dReadFOE', ...
        'ReadFoV', 'RoFOV'});
    geom.phaseFov_mm = firstFieldValue(block, {'dPhaseFOV', 'dPhaseFOE', ...
        'PhaseFoV', 'PeFOV'});
    geom.sliceThickness_mm = firstFieldValue(block, {'dThickness', ...
        'thickness', 'SliceThickness'});
    geom.inPlaneRotation_rad = firstFieldValue(block, {'dInPlaneRot', ...
        'InPlaneRotAngle'});

    pos = getStructField(block, 'sPosition');
    geom.position_SagCorTra_mm = [ ...
        firstFieldValue(pos, {'dSag'}), ...
        firstFieldValue(pos, {'dCor'}), ...
        firstFieldValue(pos, {'dTra'})];

    normal = getStructField(block, 'sNormal');
    geom.normal_SagCorTra = [ ...
        firstFieldValue(normal, {'dSag'}), ...
        firstFieldValue(normal, {'dCor'}), ...
        firstFieldValue(normal, {'dTra'})];
end

function kSpace = getKSpaceGeometry(twixHead)
    kSpace = struct();
    kSpace.phoenix = getKSpaceBlock(twixHead, {'Phoenix', 'sKSpace'});
    kSpace.measYaps = getKSpaceBlock(twixHead, {'MeasYaps', 'sKSpace'});
end

function kSpace = getKSpaceBlock(twixHead, blockPath)
    block = getByPath(twixHead, blockPath);
    kSpace = struct();
    kSpace.source = pathToText(blockPath);

    if isempty(block) || ~isstruct(block)
        return
    end

    kSpace.baseResolution = firstFieldValue(block, {'lBaseResolution', ...
        'BaseResolution'});
    kSpace.phaseEncodingLines = firstFieldValue(block, ...
        {'lPhaseEncodingLines', 'PhaseEncodingLines'});
    kSpace.partitions = firstFieldValue(block, {'lPartitions', 'Partitions'});
    kSpace.imagesPerSlab = firstFieldValue(block, {'lImagesPerSlab'});
    kSpace.phaseResolution = firstFieldValue(block, {'dPhaseResolution'});
    kSpace.sliceResolution = firstFieldValue(block, {'dSliceResolution'});
    kSpace.dimension = firstFieldValue(block, {'ucDimension', 'KSpace_Dimension'});
    kSpace.trajectory = firstFieldValue(block, {'ucTrajectory'});
end

function geom = emptyGeometryStruct()
    geom = struct();
    geom.source = '';
    geom.readFov_mm = NaN;
    geom.phaseFov_mm = NaN;
    geom.sliceThickness_mm = NaN;
    geom.inPlaneRotation_rad = NaN;
    geom.position_SagCorTra_mm = [NaN NaN NaN];
    geom.normal_SagCorTra = [NaN NaN NaN];
    geom.origin_XYZ_mm = [NaN NaN NaN];
    geom.tablePosition_SagCorTra_mm = [NaN NaN NaN];
end

function candidates = findSpatialHeaderCandidates(value)
    keywords = {'fov', 'foe', 'position', 'normal', 'inplanerot', ...
        'rottrans', 'sag', 'cor', 'tra', 'slicearray', 'svoi', ...
        'thick', 'spatial', 'origin', 'tablepos', 'readoutfov', ...
        'phasefov', 'kspace'};

    rows = collectSpatialCandidates(value, '', keywords);
    if isempty(rows)
        candidates = cell2table(cell(0, 3), ...
            'VariableNames', {'Path', 'Class', 'Value'});
    else
        candidates = cell2table(rows, ...
            'VariableNames', {'Path', 'Class', 'Value'});
    end
end

function rows = collectSpatialCandidates(value, path, keywords)
    rows = {};

    if isstruct(value)
        fields = fieldnames(value);
        for idx = 1:numel(value)
            for f = 1:numel(fields)
                fieldName = fields{f};
                nextPath = appendPath(path, fieldName, idx, numel(value));
                fieldValue = value(idx).(fieldName);
                if isSpatialPath(nextPath, keywords) && isPrintableScalar(fieldValue)
                    rows(end + 1, :) = {nextPath, class(fieldValue), ...
                        valueToText(fieldValue)}; %#ok<AGROW>
                end
                rows = [rows; collectSpatialCandidates(fieldValue, ...
                    nextPath, keywords)]; %#ok<AGROW>
            end
        end
    elseif iscell(value)
        for idx = 1:numel(value)
            nextPath = sprintf('%s{%d}', path, idx);
            rows = [rows; collectSpatialCandidates(value{idx}, ...
                nextPath, keywords)]; %#ok<AGROW>
        end
    end
end

function tf = isSpatialPath(path, keywords)
    pathLower = lower(path);
    tf = false;
    for k = 1:numel(keywords)
        if ~isempty(strfind(pathLower, keywords{k}))
            tf = true;
            return
        end
    end
end

function tf = isPrintableScalar(value)
    tf = (isnumeric(value) || islogical(value) || ischar(value) || ...
        isstring(value)) && numel(value) <= 12;
end

function text = valueToText(value)
    if isnumeric(value) || islogical(value)
        text = mat2str(value);
    elseif ischar(value)
        text = value;
    else
        text = char(join(string(value), ', '));
    end
end

function path = appendPath(basePath, fieldName, idx, nItems)
    if isempty(basePath)
        path = fieldName;
    else
        path = [basePath '.' fieldName];
    end

    if nItems > 1
        path = sprintf('%s(%d)', path, idx);
    end
end

function value = firstFieldValue(source, fieldNames)
    value = NaN;

    if isempty(source) || ~isstruct(source)
        return
    end

    for f = 1:numel(fieldNames)
        fieldName = fieldNames{f};
        if isfield(source, fieldName)
            candidate = source.(fieldName);
            if isnumeric(candidate) && isscalar(candidate)
                value = candidate;
                return
            elseif isnumeric(candidate) && isempty(candidate)
                continue
            elseif islogical(candidate) && isscalar(candidate)
                value = double(candidate);
                return
            end
        end
    end
end

function value = getStructField(source, fieldName)
    value = [];
    if isstruct(source) && isfield(source, fieldName)
        value = source.(fieldName);
    end
end

function value = getByPath(source, pathParts)
    value = source;
    for p = 1:numel(pathParts)
        part = pathParts{p};
        if isnumeric(part)
            if iscell(value)
                if numel(value) < part
                    value = [];
                    return
                end
                value = value{part};
            elseif isstruct(value)
                if numel(value) < part
                    value = [];
                    return
                end
                value = value(part);
            else
                value = [];
                return
            end
        elseif isstruct(value) && isfield(value, part)
            value = value.(part);
        else
            value = [];
            return
        end
    end
end

function text = pathToText(pathParts)
    pieces = cell(1, numel(pathParts));
    nPieces = 0;
    for p = 1:numel(pathParts)
        part = pathParts{p};
        if isnumeric(part)
            if nPieces == 0
                nPieces = 1;
                pieces{nPieces} = sprintf('{%d}', part);
            else
                pieces{nPieces} = sprintf('%s{%d}', pieces{nPieces}, part);
            end
        else
            nPieces = nPieces + 1;
            pieces{nPieces} = part;
        end
    end
    text = strjoin(pieces(1:nPieces), '.');
end
