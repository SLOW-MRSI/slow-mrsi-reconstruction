function [scannerFrequencyMHz, sourceField] = getScannerFrequencyMHzFromTwixHeader(twixHeader)
%GETSCANNERFREQUENCYMHZFROMTWIXHEADER Return scanner/Larmor frequency from a twix header.
% Siemens twix headers store lFrequency in Hz. The reconstruction code uses
% MHz for ppm-axis conversion.

scannerFrequencyMHz = [];
sourceField = '';

if nargin < 1 || isempty(twixHeader) || ~isstruct(twixHeader)
    return
end

frequencyCandidates = { ...
    {'Meas', 'lFrequency'}, ...
    {'Dicom', 'lFrequency'}, ...
    {'MeasYaps', 'sTXSPEC', 'asNucleusInfo', 1, 'lFrequency'}};

for idxCandidate = 1:numel(frequencyCandidates)
    pathParts = frequencyCandidates{idxCandidate};
    [frequencyValue, ok] = readNestedField(twixHeader, pathParts);
    if ok && isnumeric(frequencyValue) && isscalar(frequencyValue) && isfinite(frequencyValue) && frequencyValue > 0
        scannerFrequencyMHz = double(frequencyValue);
        if scannerFrequencyMHz > 1e6
            scannerFrequencyMHz = scannerFrequencyMHz / 1e6;
        end
        sourceField = formatPath(pathParts);
        return
    end
end
end

function [value, ok] = readNestedField(value, pathParts)
ok = true;
for idxPart = 1:numel(pathParts)
    part = pathParts{idxPart};
    if ischar(part) || isstring(part)
        part = char(part);
        if ~isstruct(value) || ~isfield(value, part)
            ok = false;
            value = [];
            return
        end
        value = value.(part);
    elseif isnumeric(part) && isscalar(part)
        if iscell(value)
            if numel(value) < part
                ok = false;
                value = [];
                return
            end
            value = value{part};
        else
            if numel(value) < part
                ok = false;
                value = [];
                return
            end
            value = value(part);
        end
    end
end
end

function pathText = formatPath(pathParts)
pathText = '';
for idxPart = 1:numel(pathParts)
    part = pathParts{idxPart};
    if ischar(part) || isstring(part)
        if isempty(pathText)
            pathText = char(part);
        else
            pathText = [pathText '.' char(part)]; %#ok<AGROW>
        end
    elseif isnumeric(part) && isscalar(part)
        pathText = sprintf('%s{%d}', pathText, part);
    end
end
end
