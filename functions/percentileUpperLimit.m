function upperLimit = percentileUpperLimit(mapData, percentileValue)
%PERCENTILEUPPERLIMIT Robust positive percentile for map display scaling.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

values = double(mapData(:));
values = values(isfinite(values) & values > 0);

if isempty(values)
    upperLimit = 1;
    return
end

values = sort(values);
idx = max(1, min(numel(values), ceil(percentileValue/100*numel(values))));
upperLimit = values(idx);

if upperLimit <= 0
    upperLimit = max(values);
end
if upperLimit <= 0
    upperLimit = 1;
end
end
