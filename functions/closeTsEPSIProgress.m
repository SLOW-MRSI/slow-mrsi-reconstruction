function closeTsEPSIProgress(progress)
%CLOSETSEPSIPROGRESS Close the ts-EPSI progress window if it is open.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

if isstruct(progress) && isfield(progress, 'enabled') && progress.enabled ...
        && isfield(progress, 'h') && ~isempty(progress.h) && isgraphics(progress.h)
    delete(progress.h);
end
end
