function dataOut = applyComplexSpatialFilter(dataIn, filterMode, gaussianSigma, medianKernel)
%APPLYCOMPLEXSPATIALFILTER Apply spatial filtering to each spectrum point.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

dataOut = zeros(size(dataIn), 'like', dataIn);
for kk = 1:size(dataIn,4)
    dataOut(:,:,:,kk) = applyOneVolumeFilter(dataIn(:,:,:,kk), filterMode, gaussianSigma, medianKernel);
end
end

function volOut = applyOneVolumeFilter(volIn, filterMode, gaussianSigma, medianKernel)
switch lower(filterMode)
    case 'gaussian'
        realPart = imgaussfilt3(real(volIn), gaussianSigma);
        imagPart = imgaussfilt3(imag(volIn), gaussianSigma);
    case {'median', 'moving-median', 'moving_median'}
        realPart = medfilt3(real(volIn), medianKernel, 'symmetric');
        imagPart = medfilt3(imag(volIn), medianKernel, 'symmetric');
    otherwise
        error('Unknown spatial filter mode "%s". Use "gaussian" or "median".', filterMode);
end
volOut = realPart + 1i*imagPart;
end
