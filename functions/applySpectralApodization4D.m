function dataOut = applySpectralApodization4D(dataIn, dwellTimeUs, linewidthHz)
%APPLYSPECTRALAPODIZATION4D Apply exponential linewidth apodization along dim 4.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-08

nPts = size(dataIn, 4);
t = linspace(0, dwellTimeUs/1000/1000*nPts, nPts);
weight = reshape(exp(-linewidthHz.*t), [1, 1, 1, nPts]);
if isa(dataIn, 'single')
    weight = single(weight);
end

timeData = ifft(ifftshift(dataIn, 4), [], 4) .* weight;
dataOut = fftshift(fft(timeData, [], 4), 4);
end
