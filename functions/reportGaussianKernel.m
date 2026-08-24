function reportGaussianKernel(labelText, sigmaVec)
%REPORTGAUSSIANKERNEL Print Gaussian support and effective voxel count.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

info = gaussianKernelInfo(sigmaVec);
fprintf('%s: sigma = [%g %g %g], support = [%d %d %d], effective voxels = %.3f\n', ...
    labelText, sigmaVec(1), sigmaVec(2), sigmaVec(3), ...
    info.supportSize(1), info.supportSize(2), info.supportSize(3), info.effectiveVoxels);
end

function info = gaussianKernelInfo(sigmaVec)
filterSize = 2*ceil(2*sigmaVec) + 1;
x = -floor(filterSize(1)/2):floor(filterSize(1)/2);
y = -floor(filterSize(2)/2):floor(filterSize(2)/2);
z = -floor(filterSize(3)/2):floor(filterSize(3)/2);
[xx, yy, zz] = ndgrid(x, y, z);

kernel = exp(-((xx.^2)/(2*sigmaVec(1)^2) + ...
    (yy.^2)/(2*sigmaVec(2)^2) + ...
    (zz.^2)/(2*sigmaVec(3)^2)));
kernel = kernel / max(kernel(:));

info.supportSize = filterSize;
info.effectiveVoxels = sum(kernel(:));
end
