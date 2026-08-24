function data_out = remove_kx_oversampling_5D(data_in, Nx_trunc)
%REMOVE_KX_OVERSAMPLING_5D Remove readout oversampling from 5-D EPSI data.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

    % Transform from kx to x
    data_x = ifftshift(ifft(fftshift(data_in, 1), [], 1), 1);

    % Center crop along first dimension
    Nx = size(data_x, 1);
    i1 = Nx/2 + 1 - Nx_trunc/2;
    i2 = Nx/2 + Nx_trunc/2;
    data_x = data_x(i1:i2, :, :, :, :);

    % Transform back to kx
    data_out = fftshift(fft(ifftshift(data_x, 1), [], 1), 1);
end
