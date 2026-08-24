function ktData_corrected = correctTsEPSIKx(ktData_shifted, reconInfo, progress)
%CORRECTTSEPSIKX Apply EPSI kx phase correction and remove readout oversampling.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

Nx = reconInfo.Nx;
Nkx = size(ktData_shifted.wat_e1,1);

if Nx == 48
    sw = 1/(500e-6);
    sw_fTop = 1/(240e-6);
elseif Nx == 40
    sw = 1/(800e-6);
    sw_fTop = 1/(640e-6);
elseif Nx == 60
    sw = 1/(540e-6);
    sw_fTop = 1/(360e-6);
elseif Nx == 15
    sw = 1/(260e-6);
    sw_fTop = 1/(120e-6);
else
    error('Unsupported Nx=%d. Add the EPSI sweep-width constants for this matrix size.', Nx);
end

dt_kx = make_dt_kx_from_sw(Nkx, sw_fTop, false);
sc1 = 1;

updateTsEPSIProgress(progress, 0.42, 'kx phase correction', 0/6, 'Water echo 1');
wat_e1_shp = kx_dt_phase_correct(ktData_shifted.wat_e1, sw, dt_kx, 'TimeDim', 4, 'KxDim', 1);
updateTsEPSIProgress(progress, 0.43, 'kx phase correction', 1/6, 'Water echo 2');
wat_e2_shp = kx_dt_phase_correct(ktData_shifted.wat_e2, sw, flip(dt_kx) + 1/sw*sc1*0, 'TimeDim', 4, 'KxDim', 1);
updateTsEPSIProgress(progress, 0.44, 'kx phase correction', 2/6, 'SLOW-ful echo 1');
met_s1_e1_shp = kx_dt_phase_correct(ktData_shifted.met_s1_e1, sw, dt_kx, 'TimeDim', 4, 'KxDim', 1);
updateTsEPSIProgress(progress, 0.45, 'kx phase correction', 3/6, 'SLOW-ful echo 2');
met_s1_e2_shp = kx_dt_phase_correct(ktData_shifted.met_s1_e2, sw, flip(dt_kx) + 1/sw*sc1, 'TimeDim', 4, 'KxDim', 1);
updateTsEPSIProgress(progress, 0.46, 'kx phase correction', 4/6, 'SLOW-par echo 1');
met_s2_e1_shp = kx_dt_phase_correct(ktData_shifted.met_s2_e1, sw, dt_kx, 'TimeDim', 4, 'KxDim', 1);
updateTsEPSIProgress(progress, 0.47, 'kx phase correction', 5/6, 'SLOW-par echo 2');
met_s2_e2_shp = kx_dt_phase_correct(ktData_shifted.met_s2_e2, sw, flip(dt_kx) + 1/sw*sc1, 'TimeDim', 4, 'KxDim', 1);

updateTsEPSIProgress(progress, 0.48, 'Remove oversampling', 0/6, 'Water echo 1');
Nx_trunc = size(wat_e1_shp,1)/2;
ktData_corrected.wat_e1 = remove_kx_oversampling_5D(wat_e1_shp, Nx_trunc);
updateTsEPSIProgress(progress, 0.50, 'Remove oversampling', 1/6, 'Water echo 2');
ktData_corrected.wat_e2 = remove_kx_oversampling_5D(wat_e2_shp, Nx_trunc);
updateTsEPSIProgress(progress, 0.52, 'Remove oversampling', 2/6, 'SLOW-ful echo 1');
ktData_corrected.met_s1_e1 = remove_kx_oversampling_5D(met_s1_e1_shp, Nx_trunc);
updateTsEPSIProgress(progress, 0.54, 'Remove oversampling', 3/6, 'SLOW-ful echo 2');
ktData_corrected.met_s1_e2 = remove_kx_oversampling_5D(met_s1_e2_shp, Nx_trunc);
updateTsEPSIProgress(progress, 0.56, 'Remove oversampling', 4/6, 'SLOW-par echo 1');
ktData_corrected.met_s2_e1 = remove_kx_oversampling_5D(met_s2_e1_shp, Nx_trunc);
updateTsEPSIProgress(progress, 0.58, 'Remove oversampling', 5/6, 'SLOW-par echo 2');
ktData_corrected.met_s2_e2 = remove_kx_oversampling_5D(met_s2_e2_shp, Nx_trunc);

updateTsEPSIProgress(progress, 0.60, 'Remove oversampling', 1, 'kx correction complete');
end
