function dt_kx = make_dt_kx_from_sw(Nkx, sw, reverse)
%MAKE_DT_KX_FROM_SW Build EPSI kx timing offsets from sweep width.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

% reverse=false reproduces MIDAS default atr = -atr
ts = 1/sw;
atr = (0:Nkx-1) * ts / Nkx;
atr = atr - ts/2;
if ~reverse
    atr = -atr;
end
dt_kx = -atr;            % matches MIDAS dt_o = -atr
dt_kx = dt_kx(:);
end
