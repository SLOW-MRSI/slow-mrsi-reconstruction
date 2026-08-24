function plotTsEPSIPreviewSpectra(xfDataGau, Pha_sel, Pha_s, X_ppm_met, cfg)
%PLOTTSEPSIPREVIEWSPECTRA Plot a small user-selected voxel grid.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

Pha = 0.0;
X = cfg.previewX;
Y = cfg.previewY;
Z = cfg.previewZ;
EPSIDurUs = 1e6 / cfg.spectralBandwidthHz;
linewidthFullHz = 10;
linewidthEditedHz = 10;
linewidthDifferenceHz = 15;

figure(1)
dn = 1;
for n2 = 1:length(Y)
    for n1 = 1:length(X)
        subplot(length(Y),length(X),dn)
        dn = dn + 1;

        pha = exp(1i*pi*Pha_s(Pha_sel(Y(n2),X(n1),Z)));
        Y1 = squeeze(xfDataGau.ful(Y(n2),X(n1),Z,:)) * pha;
        Y2 = squeeze(xfDataGau.par(Y(n2),X(n1),Z,:)) * pha;
        Y1b = apodizeGW_v2(Y1, EPSIDurUs, linewidthFullHz, 0);
        Y2b = apodizeGW_v2(Y2, EPSIDurUs, linewidthEditedHz, 0);

        plot(X_ppm_met, 1e4*real(Y1b.*exp(1i*pi*Pha)), X_ppm_met, 1e4*real(Y2b.*exp(1i*pi*Pha)));
        hold off
        xlim([min(X_ppm_met),max(X_ppm_met)])
        set(gca, 'XDir','reverse')
    end
end

figure(2)
dn = 1;
for n2 = 1:length(Y)
    for n1 = 1:length(X)
        subplot(length(Y),length(X),dn)
        dn = dn + 1;

        pha = exp(1i*pi*Pha_s(Pha_sel(Y(n2),X(n1),Z)));
        Y3 = squeeze(xfDataGau.dif(Y(n2),X(n1),Z,:)) * pha;
        Y3b = apodizeGW_v2(Y3, EPSIDurUs, linewidthDifferenceHz, 0);

        plot(X_ppm_met, 1e4*real(Y3b.*exp(1i*pi*Pha)) ,'m');
        hold off
        xlim([min(X_ppm_met),max(X_ppm_met)])
        set(gca, 'XDir','reverse')
    end
end
end
