function plotTsEPSIMaps(xfData, Nz)
%PLOTTSEPSIMAPS Plot summed full, edited, and water maps.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

figure(3)
clf
dn = 0;
for nz = 1:Nz
    dn = dn + 1;

    subplot(3,Nz,0+dn)
    imagesc(abs(sum((xfData.ful(:,:,nz,:)), 4)))
    title('SLOW-ful')

    subplot(3,Nz,Nz+dn)
    imagesc(abs(sum(xfData.par(:,:,nz,:), 4)))
    title('SLOW-par')

    subplot(3,Nz,2*Nz+dn)
    imagesc(abs(sum(xfData.wat(:,:,nz,:), 4)))
    title('SLOW-wat')
end
end
