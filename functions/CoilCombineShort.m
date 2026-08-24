function  [xtData_comb,xtData_comb_2,xtData_comb_3] = CoilCombineShort(xtData,xtData_2,xtData_3,opt,opt2)
%COILCOMBINESHORT Coil combination utility used by the SLOW-MRSI workflow.
% SLOW-MRSI package maintenance: Dr. Guodong Weng, University of Bern, 2026-07-04


weightSize = [size(xtData,1), size(xtData,2), size(xtData,3), size(xtData,5)];
Signal = zeros(weightSize);
Pha = zeros(weightSize, 'like', xtData);
Noise = zeros(weightSize);

for n1 = 1:size(xtData,1)
    %disp(n1)
    for n2 = 1:size(xtData,2)
        for n3 = 1:size(xtData,3)
            for n5 = 1:size(xtData,5)

                T(1,:)               = xtData(n1,n2,n3,1:1:end,n5);
                T2                   = T*exp(-1i*angle(T(1))*1);

                T3(1,:)              = xtData_2(n1,n2,n3,1:1:end,n5);
                T3                   = T3*exp(-1i*angle(T(1))*1);
                F3                   = fftshift(ifft(T3));
                
                switch opt2
                    case 'max'
                        Signal(n1,n2,n3,n5)  = max(abs(T2));
                    case 'linear'
                        X                    = [3,6];
                        Y                    = abs(T2(X));
                        f                    = polyfit(X-1,Y,1);
                        Signal(n1,n2,n3,n5)  = f(2);
                end

                Pha(n1,n2,n3,n5)     = -1i*angle(T(1));
                Noise(n1,n2,n3,n5)   = std(F3(1:100));


            end
        end
    end
end

Noise_sum = squeeze(sum(Noise,[1,2,3]));

xtData_corr = zeros(size(xtData));
xtData_corr_2 = zeros(size(xtData_2));
xtData_corr_3 = zeros(size(xtData_2));
for n1 = 1:size(xtData,1)
    %disp(n1)
    for n2 = 1:size(xtData,2)
        for n3 = 1:size(xtData,3)
            for n5 = 1:size(xtData,5)

                T(1,:)               = xtData(n1,n2,n3,1:1:end,n5);
                T_2(1,:)             = xtData_2(n1,n2,n3,1:1:end,n5);
                T_3(1,:)             = xtData_3(n1,n2,n3,1:1:end,n5);

                switch opt
                    case 'SNR'
                        clear Nor
                        Nor =  sum( (squeeze( Signal(n1,n2,n3,:)./Noise(n1,n2,n3,:) ) ).^2 );
                        xtData_corr(n1,n2,n3,:,n5)   = T.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)  /Noise(n1,n2,n3,n5)/sqrt(Nor);
                        xtData_corr_2(n1,n2,n3,:,n5) = T_2.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise(n1,n2,n3,n5)/sqrt(Nor);
                        xtData_corr_3(n1,n2,n3,:,n5) = T_3.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise(n1,n2,n3,n5)/sqrt(Nor);
                    case 'SN2'
                        clear Nor
                        Nor                          = sum( (squeeze( Signal(n1,n2,n3,:)./Noise(n1,n2,n3,:).^2 ) ).^2 );
                        xtData_corr(n1,n2,n3,:,n5)   = T.*  exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise(n1,n2,n3,n5).^2/sqrt(Nor);
                        xtData_corr_2(n1,n2,n3,:,n5) = T_2.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise(n1,n2,n3,n5).^2/sqrt(Nor);
                        xtData_corr_3(n1,n2,n3,:,n5) = T_3.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise(n1,n2,n3,n5).^2/sqrt(Nor);
                    case 'S'
                        xtData_corr(n1,n2,n3,:,n5) = T.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5);
                    case 'N'
                        xtData_corr(n1,n2,n3,:,n5) = T.*exp(Pha(n1,n2,n3,n5))/Noise(n1,n2,n3,n5);
                    case 'qSNR'
                        clear Nor
                        Nor                            = sum( (squeeze(Signal(n1,n2,n3,:))./Noise_sum).^2 );    %Noise(n1,n2,n3,:)).^2);
                        xtData_corr(n1,n2,n3,:,n5)     = T.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise_sum(n5)/sqrt(Nor);
                        xtData_corr_2(n1,n2,n3,:,n5) = T_2.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise_sum(n5)/sqrt(Nor);
                        xtData_corr_3(n1,n2,n3,:,n5) = T_3.*exp(Pha(n1,n2,n3,n5))*Signal(n1,n2,n3,n5)/Noise_sum(n5)/sqrt(Nor);
                end


            end
        end
    end
end

xtData_comb = sum(xtData_corr,5);
xtData_comb_2 = sum(xtData_corr_2,5);
xtData_comb_3 = sum(xtData_corr_3,5);
