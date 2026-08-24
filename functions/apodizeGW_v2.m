function [s_t_f2,s_t] = apodizeGW_v2(s,Dur,Hz,Hz_shf)

s = (squeeze(s));
if size(s,1) ~= 1
    s = transpose((s));
end
points = length(s);
X = 1:points;
t  =  linspace(0,Dur*1e-6*points,points);
s_t    = ifft(ifftshift((s))).*exp(-Hz.*t.^1).*exp(1i.*X/points*2*pi*(Hz_shf));
s_t_f    = fftshift(fft((s_t)));

if max(real(s_t_f)) == 0
    sc1 = 0;
else
    sc1 = 1;%max(real(s))/max(real(s_t_f));
end

s_t_f2 = s_t_f.*sc1/1;

end
