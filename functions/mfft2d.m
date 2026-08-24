function k   = mfft2d(x,norm_on,shift_on)
% perform 2d fft with correct normalization and shift
% Chao Ma
% SLOW-MRSI package maintenance: Dr. Guodong Weng, University of Bern, 2026-07-04
% 

temp         = size(x);
Ny           = temp(1);
Nx           = temp(2);

% by default, we do normalization
if nargin<2
	norm_on  = 1;
end
if nargin<3
	shift_on = 1;
end

if shift_on
	% fft with correct shift
	k        = fftshift(fft(ifftshift(...
	    	   fftshift(fft(ifftshift(...
	    	   x,...
	    	   1),[],1),1),...
	    	   2),[],2),2);
else
	k        = fft(fft(x,[],1),[],2);
end

% fft with normalization
if norm_on   == 1
	k        = k./sqrt(Ny)./sqrt(Nx);
end
