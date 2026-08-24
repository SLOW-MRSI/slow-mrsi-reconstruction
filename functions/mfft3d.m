function k   = mfft3d(x,norm_on,shift_on)
% perform 3d fft with correct normalization and shift
% Chao Ma
% SLOW-MRSI package maintenance: Dr. Guodong Weng, University of Bern, 2026-07-04
% 

Ny           = size(x,1);
Nx           = size(x,2);
Nz           = size(x,3);

% by default, we do normalization
if nargin<2
	norm_on  = 1;
end
if nargin<3
	shift_on = 1;
end

if Nz > 1
	if shift_on
		% fft with correct shift
		k    = fftshift(fft(ifftshift(...
	    	   fftshift(fft(ifftshift(...
	    	   fftshift(fft(ifftshift(...
	    	   x,...
	    	   1),[],1),1),...
	    	   2),[],2),2),...
	    	   3),[],3),3);
	else
        k    = fft(fft(fft(x,[],1),[],2),[],3);
	end

	% fft with normalization
	if norm_on   == 1
		k    = k./sqrt(Ny)./sqrt(Nx)./sqrt(Nz);
	end
else
	k        = mfft2d(x,norm_on,shift_on);
end
