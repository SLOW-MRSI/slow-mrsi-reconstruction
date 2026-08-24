function Dcorr = kx_dt_phase_correct(D, sw, dt_kx, varargin)
%KX_DT_PHASE_CORRECT  First-order phase correction due to dt(kx) sampling time offsets.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04
%
% Implements:  Dcorr(kx,...,t,...) = IFFT_t( FFT_t(D) .* exp(1i * omega * dt_kx(kx)) )
%
% INPUTS
%   D      : complex array, e.g. [Nkx Ny Nz Nt Nc] (or with extra dims like Nave)
%   sw     : sweep width along t (Hz). In MIDAS this is block.sw.
%   dt_kx  : [Nkx x 1] time offsets in seconds for each kx sample
%
% NAME-VALUE
%   'TimeDim' : which dim is t (default 4, matches your code: fft(...,[],4))
%   'KxDim'   : which dim is kx (default 1)
%   'PadToPow2' : true/false, zero-pad time to next pow2 for FFT (default false)
%
% OUTPUT
%   Dcorr : same size as D, phase-corrected
%
% Notes
%   - dt_kx sign matters. If correction goes the wrong way, flip its sign.
%   - This is the "phase correction in frequency domain then IFFT" you described.

p = inputParser;
p.addParameter('TimeDim', 4, @(x)isscalar(x)&&x>=1);
p.addParameter('KxDim', 1, @(x)isscalar(x)&&x>=1);
p.addParameter('PadToPow2', false, @(x)islogical(x)||isscalar(x));
p.parse(varargin{:});
tDim = p.Results.TimeDim;
kxDim = p.Results.KxDim;
doPad = logical(p.Results.PadToPow2);

sz = size(D);
nd = ndims(D);

if tDim > nd || kxDim > nd
    error('TimeDim or KxDim exceeds ndims(D).');
end

Nkx = sz(kxDim);
Nt  = sz(tDim);

dt_kx = dt_kx(:);
if numel(dt_kx) ~= Nkx
    error('dt_kx length (%d) must match Nkx=size(D,KxDim)=%d.', numel(dt_kx), Nkx);
end

% --- build frequency axis (Hz) consistent with fftshift(fft(...)) like your code ---
% MIDAS uses: ww = ((n/nt3)-0.5)*sw*pi  (rad/s)
% Here: f = ((n/N)-0.5)*sw  (Hz), omega = 2*pi*f (rad/s)
if doPad
    Nfft = 2^nextpow2(Nt);
else
    Nfft = Nt;
end

f = ((0:Nfft-1)/Nfft - 0.5) * sw/2;      % Hz, centered (fftshift convention)
omega = 2*pi*f(:);                     % rad/s, column [Nfft x 1]

% reshape omega so it broadcasts along the time dimension
omegaShape = ones(1, nd);
omegaShape(tDim) = Nfft;
omegaArr = reshape(omega, omegaShape);

% reshape dt_kx so it broadcasts along the kx dimension
dtShape = ones(1, nd);
dtShape(kxDim) = Nkx;
dtArr = reshape(dt_kx, dtShape);

% phase term exp(i * omega * dt(kx))
phase = exp(-1i * omegaArr .* dtArr);

%%
% figure(111)
% clear temp
% temp(1,:) = phase(40,1,1,:);
% 
% plot(real(temp))
% ylim([-1,1])
% f(1)*dt_kx(40)
%%

% --- FFT along time ---
if doPad
    % pad along time dim
    padSz = sz;
    padSz(tDim) = Nfft;
    Dpad = complex(zeros(padSz, 'like', D));
    idx = repmat({':'}, 1, nd);
    idx{tDim} = 1:Nt;
    Dpad(idx{:}) = D;
    Sf = fftshift(fft(Dpad, [], tDim), tDim);
    Sf = Sf .* phase;
    Difft = ifft(ifftshift(Sf, tDim), [], tDim);
    % crop back
    Dcorr = Difft(idx{:});
else
    Sf = fftshift(fft(D, [], tDim), tDim);
    Sf = Sf .* phase;
    Dcorr = ifft(ifftshift(Sf, tDim), [], tDim);
end
end
