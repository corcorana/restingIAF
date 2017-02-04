function [d0, d1, d2] = savGolDiff(x, poly, Fw, Fs, hWin)
% Savitzky-Golay Smoothing and Differentiation Filter for extracting
% estimates of the 0th (smoothed), 1st, & 2nd derivative of an input PSD.
%
% "Savitzky-Golay filters are optimal in the sense that they minimize the 
% least-squares error in fitting a polynomial to frames of noisy data."
%
% Depends on MATLAB Signal Processing Toolbox function to implement S-G 
% filter. 
%
% NB: sgolay preferred over sgolayfilt as furnishes coefficients that
% enable extraction of higher order derivatives (will come in handy later)
%
% Last modified AC, 31/01/2017.
%%
% Outputs:
%   d0 = smoothed PSD
%   d1 = 1st derivative of PSD
%   d2 = 2nd derivative of PSD
%
% Inputs:
%   x = spectral data to be filtered/differentiated
%   poly = polynomial order (for sgolay function ** must be < Fw **)
%   Fw = frame width (i.e. number of samples, for sgolay function ** must be odd **)
%   Fs = sampling rate
%   hWin = size (i.e. number of samples) of pwelch window
%%

[~, g] = sgolay(poly, Fw);      

% loop adapted from sgolay documentation
dt = Fs/hWin;
dx = zeros(length(x),3);
for p = 0:2                 % p determines order of estimated derivatives
    dx(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end

d0 = dx(:,1);        % smoothed signal post S-G diff filt
d1 = dx(:,2);        % 1st derivative
d2 = dx(:,3);        % 2nd derivative

end
