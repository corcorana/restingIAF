function [d0, d1, d2] = sgfDiff(x, Fw, poly, Fs, tlen)
% Savitzky-Golay Smoothing and Differentiation Filter for extracting
% estimates of the 0th (smoothed), 1st, & 2nd derivative function of an 
% input PSD.
%
% "Savitzky-Golay filters are optimal in the sense that they minimize the 
% least-squares error in fitting a polynomial to frames of noisy data."
%
% Depends on MATLAB Signal Processing Toolbox function to implement S-G 
% filter.
%%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
% Outputs:
%   d0 = smoothed PSD estimates
%   d1 = 1st derivative of d0
%   d2 = 2nd derivative of d0
%
% Inputs:
%   x = spectral data to be filtered/differentiated (vector)
%   poly = polynomial order (integer, must be < Fw)
%   Fw = frame width (i.e. number of samples, must be an odd integer)
%   Fs = sampling rate (integer)
%   tlen = taper length (i.e. number of samples of pwelch window, integer)

%% setup inputParser
p = inputParser;
p.addRequired('x',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('Fw',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive', 'odd'}));
p.addRequired('poly',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive', '<', Fw }));
p.addRequired('Fs',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
p.addRequired('tlen',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));  
p.parse(x, Fw, poly, Fs, tlen)
%%

[~, g] = sgolay(poly, Fw);      

dt = Fs/tlen;
dx = zeros(length(x),3);
for p = 0:2                 % p determines order of estimated derivatives
    dx(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end

d0 = dx(:,1);        % smoothed signal post S-G diff filt
d1 = dx(:,2);        % 1st derivative
d2 = dx(:,3);        % 2nd derivative

end
