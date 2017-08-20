function [peakF, posZ1, posZ2, f1, f2, inf1, inf2, Q, Qf] = peakBounds(d0, d1, d2, f, w, minPow, minDiff, fres)
% Take derivatives from Savitzky-Golay curve-fitting and differentiation
% function sgfDiff, pump out estimates of alpha-band peak & bounds.
% Also calculates primary peak area Qf via integration between inflections.
%
% Depends on `findF1`, `findF2`, and `lessThan1` functions to locate 
% bounds of individual alpha band.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
% Outputs:
%   peakF = peak frequency estimate
%   posZ1 = freq of 1st positive zero-crossing (lower bound alpha interval)
%   posZ2 = freq of 2nd positive zero-crossing (upper bound alpha interval)
%   f1 = freq bin for posZ1
%   f2 = freq bin for posZ2
%   inf1 = inflection point, ascending edge
%   inf2 = inflection point, descending edge
%   Q = area under peak between inf1 & inf2
%   Qf = Q divided by bandwidth of Q
%
% Required inputs:
%   d0 = smoothed PSD estimate vector
%   d1 = 1st derivative vector
%   d2 = 2nd derivative vector
%   f = frequency bin vector
%   w = bounds of initial alpha window
%   minPow = vector of minimum power threshold values defining candidate peaks (regression fit of background spectral activity)
%   minDiff = minimum difference required to distinguish peak as dominant (proportion of primary peak height)
%   fres = frequency resolution (determine how many bins to search to establish shallow rolloff in d1)

%% setup inputParser
p = inputParser;
p.addRequired('d0',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('d1',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('d2',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('f',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'nonnegative', 'increasing'}));
p.addRequired('w',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'nonnegative', 'increasing', 'size', [1 2]}));    
p.addRequired('minPow',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('minDiff',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', '>=', 0, '<=', 1}));
p.addRequired('fres',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'positive'}));  
p.parse(d0, d1, d2, f, w, minPow, minDiff, fres)
%%

% evaluate derivative for zero-crossings

[~, lower_alpha] = min(abs(f-w(1)));      % set lower bound for alpha band
[~, upper_alpha] = min(abs(f-w(2)));      % set upper bound for alpha band


negZ = zeros(1,4);                              % initialise for zero-crossing count & frequency bin
cnt = 0;                                        % start counter at 0
for k = lower_alpha-1:upper_alpha+1             % step through frequency bins in alpha band (start/end at bound -/+ 1 to make sure don't miss switch)
    if sign(d1(k)) > sign(d1(k+1))              % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
       	[~, maxk] = max([d0(k), d0(k+1)]);      % ensure correct frequency bin is picked out (find larger of two values either side of crossing (in the smoothed signal))
       	if maxk == 1
        	maxim = k;
        elseif maxk == 2
            maxim = k+1;
        end
        cnt = cnt+1;                % advance counter by 1
      	negZ(cnt,1) = cnt;          % zero-crossing (i.e. peak) count
        negZ(cnt,2) = maxim;        % keep bin index for later
      	negZ(cnt,3) = f(maxim);     % zero-crossing frequency            
     	negZ(cnt,4) = d0(maxim);    % power estimate
    end
end
    
% sort out appropriate estimates for output
if negZ(1,1) == 0                   % if no zero-crossing detected --> report NaNs
    peakF = NaN;
    subBin = NaN;
elseif size(negZ, 1) == 1           % if singular crossing...
    if log10(negZ(1, 4)) > minPow(negZ(1,2))      % ...and peak power is > minimum threshold --> report frequency
        peakBin = negZ(1, 2);
        peakF = negZ(1, 3);
    else
        peakF = NaN;                % ...otherwise, report NaNs
        subBin = NaN;        
    end
else negZ = sortrows(negZ, -4);     % if >1 crossing, re-sort from largest to smallest peak...
    if log10(negZ(1, 4)) > minPow(negZ(1,2))        % ...if highest peak exceeds min threshold...
        if negZ(1, 4)*(1-minDiff) > negZ(2, 4)      % ...report frequency of this peak.
        	peakBin = negZ(1, 2);
            peakF = negZ(1, 3); 
        else                        % ...if not...
            peakF = NaN;
            subBin = negZ(1, 2);                    % ... index as a subpeak for starting alpha bound search.
        end
    else
        peakF = NaN;                % ...otherwise, report NaNs
        subBin = NaN;
    end
end


%% search for positive (upward going) zero-crossings (minima / valleys) either side of peak/subpeak(s)
slen = round(1/fres);               % define number of bins included in shollow slope search (approximate span = 1 Hz)

if isnan(peakF) && isnan(subBin);       % if no evidence of peak activity, no parameter estimation indicated
    
    posZ1 = NaN;
    posZ2 = NaN;
    f1 = NaN;
    f2 = NaN;
    inf1 = NaN;
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

elseif isnan(peakF)     % deal with spectra lacking a clear primary peak (similar strategy to peak; take highest subpeak as start point, look for minima)
    
    [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, subBin); 
    [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, subBin); 
        
    % inflections / Q values not calculated as these spectra won't be included in averaged channel peak analyses 
    inf1 = NaN;     
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

else            % now for the primary peak spectra
    
    [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, peakBin);  
    [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, peakBin); 
    
    % define boundaries by inflection points (requires 2nd derivative of smoothed signal)
    inf1 = zeros(1,2);                  % initialise for zero-crossing count & frequency
    cnt = 0;                            % start counter at 0
    for k = 1:peakBin-1                 % step through frequency bins prior peak
        if sign(d2(k)) > sign(d2(k+1))                  % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
            [~, mink] = min(abs([d2(k), d2(k+1)]));     % ensure correct frequency bin is picked out (find smaller of two values either side of crossing)
            if mink == 1
                min1 = k;
            else
                min1 = k+1;
            end
            cnt = cnt+1;                % advance counter by 1
            inf1(cnt,1) = cnt;          % zero-crossing count
            inf1(cnt,2) = f(min1);      % zero-crossing frequency
        end
    end

    % sort out appropriate estimates for output
    if size(inf1, 1) == 1               % if singular crossing --> report frequency
        inf1 = inf1(1, 2);
    else
        inf1 = sortrows(inf1, -2);      % sort by frequency values (descending)...
        inf1 = inf1(1, 2);              % take highest frequency (bin nearest to peak)
    end
    
    for k = peakBin+1:length(d2)-1                      % step through frequency bins post peak
        if sign(d2(k)) < sign(d2(k+1))                  % look for upward zero-crossing
            [~, mink] = min(abs([d2(k), d2(k+1)]));     % ensure frequency bin nearest zero-crossing point picked out (find smaller of two values either side of crossing)
                if mink == 1
                    min2 = k;
                else
                    min2 = k+1;
                end
                inf2 = f(min2);         % zero-crossing frequency
                break                   % break loop (only need to record first crossing)
                
        end

    end

    % estimate approx. area under curve between inflection points either
    % side of peak, scale by inflection band width 
    Q = trapz(f(min1:min2), d0(min1:min2));
    Qf = Q / (min2-min1);

end 

end

