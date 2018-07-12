function [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, bin)
% Searches 1st derivative for evidence of local minima or near horizontal 
% function prior to alpha peak. This location will be taken as the lower 
% bound of the individual alpha band used to calculate CoG (f1).
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Outputs:
%   f1 = frequency bin indexing f1
%   posZ1 = frequency of f1
%
%% Required inputs:
%   f = frequency bin vector
%   d0 = smoothed PSD estimate vector
%   d1 = 1st derivative vector
%   negZ = vector / matrix of negative zero-crossings detected within alpha window
%   minPow = vector of minimum power threshold values defining candidate peaks
%   slen = number of bins examined for evidence of shallow slope (~1 Hz interval)
%   bin = frequency bin indexing peak / subpeak

%% setup inputParser
p = inputParser;
p.addRequired('f',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'nonnegative', 'increasing'}));
p.addRequired('d0',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('d1',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('negZ',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'2d', 'nonempty'}));
p.addRequired('minPow',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('slen',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive', '>', 1}));
p.addRequired('bin',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
p.parse(f, d0, d1, negZ, minPow, slen, bin)
%%
            
posZ1 = zeros(1,4);
    
% contingency for multiple peaks - aim to identify left-most peak in range for upper bound of k in next loop (avoid falling into local minimum)
if size(negZ, 1) >1                 
 	negZ = sortrows(negZ, 3);       % sort by frequency (ascending) 
  	for z = 1:size(negZ, 1)                 
        if log10(negZ(z, 4)) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(bin)) % relax power constraint, as may result in excessively narrow alpha window in noisy spectra with shallow peakF (i.e. precisely where we want CoG to take breadth into account)
         	leftPeak = negZ(z, 2);
          	break                   % break off search when conditions satisfied
      	else leftPeak = bin;        % if fail to satisfy search conditions, default to bin (sub)peak
        end
    end
else leftPeak = bin;                % if no other peaks were identified, take bin (sub)peak as boundary   
end
    
cnt = 0;                        % start counter at 0
for k = 2:leftPeak-1            % step through frequency bins up to left-most peak in search window
	if sign(d1(k)) < sign(d1(k+1))        % look for switch from negative to positive derivative values (i.e. upward/positive zero-crossing)
    	[~, mink] = min(abs([d0(k-1), d0(k), d0(k+1)]));    % search around crossing for local minimum in d0 (indexing 1st derivative sometimes results in small errors)
        if mink == 1
         	minim = k-1;
       	elseif mink == 2
          	minim = k;
        else
           	minim = k+1;
        end
            
     	cnt = cnt+1;                % advance counter by 1
      	posZ1(cnt,1) = cnt;         % zero-crossing count
       	posZ1(cnt,2) = minim;       % zero-crossing frequency bin
      	posZ1(cnt,3) = f(minim);    % zero-crossing frequency

   	% look for consistent low d1 values for signs of shallow slope (levelling off)
  	elseif abs(d1(k)) < 1 && lessThan1(d1(k+1:k+slen))
     	minim = k;
     	cnt = cnt+1;                % advance counter by 1
     	posZ1(cnt,1) = cnt;         % zero-crossing count
       	posZ1(cnt,2) = minim;       % zero-crossing frequency bin
       	posZ1(cnt,3) = f(minim);    % zero-crossing frequency
    end
    
end

% sort out appropriate estimates for output
if size(posZ1, 1) == 1              % if singular crossing --> report frequency
	f1 = posZ1(1, 2);
	posZ1 = posZ1(1, 3);
else                                % else sort by frequency values (descending), take highest frequency (bin nearest to peak)
	posZ1 = sortrows(posZ1, -3);    
	f1 = posZ1(1, 2);
	posZ1 = posZ1(1, 3);
end
    
end
