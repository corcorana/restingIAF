function [peakF, posZ1, posZ2, f1, f2, inf1, inf2, Q, Qf] = alphaParams(d0, d1, d2, f, w, minPow, minDiff, t)
% Take derivatives from Savitzky-Golay curve-fitting and differentiation
% function savGolDiff, pump out estimates of alpha band peak & bounds.
% Also calculate primary peak area Qf via integration between inflections.
%
% Last modified AC 27/02/2017 - correct sign error on 2nd deriv inflection
% search (should be upward not downward sign change)
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
% Inputs:
%   d0 = smoothed derivative vector
%   d1 = 1st derivative vector
%   d2 = 2nd derivative vector
%   f = freq vector
%   w = bounds of initial alpha window
%   minPow = minimum threshold vector (regression fit of background
%   spectral activity)
%   minDiff = minimum difference required to distinguish peak as be primary (proportion of primary peak height)
%   t = partially depricated, currently only defines limits of minima search window
%%
% evaluate derivative for zero-crossings

[~, lower_alpha] = min(abs(f-w(1)));      % set lower bound for alpha band
[~, upper_alpha] = min(abs(f-w(2)));      % set upper bound for alpha band


    negZ = zeros(1,4);                     % initialise for zero-crossing count & frequency bin (?? is this sort of thing helpful in functions)
    cnt = 0;                            % start counter at 0
    % NB: following doesn't look very efficient, but good for coping with
    % varying numbers of zero-crossings (and ensuring bin containing value closest to crossing is indexed)
    for k = lower_alpha:upper_alpha-1     % step through frequency bins in alpha band
        if sign(d1(k)) > sign(d1(k+1))            % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
            %if range(sign(d1sgf(k+1:k+t)))==0           % must remain
            %downgoing for next t frequency bins (exclude noisy
            %fluctuations about zero) - ? unnecessary due to smoothing
                [~, maxk] = max([d0(k), d0(k+1)]);    % ensure correct frequency bin is picked out (find larger of two values either side of crossing (in the smoothed signal))
                if maxk == 1
                    maxim = k;
                elseif maxk == 2
                    maxim = k+1;
                end
                cnt = cnt+1;                % advance counter by 1
                negZ(cnt,1) = cnt;          % zero-crossing (i.e. peak) count
                negZ(cnt,2) = maxim;        % keep bin index for later
                negZ(cnt,3) = f(maxim);     % zero-crossing frequency            
                negZ(cnt,4) = d0(maxim);     % power estimate (original signal - as smoothing will diminish)
            %end
        end
    end
    
% sort out appropriate estimates for output
peakBin = 0;
peakF = 0;
subBin = 0;
if negZ(1,1) == 0               % if no zero-crossing detected --> report NaNs
    peakF = NaN;
    subBin = NaN;
elseif size(negZ, 1) == 1       % if singular crossing...
    if negZ(1, 4) > minPow(negZ(1,2))      %       ...and peak power is > minimum threshold --> report frequency
        peakBin = negZ(1, 2);
        peakF = negZ(1, 3);
    else
        peakF = NaN;           %       ...otherwise, report NaNs
        subBin = NaN;        
    end
else negZ = sortrows(negZ, -4);     % if >1 crossing, re-sort from largest to smallest peak...
    if negZ(1, 4) > minPow(negZ(1,2))          %       ...if highest peak exceeds min threshold...
        if negZ(1, 4)*(1-minDiff) > negZ(2, 4)      %...report frequency of this peak.
        	peakBin = negZ(1, 2);
            peakF = negZ(1, 3); 
        else                        %       ...if not...
            peakF = NaN;
            subBin = negZ(1, 2);                    %... index as a subpeak for starting alpha bound search.
        end
    else
        peakF = NaN;               %   ...otherwise, report NaNs
        subBin = NaN;
    end
end


%% search for positive (upward going) zero-crossings (minima / valleys) either side of peak/subpeak(s)
% altered approach 11/1/17 - don't accept subpeaks < minPow, as likely just
% noise (according to assumptions of mPow). can't compare with PAF either. CoG
% primarily of interest in cases of multiple substantial peaks, or
% broadband alpha (but could find f1/f2 in better chans, calculate CoG all chans)
% 18/1 - but can still set bounds of search between a small alpha band subpeak and end
% of spectral range

posZ1 = zeros(1,4);                     % initialise for zero-crossing count & frequency bin
posZ2 = zeros(1,4);                     % initialise for zero-crossing count & frequency bin

if isnan(peakF) && isnan(subBin);        % no point doing anything if the PSD is just noise, give up and move on
    
    posZ1 = NaN;
    posZ2 = NaN;
    f1 = NaN;
    f2 = NaN;
    inf1 = NaN;     % inflections / Q values not calculated as these spectra aren't included in averaged channel peak analyses 
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

elseif isnan(peakF)         % first, deal with spectra lacking a clear primary peak (similar strategy to peak; take highest subpeak as start point, look for minima)
    % contingency for multiple peaks - aim to identify left-most peak in range for upper bound of k in next loop (avoid falling into local minima)

    inf1 = NaN;     % inflections / Q values not calculated as these spectra aren't included in averaged channel peak analyses 
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;
    
    if size(negZ, 1) >1                 
        negZ = sortrows(negZ, 3);           % sort by frequency (ascending) 
        for z = 1:size(negZ, 1)                 
            if negZ(z, 4) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(subBin)) %&& range([d1(negZ(z, 2)-1), d1(negZ(z, 2)+1)]) > 1      % search for lowerst frequency bin above minimal power threshold (but not just part of shallow rolloff NB need to assess around peak, as peak d1 will be close to 0)
                leftPeak = negZ(z, 2);      % NB: minPow threshold still relevant as may have failed to identify primary peak due to presence of >1 distinct peaks 
                break                       % break off search when conditions satisfied
            else leftPeak = subBin;        % if fail to satisfy search conditions, default to subPeak
            end
            
        end
    else leftPeak = subBin;
    end
        
  	cnt = 0;                    % start counter at 0
    for k = t:leftPeak-1     % step through frequency bins up to left-most peak
      	if sign(d1(k)) < sign(d1(k+1)) %&& range(abs(d1(k-3:k-1))) < 1         % look for switch from negative to positive derivative values (i.e. upward zero-crossing) [ && stay low (?? necessary with leftPeak)]
             %   if range(sign(d1(k-t:k-1)))==0           % must remain downgoing for preceeding t frequency bins (exclude noisy fluctuations about zero)

        	[~, mink] = min(abs([d0(k-1), d0(k), d0(k+1)]));    % search around crossing for local minimum in d0 - indexing 1st derivative sometimes results in small errors
            if mink == 1
             	minim1 = k-1;
          	elseif mink == 2
              	minim1 = k;
         	elseif mink == 3
              	minim1 = k+1;
            end
            
            cnt = cnt+1;            % advance counter by 1
          	posZ1(cnt,1) = cnt;         % zero-crossing count
         	posZ1(cnt,2) = minim1;      % bin for CoG
           	posZ1(cnt,3) = f(minim1);        % zero-crossing frequency
          	posZ1(cnt,4) = d0(minim1);        % power estimate
             %   end
            % try diminished d1 for shallow slope
        elseif abs(d1(k)) < 1 && d1(k) < d1(k+1) && d1(k+1) < d1(k+2) && d1(k+2) < d1(k+3) && d1(k+3) < d1(k+4)  % for shallow d1, values should be getting more positive as approach leading/ascending edge
            %abs(d1(k)) < 1 & abs(d1(k:k+t)) < 1;   %range(d1(k-4:k)) < 1 % need 2nd condition to guard against exclusion of shoulder portion of curve (that fall short of registering as secondary peaks)
            minim1 = k;
            cnt = cnt+1;            % advance counter by 1
            posZ1(cnt,1) = cnt;         % zero-crossing count
            posZ1(cnt,2) = minim1;
            posZ1(cnt,3) = f(minim1);        % zero-crossing frequency
            posZ1(cnt,4) = d0(minim1);        % power estimate

        end
    end
        
        % sort out appropriate estimates for output
  	if size(posZ1, 1) == 1      % if singular crossing --> report frequency
        f1 = posZ1(1, 2);
        posZ1 = posZ1(1, 3);
  	elseif size(posZ1, 1) >1        % if > 1 crossing detected...
      	posZ1 = sortrows(posZ1, -3);        % sort by frequency values (descending)...
      	f1 = posZ1(1, 2);
        posZ1 = posZ1(1, 3);                  % take highest frequency (bin nearest to peak) %% PROBLEM - might choose adjacent bin if shallow peak ??impose conditional, must be intervening d1 acceleration, otherwise just search for minimum/zcross
    end

    
    %% look for end of alpha peak
        
    % contingency for multiple peaks - try to identify left-most peak in range for upper bound of k in next loop (avoid falling into local minima)
    if size(negZ, 1) >1                 
        negZ = sortrows(negZ, -3);          % sort by frequency (descending)
        for z = 1:size(negZ, 1)
            if negZ(z, 4) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(subBin))  % && range([d1(negZ(z, 2)-1), d1(negZ(z, 2)+1)]) > 1          % look for first freq bin with power > min threshold and adjacent d1 > 1
                rightPeak = negZ(z, 2);     
                break                       % break search when conditions satisfied
            else rightPeak = subBin;       % if fail to satisfy search conditions, default to subBin
            end
        end
    else rightPeak = subBin;                % if no other peaks were identified, take subBin as boundary
    end
    
        
	cnt = 0;                            % start counter at 0
  	for k = rightPeak+1:length(d1)-t     % step through frequency bins following right-most peak (cut of end of range to allow for conditional search ahead of d1 < 1)
      	if sign(d1(k)) < sign(d1(k+1))            % look for switch from negative to positive derivative values (i.e. upward zero-crossing)
             %   if range(sign(d1(k+1:k+t)))==0           % must remain downgoing for next t frequency bins (exclude noisy fluctuations about zero)

         	[~, mink] = min(abs([d0(k), d0(k+1)]));    % ensure correct frequency bin is picked out (find smaller of two values either side of crossing (in original pxx signal))
            if mink == 1
             	minim2 = k;
            elseif mink == 2
            	minim2 = k+1;
            end
            
            cnt = cnt+1;            % advance counter by 1
          	posZ2(cnt,1) = cnt;         % zero-crossing count
          	posZ2(cnt,2) = minim2;        % bin for CoG
           	posZ2(cnt,3) = f(minim2);        % zero-crossing frequency
           	posZ2(cnt,4) = d0(minim2);        % power estimate (from original signal)
             %   end
            % try to make use of shallow slope in d1 as (early) indicator of peak boundary
        elseif abs(d1(k)) < 1 && d1(k) < d1(k+1) && d1(k+1) < d1(k+2) && d1(k+2) < d1(k+3) && d1(k+3) < d1(k+4)     % descending limb, 1st d values are -ve but approaching zero
            %abs(d1(k)) < 1 & abs(d1(k:k+t)) < 1        %range(d1(k:k+4)) < 1
               %% ??need a check that there isn't a more prominant drop off (or indeed, secondary peak) later (i.e confirm that decrease is maintained to some extent) - have cake (cut off shallow slope sans zcross) and eat (don't cut off too soon when there's a peak still to come) 
                minim2 = k;
                cnt = cnt+1;            % advance counter by 1
                posZ2(cnt,1) = cnt;         % zero-crossing count
                posZ2(cnt,2) = minim2;
                posZ2(cnt,3) = f(minim2);        % zero-crossing frequency
                posZ2(cnt,4) = d0(minim2);        % power estimate (from original signal)
        end
    end
    
    f2 = posZ2(1, 2);
    posZ2 = posZ2(1, 3);            % just take first estimate for output
    
else            % now for the primary peak spectra
    %% look for beginning of alpha peak
    
    % contingency for multiple peaks - aim to identify left-most peak in range for upper bound of k in next loop (avoid falling into local minima)
    if size(negZ, 1) >1                 
        negZ = sortrows(negZ, 3);           % sort by frequency (ascending) 
        for z = 1:size(negZ, 1)                 
            if negZ(z, 4) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(peakBin)) % relax power constraint, as may result in excessively narrow alpha window in noisy spectra with shallow peakF (i.e. precisely where we want CoG to take breadth into account)
                leftPeak = negZ(z, 2);
                break                       % break off search when conditions satisfied
            else leftPeak = peakBin;        % if fail to satisfy search conditions, default to peakBin
            end
        end
    else leftPeak = peakBin;                % if no other peaks were identified, take peakF as boundary   
    end
    
    cnt = 0;                    % start counter at 0
    for k = t:leftPeak-1     % step through frequency bins up to left-most peak
        if sign(d1(k)) < sign(d1(k+1)) %&& range(abs(d1(k-3:k-1))) < 1         % look for switch from negative to positive derivative values (i.e. upward zero-crossing) [ && stay low (?? necessary with leftPeak)]
           % if range(sign(d1(k-t:k-1)))==0           % must remain downgoing for next t frequency bins (exclude noisy fluctuations about zero)

         	[~, mink] = min(abs([d0(k-1), d0(k), d0(k+1)]));    % search around crossing for local minimum in d0 - indexing 1st derivative sometimes results in small errors
           	if mink == 1
            	minim1 = k-1;
         	elseif mink == 2
               	minim1 = k;
           	elseif mink == 3
              	minim1 = k+1;
            end
            
            cnt = cnt+1;            % advance counter by 1
           	posZ1(cnt,1) = cnt;         % zero-crossing count
           	posZ1(cnt,2) = minim1;        % bin for CoG
           	posZ1(cnt,3) = f(minim1);        % zero-crossing frequency
           	posZ1(cnt,4) = d0(minim1);        % power estimate
          %  end
        % try diminished d1 for shallow slope
        elseif abs(d1(k)) < 1 && d1(k) < d1(k+1) && d1(k+1) < d1(k+2) && d1(k+2) < d1(k+3) && d1(k+3) < d1(k+4)
            %abs(d1(k)) < 1 & abs(d1(k:k+t)) < 1   % range(d1(k-4:k)) <1 % need 2nd condition to guard against exclusion of shoulder portion of curve (that fall short of registering as secondary peaks)
            minim1 = k;
            cnt = cnt+1;            % advance counter by 1
            posZ1(cnt,1) = cnt;         % zero-crossing count
            posZ1(cnt,2) = minim1;
            posZ1(cnt,3) = f(minim1);        % zero-crossing frequency
            posZ1(cnt,4) = d0(minim1);        % power estimate
            
        end
    end


    % sort out appropriate estimates for output
    if size(posZ1, 1) == 1      % if singular crossing --> report frequency
        f1 = posZ1(1, 2);
        posZ1 = posZ1(1, 3);
    elseif size(posZ1, 1) >1        % if > 1 crossing detected...
        posZ1 = sortrows(posZ1, -3);        % sort by frequency values (descending)...
        f1 = posZ1(1, 2);
        posZ1 = posZ1(1, 3);                  % take highest frequency (bin nearest to peak) %% PROBLEM - might choose adjacent bin if shallow peak ??impose conditional, must be intervening d1 acceleration, otherwise just search for minimum/zcross
    end
    
    
    %% look for end of alpha peak
        
    % contingency for multiple peaks - try to identify left-most peak in range for upper bound of k in next loop (avoid falling into local minima)
    if size(negZ, 1) >1                 
        negZ = sortrows(negZ, -3);          % sort by frequency (descending)
        for z = 1:size(negZ, 1)
            if negZ(z, 4) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(peakBin)) % && range([d1(negZ(z, 2)-1), d1(negZ(z, 2)+1)]) > 1          % look for first freq bin with power > min threshold and adjacent d1 > 1
                rightPeak = negZ(z, 2);     
                break                       % break search when conditions satisfied
            else rightPeak = peakBin;       % if fail to satisfy search conditions, default to peakBin
            end
         end
    else rightPeak = peakBin;                % if no other peaks were identified, take peakF as boundary
    end
    
    cnt = 0;                            % start counter at 0
    for k = rightPeak+1:length(d1)-t     % step through frequency bins following right-most peak
        if sign(d1(k)) < sign(d1(k+1))            % look for switch from negative to positive derivative values (i.e. upward zero-crossing)
          %  if range(sign(d1(k+1:k+t)))==0           % must remain downgoing for next t frequency bins (exclude noisy fluctuations about zero)
         
          	[~, mink] = min(abs([d0(k), d0(k+1)]));    % ensure correct frequency bin is picked out (find smaller of two values either side of crossing (in original pxx signal))
          	if mink == 1
             	minim2 = k;
           	elseif mink == 2
              	minim2 = k+1;
            end
            
            cnt = cnt+1;            % advance counter by 1
          	posZ2(cnt,1) = cnt;         % zero-crossing count
           	posZ2(cnt,2) = minim2;          % bin for CoG
            posZ2(cnt,3) = f(minim2);        % zero-crossing frequency
          	posZ2(cnt,4) = d0(minim2);        % power estimate (from original signal)
          %  end
        % try to make use of shallow slope in d1 as (early) indicator of peak boundary
        elseif abs(d1(k)) < 1 && d1(k) > d1(k+1) && d1(k+1) > d1(k+2) && d1(k+2) > d1(k+3) && d1(k+3) > d1(k+4)
            %abs(d1(k)) < 1 & abs(d1(k:k+t)) < 1     % range(d1(k:k+4)) < 1
           %% ??need a check that there isn't a more prominant drop off (or indeed, secondary peak) later (i.e confirm that decrease is maintained to some extent) - have cake (cut off shallow slope sans zcross) and eat (don't cut off too soon when there's a peak still to come) 
            minim2 = k;
            cnt = cnt+1;            % advance counter by 1
            posZ2(cnt,1) = cnt;         % zero-crossing count
            posZ2(cnt,2) = minim2;
            posZ2(cnt,3) = f(minim2);        % zero-crossing frequency
            posZ2(cnt,4) = d0(minim2);        % power estimate (from original signal)
        end
        
    end
    
    f2 = posZ2(1, 2);
    posZ2 = posZ2(1, 3);            % just take first estimate for output
    
    

%% define boundaries by inflection points (2nd derivative of smoothed signal) - only for included channels (primary peaks)
    inf1 = zeros(1,3);                     % initialise for zero-crossing count & frequency bin
    cnt = 0;                            % start counter at 0
    for k = 1:peakBin-1     % step through frequency bins in alpha band
        if sign(d2(k)) > sign(d2(k+1))            % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
            %if range(sign(d1sgf(k+1:k+t)))==0           % must remain downgoing for next t frequency bins (exclude noisy fluctuations about zero)
                [~, mink] = min(abs([d2(k), d2(k+1)]));      % ensure correct frequency bin is picked out (find smaller of two values either side of crossing)
                if mink == 1
                    min1 = k;
                elseif mink == 2;
                    min1 = k+1;
                end
                cnt = cnt+1;            % advance counter by 1
                inf1(cnt,1) = cnt;         % zero-crossing count
                inf1(cnt,2) = f(min1);        % zero-crossing frequency
                inf1(cnt,3) = d0(min1);        % power estimate (from original signal)
            %end
        end
    end

    
    % sort out appropriate estimates for output
    if size(inf1, 1) == 1      % if singular crossing --> report frequency
        inf1 = inf1(1, 2);
    elseif size(inf1, 1) >1        % if > 1 crossing detected...
        inf1 = sortrows(inf1, -2);        % sort by frequency values (descending)...
        inf1 = inf1(1, 2);                  % take highest frequency (bin nearest to peak)
    end
    
    
    % look for end of alpha peak
    inf2 = zeros(1,1);                     % initialise for zero-crossing count & frequency bin
    for k = peakBin+1:length(d2)-1     % step through frequency bins in alpha band
        if sign(d2(k)) < sign(d2(k+1))            % look for upward zero-crossing
            %if range(sign(d1sgf(k+1:k+t)))==0           % must remain downgoing for next t frequency bins (exclude noisy fluctuations about zero)
            [~, mink] = min(abs([d2(k), d2(k+1)]));    % ensure frequency bin nearest zero-crossing point picked out (find smaller of two values either side of crossing)
                if mink == 1
                    min2 = k;
                elseif mink == 2
                    min2 = k+1;
                end
                inf2 = f(min2);        % zero-crossing frequency
                break                   % break loop as only need to record first instance of crossing
                
            %end
        end

    end

    
    % estimate approx. area under curve between inflection points either
    % side of peak, scale by inflection band width (similar to G et al rho)
    
   
    Q = trapz(f(min1:min2), d0(min1:min2));
    Qf = Q / (min2-min1);       % narrower bands privileged over broad

end 

end

