function [paf, cog] = meanIAF(sums, nchan)
% Average 2 repeated-measures PAF and CoG means in proportion to number of
% channels that contributed to calculation of each respective mean.
% No averaging performed if only one valid estimate provided per subject.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
% Outputs:
%   paf = grand average peak alpha frequency
%   cog = grand average centre of gravity
%
% Required inputs:
%   sums = structure containing mean PAF/CoG estimates and number of
%          channels from which estimates derived
%   nchan = number of channels included in analysis (all available channels)

%% setup variable check
if ~exist('sums', 'var')
    error('Provide structure containing PAF and CoG estimate fields');
end
if ~exist('nchan', 'var')
    error('Provide vector containing number of channels in each recording')
elseif length(nchan) ~= 2
    error('Length of nchan vector ~= 2');
end
%%

if length(sums) == 2
    
    % do peaks
    if isnan(sums(1).paf) && isnan(sums(2).paf)
        paf = NaN;
    elseif isnan(sums(1).paf) || isnan(sums(2).paf)
        paf = nansum([sums(1).paf, sums(2).paf]);
    else
        paf = sum([sums(1).paf * (sums(1).pSel/nchan(1)), sums(2).paf * (sums(2).pSel/nchan(2))]) / ((sums(1).pSel/nchan(1)) + (sums(2).pSel/nchan(2)));
    end
    
    % do gravs
    if isnan(sums(1).cog) && isnan(sums(2).cog)
        cog = NaN;
    elseif isnan(sums(1).cog) || isnan(sums(2).cog)
        cog = nansum([sums(1).cog, sums(2).cog]);
    else
        cog = sum([sums(1).cog * (sums(1).gSel/nchan(1)), sums(2).cog * (sums(2).gSel/nchan(2))]) / ((sums(1).gSel/nchan(1)) + (sums(2).gSel/nchan(2)));
    end
    
elseif length(sums) == 1        % if only one set of estimates, cannot average
    paf = sums.paf;
    cog = sums.cog;

else
    paf = NaN;
    cog = NaN;
end

end

