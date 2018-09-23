function [paf, cog] = meanIAF(sums, nchan, cmin)
% Average repeated-measures PAF and CoG means in proportion to number of
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
%   nchan = vector containing number of channels per recording
%   cmin = minimum number of channel estimates required for inclusion

%% setup variable check
if ~exist('sums', 'var')
    error('Provide structure containing PAF and CoG estimate fields');
end
if ~exist('nchan', 'var')
    error('Provide vector containing number of channels in each recording')
end
%%
  
% peaks
select = [sums.pSel] >= cmin;
if sum(select) == 0
    paf = NaN;
else
    paf = sum([sums(select).paf].*([sums(select).pSel]./nchan(select)))...
        / sum([sums(select).pSel]./nchan(select));
end
    
% gravs
select = [sums.gSel] >= cmin;
if sum(select) == 0
    cog = NaN;
else
    cog = sum([sums(select).cog].*([sums(select).gSel]./nchan(select)))...
        / sum([sums(select).gSel]./nchan(select));
end


end

