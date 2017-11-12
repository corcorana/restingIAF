function tval = lessThan1(d1)
% When encounter 1st derivative absolute value < 1, eval whether following
% values (within segment) remain < +/- 1 
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Output:
%   tval = logical
%
%% Required input:
%   d1 = segment of 1st derivative spanning approx. 1 Hz

%% setup variable check
if ~exist('d1', 'var')
    error('Provide vector containing series of 1st derivative values for slope analysis')
elseif length(d1) < 2
    error('Length of 1st derivative segment < 2');
end
%%

t = zeros(1, length(d1));
for kx = 1:length(d1)
    t(kx) = abs(d1(kx) < 1);
end

if all(t == 1)
    tval = 1;
else
    tval = 0;
end

end