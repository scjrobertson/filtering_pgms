function [dist, varargout] = hellingerOspa(muX, Sxx, muY, Syy, c, p)
% OSPA - Optimal Subpattern Assigment
%   [dist, varargout] = ospa(X, Y, c, p)
%
%   Compute the distance between two finite sets using the
%   Schumacher metric. This is a carbon copy of B. Vo's original 
%   implementation.
%
%   See also omat, run_gmm_phd.
%
%   Inputs
%       muX - (d, n) matrix. A finite set of Gaussian means.
%       Sxx - (d, d, n) matrix. A finite set of Gaussian covariance, in
%           muX's respective order.
%       muy - (d, n) matrix. A finite set of Gaussian means.
%       Syy - (d, d, n) matrix. A finite set of Gaussian covariance, in
%           muY's respective order.
%       c - double. Cut-off parameter.
%       p - double. p-parameter of the metric.
%
%   Outputs
%       dist - double. The measured distance between two sets.

%% Admin - Valid output number
if nargout ~= 1 && nargout ~=3
   error('Incorrect number of outputs'); 
end
%% Case 1 - Both sets are empty
if isempty(muX) && isempty(muY)
    dist = 0;

    if nargout == 3
        varargout(1)= {0};
        varargout(2)= {0};
    end
    
    return;
end
%% Case 2 - A single set is empty
if isempty(muX) || isempty(muY)
    dist = c;

    if nargout == 3
        varargout(1)= {0};
        varargout(2)= {c};
    end
    
    return;
end
%% Case 3 - Non-empty sets
%Calculate sizes of the input point patterns
n = size(muX, 2);
m = size(muY, 2);
distance = zeros(n, m);

for i = 1:n
   for j = 1:m
        distance = hellingerDistance(muX(:, i), Sxx(:, :, i), muY(:, j), Syy(:, :, j));
   end
end
distance = min(c, distance).^p;

%Compute optimal assignment and cost using the Hungarian algorithm
[~, cost] = Hungarian(distance);

%Calculate final distance
dist = ( 1/max(m,n)*( c^p*abs(m-n) + cost ) ) ^(1/p);

%Output components if called for in varargout
if nargout == 3
    varargout(1)= {(1/max(m,n)*cost)^(1/p)};
    varargout(2)= {(1/max(m,n)*c^p*abs(m-n))^(1/p)};
end
end