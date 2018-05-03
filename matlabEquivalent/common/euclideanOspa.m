function [dist, varargout] = euclideanOspa(X, Y, c, p)
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
%       X - (K, N) matrix. A finite set of theoretical target locations.
%       Y - (K, M) matrix. A finite set of measured target locations.
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
if isempty(X) && isempty(Y)
    dist = 0;

    if nargout == 3
        varargout(1)= {0};
        varargout(2)= {0};
    end
    
    return;
end
%% Case 2 - A single set is empty
if isempty(X) || isempty(Y)
    dist = c;

    if nargout == 3
        varargout(1)= {0};
        varargout(2)= {c};
    end
    
    return;
end
%% Case 3 - Non-empty sets
%Calculate sizes of the input point patterns
n = size(X, 2);
m = size(Y, 2);

%Determine cost matrices
XX = repmat(X,[1 m]);
YY = reshape(repmat(Y,[n 1]),[size(Y,1) n*m]);
D = reshape(sqrt(sum((XX-YY).^2)),[n m]);
D = min(c,D).^p;

%Compute optimal assignment and cost using the Hungarian algorithm
[~, cost] = Hungarian(D);

%Calculate final distance
dist = ( 1/max(m,n)*( c^p*abs(m-n) + cost ) ) ^(1/p);

%Output components if called for in varargout
if nargout == 3
    varargout(1)= {(1/max(m,n)*cost)^(1/p)};
    varargout(2)= {(1/max(m,n)*c^p*abs(m-n))^(1/p)};
end
end