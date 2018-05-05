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
% Calculate sizes of the input point patterns
[d, n] = size(muX);
m = size(muY, 2);
% Vectorise some operations
muXX = repmat(muX, [1 m]);
muYY = reshape(repmat(muY,[n 1]), [d n*m]);
SXX = repmat(Sxx, [1 1 m]);
SYY = reshape(repmat(reshape(Syy, [d d*n]), [1 m]), [d d n*m]);
muZZ = reshape(muXX - muYY, [d 1 n m]);
SZZ = SXX + SYY;
% To speed things up determinants are calculated by "hand".
detX = reshape(fourByFourDeterminant(SXX), [n m]);
detY = reshape(fourByFourDeterminant(SYY), [n m]);
detZ = reshape(fourByFourDeterminant(SZZ), [n m]);
% Calculate inverses
SZZInv = reshape(simplifiedMultinv(SZZ, d, n*m), [d d n m]); % Inverse of SZZ
leftProduct = sum(bsxfun(@times, SZZInv, muZZ), 2); % Next two lines calculate muZ'*SZZ*muZ
rightProduct = reshape(sum(muZZ.*leftProduct, 1), [n m]); 
% Calculate Hellinger distance components
delta = sqrt(sqrt(detX.*detY)./((0.5.^d)*detZ));
epsilon = exp(-0.25*rightProduct);
distance = real(sqrt(1 - delta.*epsilon)); %Hellinger distance
%Compute optimal assignment and cost using the Hungarian algorithm
[~, cost] = Hungarian(distance);

%Calculate final distance
dist = ( 1/max(m,n)*( c^p*abs(m-n) + cost ) ) ^(1/p);

%Output components if called for in varargout
if nargout == 3
    varargout(1)= {(1/max(m,n)*cost)^(1/p)};
    varargout(2)= {(1/max(m,n)*c^p*abs(m-n))^(1/p)};
end

    function determinant = fourByFourDeterminant(S)
        % FOURBYFOURDETERMINANT - Determines a 4x4xn array's determinant
        explicit = S(1, 1, :).*S(2, 2, :).*S(3, 3, :).*S(4, 4, :) ... % New block
            - S(1, 1, :).*S(2, 2, :).*S(3, 4, :).*S(4, 3, :) ...
            - S(1, 1, :).*S(3, 2, :).*S(2, 3, :).*S(4, 4, :) ...
            + S(1, 1, :).*S(3, 2, :).*S(2, 4, :).*S(4, 3, :) ...
            + S(1, 1, :).*S(4, 2, :).*S(2, 3, :).*S(3, 4, :) ...
            - S(1, 1, :).*S(4, 2, :).*S(2, 4, :).*S(3, 3, :) ...
            - S(2, 1, :).*S(1, 2, :).*S(3, 3, :).*S(4, 4, :) ... % New block
            + S(2, 1, :).*S(1, 2, :).*S(3, 4, :).*S(4, 3, :) ...
            + S(2, 1, :).*S(3, 2, :).*S(1, 3, :).*S(4, 4, :) ...
            - S(2, 1, :).*S(3, 2, :).*S(1, 4, :).*S(4, 3, :) ...
            - S(2, 1, :).*S(4, 2, :).*S(1, 3, :).*S(3, 4, :) ...
            + S(2, 1, :).*S(4, 2, :).*S(1, 4, :).*S(3, 3, :) ...
            + S(3, 1, :).*S(1, 2, :).*S(2, 3, :).*S(4, 4, :) ... % New block
            - S(3, 1, :).*S(1, 2, :).*S(2, 4, :).*S(4, 3, :) ...
            - S(3, 1, :).*S(2, 2, :).*S(1, 3, :).*S(4, 4, :) ...
            + S(3, 1, :).*S(2, 2, :).*S(1, 4, :).*S(4, 3, :) ...
            + S(3, 1, :).*S(4, 2, :).*S(1, 3, :).*S(2, 4, :) ...
            - S(3, 1, :).*S(4, 2, :).*S(1, 4, :).*S(2, 3, :) ...
            - S(4, 1, :).*S(1, 2, :).*S(1, 3, :).*S(2, 4, :) ... % New block
            + S(4, 1, :).*S(1, 2, :).*S(1, 4, :).*S(2, 3, :) ...
            + S(4, 1, :).*S(2, 2, :).*S(1, 3, :).*S(3, 4, :) ...
            - S(4, 1, :).*S(2, 2, :).*S(1, 4, :).*S(3, 3, :) ...
            - S(4, 1, :).*S(3, 2, :).*S(1, 3, :).*S(2, 4, :) ...
            + S(4, 1, :).*S(3, 2, :).*S(1, 4, :).*S(2, 3, :);
        determinant = squeeze(explicit);
    end
end