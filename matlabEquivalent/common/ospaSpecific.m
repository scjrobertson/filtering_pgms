function [euclideanOspa, hellingerOspa] = ospaSpecific(muX, Sxx, muY, Syy, euclideanC, euclideanP, hellingerC, hellingerP)
% OSPA - Optimal Subpattern Assigment
%   [euclideanOspa, hellingerOspa] = hellingerOspaSpecific(muX, Sxx, muY, Syy, euclideanC, euclideanP, hellingerC, hellingerP)
%
%   Compute the distance between two finite sets using the
%   Schumacher metric. This function calculates both the
%   Euclidean-OSPA and Hellinger-OSPA betwen a two finite sets of
%   Gaussian distrubtions. This is not general, it only works for a 
%   4x4 covariance matrix, but hopefully it is fast.
%
%   See also euclideanOspa, hellingerOspa, Hungarian.
%
%   Inputs
%       muX - (d, n) matrix. A finite set of Gaussian means.
%       Sxx - (d, d, n) matrix. A finite set of Gaussian covariance, in
%           muX's respective order.
%       muy - (d, n) matrix. A finite set of Gaussian means.
%       Syy - (d, d, n) matrix. A finite set of Gaussian covariance, in
%           muY's respective order.
%       euclideanC - double. Cut-off parameter for Euclidean-OSPA.
%       euclideanP - double. p-parameter for the Euclidean-OSPA.
%       hellingerC - double. Cut-off parameter for Hellinger-OSPA.
%       hellingerP - double. p-parameter for the Hellinger-OSPA.
%
%   Outputs
%       euclideanOSPA - (3, 1) array. The Euclidean-OSPA components.
%       hellingerOSPA - (3, 1) array. The Hellinger OSPA components.
%% Case 1 - Both sets are empty
if isempty(muX) && isempty(muY)
    euclideanOspa = zeros(3, 1);
    hellingerOspa = zeros(3, 1);
    return;
end
%% Case 2 - A single set is empty
if isempty(muX) || isempty(muY)
    euclideanOspa = [euclideanC; 0; euclideanC];
    hellingerOspa = [hellingerC; 0; hellingerC];   
    return;
end
%% Case 3 - Non-empty sets
% Calculate sizes of the input point patterns
[d, n] = size(muX);
m = size(muY, 2);
largestCardinality = max(m, n);
cardinalityDifference = abs(m-n);
% Vectorise some operations
muXX = repmat(muX, [1 m]);
muYY = reshape(repmat(muY,[n 1]), [d n*m]);
SXX = repmat(Sxx, [1 1 m]);
SYY = reshape(repmat(reshape(Syy, [d d*m]), [1 n]), [d d n*m]);
difference = muXX - muYY;
%% Euclidean distance
euclideanDistance = reshape(sqrt(sum((difference).^2)), [n m]);
euclideanCutOff = min(euclideanC, euclideanDistance).^euclideanP;
[~, euclideanCost] = Hungarian(euclideanCutOff);
euclideanOspa(1, 1) = ((1/largestCardinality)*( (euclideanC^euclideanP)*cardinalityDifference + euclideanCost) ) ^(1/euclideanP);
euclideanOspa(2, 1) = ((1/largestCardinality)*euclideanCost)^(1/euclideanP);
euclideanOspa(3, 1) = ((1/largestCardinality)*(euclideanC^euclideanP)*cardinalityDifference)^(1/euclideanP);
%% Hellinger distance
muZZ = reshape(difference, [d 1 n m]);
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
hellingerDistance = real(sqrt(1 - delta.*epsilon)); %Hellinger distance
hellingerCutOff = min(hellingerC, hellingerDistance).^hellingerP;
%Compute optimal assignment and cost using the Hungarian algorithm
[~, hellingerCost] = Hungarian(hellingerCutOff);
%Calculate final distance
hellingerOspa(1, 1) = ((1/largestCardinality)*( (hellingerC^hellingerP)*cardinalityDifference + hellingerCost) ) ^(1/hellingerP);
hellingerOspa(2, 1) = ((1/largestCardinality)*hellingerCost)^(1/hellingerP);
hellingerOspa(3, 1) = ((1/largestCardinality)*(hellingerC^hellingerP)*cardinalityDifference)^(1/hellingerP);

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