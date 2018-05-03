function [pdafUpdated, clutterUpdated] = loopyBeliefPropagation(pdafLikelihoods, clutterLikelihoods,...
    covergenceTolerance, maximumNumberOfIterations)
% LOOPYBELIEFPROPAGATION -- Resolves the data association by loopy belief propagation
%   [pdafUpdated, clutterUpdated] = loopyBeliefPropagation(pdafLikelihoods,
%   clutterLikelihoods, covergenceTolerance, maximumNumberOfIterations)
%
%   Code originally written by Jason Williams.
%
%   Inputs
%       pdafLikelihoods - (n, m+1) matrix. The pdaf likelihoods for each of n
%           targets casuing the m measurements. The first column is the
%           missed detection probability.
%       clutterLikelihood - (m, 1) vector. The likelihood of the clutter
%           causing the m measurements.
%       convergenceTolerance - The convergence threshold for the loopy
%           belief 
%       maximumNumberOfIteration - The maximum number of loopy belief
%           propagation iterations -- just to avoid non convergence.
%
%   Outputs
%       pdafUpdated - (n, m+1) matrix. The resolved pdaf probabilities.
%       clutterUpdated - (n, m) vector. The resolved pdaf probabilities.
%% Determine shape of the data
[numberOfTargets, numberOfUpdateOptions] = size(pdafLikelihoods);
numberOfMeasurements = numberOfUpdateOptions - 1;
%% Declare the variables
mu = ones(numberOfTargets, numberOfMeasurements);
muOld = zeros(numberOfTargets, numberOfMeasurements);
nu = zeros(numberOfTargets, numberOfMeasurements);

pdafUpdated = zeros(numberOfTargets, numberOfUpdateOptions);
clutterUpdated = zeros(numberOfMeasurements, 1);
%% Run loopy belief update propagation
counter = 1;
while( max(abs(mu(:) - muOld(:))) > covergenceTolerance && counter < maximumNumberOfIterations)
    muOld = mu;
    
    for i = 1:numberOfTargets
        detectionLikelihood = pdafLikelihoods(i, 2:end).*mu(i, :);
        normalisingConstant = pdafLikelihoods(i, 1) + sum(detectionLikelihood);
        nu(i, :) = pdafLikelihoods(i, 2:end)./(normalisingConstant - detectionLikelihood);
    end
    
    
    for i = 1:numberOfMeasurements
        normalisingConstant = clutterLikelihoods(i) + sum(nu(:, i));
        mu(:, i) = 1./(normalisingConstant - nu(:, i));
    end
    
    counter = counter + 1;
end
%% Calculate final probabilities
for i = 1:numberOfTargets
   normalisingConstant = pdafLikelihoods(i, 1) + sum(pdafLikelihoods(i, 2:end).*mu(i, :));
   pdafUpdated(i, 1) = pdafLikelihoods(i, 1)/normalisingConstant;
   pdafUpdated(i, 2:end) = pdafLikelihoods(i, 2:end).*mu(i, :)/normalisingConstant;
end

for i = 1:numberOfMeasurements
   normalisingConstant = clutterLikelihoods(i) + sum(nu(:, i)); 
   clutterUpdated(i) = pdafUpdated(i)/normalisingConstant;
end
end