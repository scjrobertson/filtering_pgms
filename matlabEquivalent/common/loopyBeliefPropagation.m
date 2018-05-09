function [pdafUpdated, clutterUpdated] = loopyBeliefPropagation(pdafLikelihoods, clutterLikelihoods,...
    covergenceTolerance, maximumNumberOfIterations)
% LOOPYBELIEFPROPAGATION -- Resolves the data association by loopy belief propagation
%   [pdafUpdated, clutterUpdated] = loopyBeliefPropagation(pdafLikelihoods,
%   clutterLikelihoods, covergenceTolerance, maximumNumberOfIterations)
%
%   Code originally written by Jason Williams; however, it has been
%   (mostly) vectorised.
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
pdafUpdated = zeros(numberOfTargets, numberOfUpdateOptions);
if numberOfTargets == 0
   clutterUpdated = clutterLikelihoods./sum(clutterLikelihoods, 2); 
   return;
end
if numberOfMeasurements == 0 && numberOfTargets > 0
   pdafUpdated(:, 1) = ones(numberOfTargets, 1);
   clutterUpdated = zeros(1, numberOfMeasurements);
   return;
end
%% Continue if there are targets
mu = ones(numberOfTargets, numberOfMeasurements);
muOld = zeros(numberOfTargets, numberOfMeasurements);
nu = zeros(numberOfTargets, numberOfMeasurements);
%% Run loopy belief update propagation
counter = 1;
while(max(abs(mu(:) - muOld(:))) > covergenceTolerance && counter < maximumNumberOfIterations)
    muOld = mu;
    %% Update detection likelihoods
    detectionLikelihoods = pdafLikelihoods(:, 2:end).*mu;
    normalisingConstants = pdafLikelihoods(:, 1) + sum(detectionLikelihoods, 2);
    nu = pdafLikelihoods(:, 2:end)./(normalisingConstants - detectionLikelihoods);
    %% Update clutter likelihoods
    totalClutterProbabilities = clutterLikelihoods + sum(nu, 1);
    mu = 1./(totalClutterProbabilities - nu);
    %% Reset counter
    counter = counter + 1;
end
%% Calculate final probabilities
detectionLikelihoods = pdafLikelihoods(:, 2:end).*mu;
normalisingConstants = pdafLikelihoods(:, 1) + sum(detectionLikelihoods, 2);
pdafUpdated(:, 1) = pdafLikelihoods(:, 1)./normalisingConstants;
pdafUpdated(:, 2:end) = detectionLikelihoods./normalisingConstants;

totalClutterProbabilities = clutterLikelihoods + sum(nu, 1);
clutterUpdated = clutterLikelihoods./totalClutterProbabilities;
end