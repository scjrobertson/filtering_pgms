function [stateEstimates, cardinalityEstimates] = runFilter(model, targetPriors, measurements)
%% Admin
simulationLength = size(measurements, 2);
I = eye(model.xDimension);
%% Add in initial targets
targetLabels = find(targetPriors.birthTimes == 1);
mu = targetPriors.means(:, targetLabels);
S = targetPriors.covariances(:, :, targetLabels);
targetNumber = size(targetLabels, 2);
%% Run the filter
for i = 2:2%simulationLength
    %% Predict the targets' states
    muPred = model.A*mu + model.u;
    SPred =  quadraticMultiprod(S, model.A, model.Atranspose, model.xDimension, targetNumber) + model.R; % SPred = A*S*A' + R
    %% Create update components
    muZ = model.C*muPred; % Predicted measurement mean
    Kxz = rightMultiprod(SPred, model.Ctranspose, model.xDimension, model.xDimension, model.zDimension, targetNumber); % Cross covariance
    Z = leftMultiprod(model.C, Kxz, model.zDimension, model.xDimension, model.zDimension, targetNumber) + model.Q; % Predicted measurement covariance
    ZInv = simplifiedMultinv(Z, model.zDimension, targetNumber); % Inverse of Z, used to determine likelihoods
    detZ = squeeze(Z(1, 1, :).*Z(2, 2, :) - Z(2, 1, :).*Z(1, 2, :)); % Determinants
    normalisingConstants = (model.detectionProbability)./sqrt((2*pi)^(model.zDimension)*detZ); % Normalising constants for likelihoods
    K = simplifiedMultiprod(Kxz, ZInv, model.xDimension, model.zDimension, model.zDimension, targetNumber); % Kalman gain
    STemp = rightMultiprod(K, model.C, model.xDimension, model.zDimension, model.xDimension, targetNumber);
    SUpdated = simplifiedMultiprod((I - STemp), SPred, model.xDimension, model.xDimension, model.xDimension, targetNumber); % Update covariance matrices
    %% Update the targets states
    numberOfMeasurements = size(measurements{i}, 2);
    associationMatrix = (1-model.detectionProbability)*ones(targetNumber, numberOfMeasurements+1);
    muZRep = repmat(muZ, [1 numberOfMeasurements]);
    zRep = reshape(repmat(measurements{i},[targetNumber 1]),[model.zDimension targetNumber*numberOfMeasurements]);
    % Calculate the likelihoods
    difference = reshape(zRep - muZRep, [model.zDimension 1 targetNumber numberOfMeasurements]);
    ZInvStacked = repmat(reshape(ZInv, [model.zDimension model.zDimension*targetNumber]), [1 numberOfMeasurements]);
    ZInvShaped = reshape(ZInvStacked, [model.zDimension model.zDimension targetNumber numberOfMeasurements]);
    leftProduct = sum(bsxfun(@times, ZInvShaped, difference), 2);
    rightProduct = sum(difference.*leftProduct, 1);
    associationMatrix(:, 2:end) = exp(-0.5*squeeze(rightProduct)).*repmat(normalisingConstants, [1 numberOfMeasurements]);
    %% Loopy Belief Propagation
    clutterLikelihoods = ones(1, numberOfMeasurements)/model.observationSpaceVolume;
    [updatedAssociationMatrix, ~] = loopyBeliefPropagation(associationMatrix, clutterLikelihoods, 10e-6, 200);
    %% Updated States
    missedDetectionProbabilities = reshape(updatedAssociationMatrix(:, 1), [1 1 targetNumber]);
    % Means
    muPredReshaped = reshape(muPred, [model.xDimension 1 targetNumber]);
    innovation = sum(bsxfun(@times, K, permute(difference, [2 1 3 4])), 2);
    muUpdated = muPredReshaped + reshape(innovation, [model.xDimension numberOfMeasurements targetNumber]);
    mu = bsxfun(@times, missedDetectionProbabilities, muPredReshaped) + sum(bsxfun(@times, updatedAssociationMatrix(:, 2:end), muUpdated), 2);
    % Covariance matrices
    postiveDefiniteComponent = bsxfun(@times, missedDetectionProbabilities, SPred) + bsxfun(@times, 1- missedDetectionProbabilities, SUpdated);
    %% Reassign
    mu = reshape(mu, [model.xDimension targetNumber]);
    S = SPred
end
stateEstimates = 0;
cardinalityEstimates = zeros(1, simulationLength);
end

%rProbPacked = reshape(rProbability, [1 1 numberOfMeasurements]);
%difference = jpdafTrackerManager.x(:, i) - jpdafTrackerManager.updatedX(:, i, :);
%differencePermuted = permute(difference, [2 1 3]);
%rankOneMatrices = sum(multiprod(multiprod(difference, differencePermuted), rProbPacked), 3);
%positiveDefiniteMatrices = rProbability(1)*jpdafTrackerManager.updatedS(:, :, i, 1) + (1-rProbability(1))*jpdafTrackerManager.updatedS(:, :, i, 2);
%jpdafTrackerManager.S(:, :, i) = rankOneMatrices + positiveDefiniteMatrices;