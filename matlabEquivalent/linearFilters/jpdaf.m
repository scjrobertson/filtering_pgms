function stateEstimates = jpdaf(model, targetPriors, measurements)
% JPDAF -- Estimates the target's states using approximate JPDAF
%   stateEstimates = jpdaf(model, targetPriors, measurements)
%
%   Runs an approximate JPDAF - it resolves the data association by way of
%   loopy belief propagation. This implementation avoid for loops and is
%   not very readable. Using matrices instead of 3D arrays throughout would
%   probably be more efficient.
%
%   See also generateModel, generateGroundTruth, jipdaf, loopyBeliefPropagation
%   and plotResults.
%
%   Inputs
%       model - struct. The struct declared in generateModel, it has fields
%           describing the linear motion and measurement models, etc.
%       targetPriors - struct. A struct with the following fields:
%           birthsTimes - (1, n) array. The targets' birth times.
%           deathTimes - (1, n) array. The targets' death times.
%           means - (d, n) array. The targets' prior means.
%           covariance - (d, d, n) array. The targets' prior
%       measurements - (1, m) cell. Contains the measurements --
%           target-generated and clutter generated measurements -- for each
%           time-step of the simulation.
%   Output
%       stateEstimates - struct. A structure with the following fields.
%           uniqueLables - (1, ) array. A set of all track labels used in
%               the simulation -- for plotting purposes.
%           labels - (1, m) cell. Contains target 'identity' for each
%               track -- for plotting purposes.
%           means - (1, m) cell. Contains the state estimates means.
%           covariances - (1, m) cell. Contains the state estimates
%               covariance matrices.
%           cardinality - (1, m) array. The estimated number of targets present at each
%               time-step of the simulation.
%% Admin
simulationLength = size(measurements, 2);
I = eye(model.xDimension);
xLimits = model.observationSpaceLimits(1, :)';
yLimits = model.observationSpaceLimits(2, :)';
%% Add in initial targets
targetLabels = find(targetPriors.birthTimes == 1);
mu = targetPriors.means(:, targetLabels);
S = targetPriors.covariances(:, :, targetLabels);
targetNumber = size(targetLabels, 2);
%% State estimates
stateEstimates.uniqueLabels = targetLabels;
stateEstimates.labels = cell(1, simulationLength); stateEstimates.labels{1} = targetLabels;
stateEstimates.means = cell(1, simulationLength); stateEstimates.means{1} = mu;
stateEstimates.covariances = cell(1, simulationLength); stateEstimates.covariances{1} = S;
stateEstimates.cardinality = zeros(1, simulationLength); stateEstimates.cardinality(1) = targetNumber;
%% Run the filter
for i = 2:simulationLength
    %% Run the filter, if there are any targets.
    if targetNumber ~= 0
        %% Predict the targets' states
        muPred = model.A*mu + model.u;
        SPred =  quadraticMultiprod(S, model.A, model.Atranspose, model.xDimension, targetNumber) + model.R; % SPred = A*S*A' + R
        %% Create update components
        muZ = model.C*muPred; % Predicted measurement mean
        Kxz = rightMultiprod(SPred, model.Ctranspose, model.xDimension, model.xDimension, model.zDimension, targetNumber); % Cross covariance
        Z = leftMultiprod(model.C, Kxz, model.zDimension, model.xDimension, model.zDimension, targetNumber) + model.Q; % Predicted measurement covariance
        ZInv = simplifiedMultinv(Z, model.zDimension, targetNumber); % Inverse of Z, used to determine likelihoods
        detZ = reshape(Z(1, 1, :).*Z(2, 2, :) - Z(2, 1, :).*Z(1, 2, :), [targetNumber 1]);
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
        difference = reshape(zRep - muZRep, [1 model.zDimension targetNumber numberOfMeasurements]);
        differencePermuted = permute(difference, [2 1 3 4]);
        ZInvShaped = repmat(ZInv, [1 1 1 numberOfMeasurements]);
        leftProduct = sum(bsxfun(@times, ZInvShaped, difference), 2);
        rightProduct = reshape(sum(bsxfun(@times, differencePermuted, leftProduct), 1), [targetNumber numberOfMeasurements]);
        associationMatrix(:, 2:end) = repmat(normalisingConstants, [1 numberOfMeasurements]).*exp(-0.5*rightProduct);
        %% Loopy Belief Propagation
        clutterLikelihoods = ones(1, numberOfMeasurements)/model.observationSpaceVolume;
        [updatedAssociationMatrix, ~] = loopyBeliefPropagation(associationMatrix, clutterLikelihoods, 10e-6, 5); % JPDAF
        %associationMatrix(:, 1) = associationMatrix(:, 1)./model.observationSpaceVolume;
        %updatedAssociationMatrix = associationMatrix./sum(associationMatrix, 2); %Multiple instance of PDAFs
        %% Update states
        % Means
        associationProbabilties =  permute(reshape(updatedAssociationMatrix, [1 targetNumber numberOfMeasurements+1]), [1 3 2]);
        missedDetectionProbabilities = associationProbabilties(1, 1, :);
        muUpdated = zeros([model.xDimension 1 targetNumber numberOfMeasurements+1]);
        muUpdated(:, :, :, 1) = reshape(muPred, [model.xDimension 1 targetNumber]);
        muUpdated(:, :, :, 2:end) = muUpdated(:, :, :, 1) + sum(bsxfun(@times, K, difference), 2);
        if targetNumber == 1; mu = sum(bsxfun(@times, squeeze(muUpdated), associationProbabilties), 2); 
        else; mu = sum(bsxfun(@times, permute(squeeze(muUpdated), [1 3 2]), associationProbabilties), 2); end
        % Covariance matrices
        muShift = permute(muUpdated - mu, [1 2 4 3]);
        outerProduct = bsxfun(@times, muShift, permute(muShift, [2 1 3 4]));
        associationProbabiltiesReshaped = reshape(associationProbabilties, [1 1 numberOfMeasurements+1 targetNumber]);
        rankOneMatrices = squeeze(sum(bsxfun(@times, associationProbabiltiesReshaped, outerProduct), 3));
        positiveDefiniteComponents = bsxfun(@times, missedDetectionProbabilities, SPred) + bsxfun(@times, 1- missedDetectionProbabilities, SUpdated);
        %% Reassign
        mu = reshape(mu, [model.xDimension targetNumber]);
        S = rankOneMatrices + positiveDefiniteComponents;
    end
    %% Add in new targets
    newTargetLabels = find(targetPriors.birthTimes == i);
    numberOfNewTargets = size(newTargetLabels, 2);
    if numberOfNewTargets ~= 0
        newIndices = (targetNumber+1):(targetNumber + numberOfNewTargets);
        targetLabels(newIndices) = newTargetLabels;
        mu(:, newIndices) = targetPriors.means(:, newTargetLabels);
        S(:, :, newIndices) = targetPriors.covariances(:, :, newTargetLabels);
        targetNumber = targetNumber + numberOfNewTargets;
    end
    %% Remove dead targets
    aliveTargetIndices = targetPriors.deathTimes(targetLabels) >= i; 
    if targetNumber ~= 0
        targetLabels = targetLabels(aliveTargetIndices);
        mu = mu(:, aliveTargetIndices);
        S = S(:, :, aliveTargetIndices);
        targetNumber = size(targetLabels, 2);
    end 
    %% Save the data
    stateEstimates.uniqueLabels = union(stateEstimates.uniqueLabels, targetLabels);
    stateEstimates.labels{i} = targetLabels;
    stateEstimates.means{i} = mu;
    stateEstimates.covariances{i} = S;
    stateEstimates.cardinality(i) = targetNumber;
end
end