function stateEstimates = jipdaf(model, measurements)
% JIPDAF -- Estimates the target's states using approximate JPDAF
%   stateEstimates = jipdaf(model, targetPriors, measurements)
%
%   Runs an approximate JIPDAF - it resolves the data association by way of
%   loopy belief propagation. This implementation avoid for loops and is
%   not very readable. Using matrices instead of 3D arrays throughout would
%   probably be more efficient.
%
%   See also generateModel, generateGroundTruth, jpdaf, loopyBeliefPropagation
%   and plotResults.
%
%   Inputs
%       model - struct. The struct declared in generateModel, it has fields
%           describing the linear motion and measurement models, etc.
%       measurements - (1, m) cell. Contains the measurements --
%           target-generated and clutter generated measurements -- for each
%           time-step of the simulation.
%   Output
%       stateEstimates - struct. A structure with the following fields.
%           means - (1, m) cell. Contains the state estimates means.
%           covariances - (1, m) cell. Contains the state estimates
%               covariance matrices.
%           cardinality - (1, m) array. The estimated number of targets present at each
%               time-step of the simulation.
%% Admin
simulationLength = size(measurements, 2);
I = eye(model.xDimension);
%% State estimates
stateEstimates.means = cell(1, simulationLength);
stateEstimates.covariances = cell(1, simulationLength);
stateEstimates.cardinality = zeros(1, simulationLength);
%% Declare state variables
r = zeros(0); mu = zeros(model.xDimension, 0); S = zeros(model.xDimension, model.xDimension, 0); % Established targets' states
lambdaU = zeros(0); xU = zeros(model.xDimension, 0); SU = zeros(model.xDimension, model.xDimension, 0); % Poisson point process parameters
targetNumber = 0;
%% Run the filter
for i = 2:simulationLength
    %% Thin, then predict the Poisson Point Process
    lambdaU = model.PoissonSurvivalProbability*lambdaU;
    survivingPoissonIndices = lambdaU > model.lambdaThreshold;
    lambdaU = lambdaU(survivingPoissonIndices);
    xU = model.A*xU(:, survivingPoissonIndices) + model.u;
    SU = quadraticMultiprod(SU(:, :, survivingPoissonIndices), model.A, model.Atranspose, model.xDimension, targetNumber) + model.R;
    %% Add tentative spawned targets to the PPP
    intensitySize = size(lambdaU, 2) + model.numberOfSpawningLocations;
    lambdaU(end+1:intensitySize) = model.newTargetProbability;
    xU(:, end+1:intensitySize) = model.spawnMeans;
    SU(:, :, end+1:intensitySize) = model.spawnCovariances;
    %% Start update
    numberOfMeasurements = size(measurements{i}, 2);
    if (targetNumber ~= 0)
        %% Predict the targets' states
        rPred = model.survivalProbability*r; % Exclude this if model.survivalProbability = 1
        muPred = model.A*mu + model.u;
        SPred =  quadraticMultiprod(S, model.A, model.Atranspose, model.xDimension, targetNumber) + model.R; % SPred = A*S*A' + R
        %% Create update components for existing targets
        muZ = model.C*muPred; % Predicted measurement mean
        Kxz = rightMultiprod(SPred, model.Ctranspose, model.xDimension, model.xDimension, model.zDimension, targetNumber); % Cross covariance
        Z = leftMultiprod(model.C, Kxz, model.zDimension, model.xDimension, model.zDimension, targetNumber) + model.Q; % Predicted measurement covariance
        ZInv = simplifiedMultinv(Z, model.zDimension, targetNumber); % Inverse of Z, used to determine likelihoods
        detZ = reshape(Z(1, 1, :).*Z(2, 2, :) - Z(2, 1, :).*Z(1, 2, :), [targetNumber 1]);
        normalisingConstants = (model.detectionProbability)./sqrt((2*pi)^(model.zDimension)*detZ); % Normalising constants for likelihoods
        K = simplifiedMultiprod(Kxz, ZInv, model.xDimension, model.zDimension, model.zDimension, targetNumber); % Kalman gain
        STemp = rightMultiprod(K, model.C, model.xDimension, model.zDimension, model.xDimension, targetNumber);
        SUpdated = simplifiedMultiprod((I - STemp), SPred, model.xDimension, model.xDimension, model.xDimension, targetNumber); % Update covariance matrices
        % Create updated state components
        associationMatrix = repmat(1 - rPred + rPred*(1-model.detectionProbability), [1 numberOfMeasurements +1]);
        rUpdated = ones(targetNumber, numberOfMeasurements+1); rUpdated(:, 1) = rPred./sum(associationMatrix(:, 1));
        muZRep = repmat(muZ, [1 numberOfMeasurements]);
        zRep = reshape(repmat(measurements{i},[targetNumber 1]),[model.zDimension targetNumber*numberOfMeasurements]);
        % Calculate the likelihoods
        difference = reshape(zRep - muZRep, [model.zDimension 1 targetNumber numberOfMeasurements]);
        ZInvStacked = repmat(reshape(ZInv, [model.zDimension model.zDimension*targetNumber]), [1 numberOfMeasurements]);
        ZInvShaped = reshape(ZInvStacked, [model.zDimension model.zDimension targetNumber numberOfMeasurements]);
        leftProduct = sum(bsxfun(@times, ZInvShaped, difference), 2);
        rightProduct = reshape(sum(difference.*leftProduct, 1), [targetNumber numberOfMeasurements]);
        associationMatrix(:, 2:end) = repmat(normalisingConstants, [1 numberOfMeasurements]).*exp(-0.5*rightProduct);
    end
    %% Create update components for the PPP
    % Create update components
    muZU = model.C*xU; % Predicted measurement mean
    KUxz = rightMultiprod(SU, model.Ctranspose, model.xDimension, model.xDimension, model.zDimension, intensitySize); % Cross covariance
    ZU = leftMultiprod(model.C, KUxz, model.zDimension, model.xDimension, model.zDimension, intensitySize) + model.Q; % Predicted measurement covariance
    ZUInv = simplifiedMultinv(ZU, model.zDimension, intensitySize); % Inverse of Z, used to determine likelihoods
    detZU = reshape(ZU(1, 1, :).*ZU(2, 2, :) - ZU(2, 1, :).*ZU(1, 2, :), [intensitySize 1]);
    normalisingConstants = (model.detectionProbability)./sqrt((2*pi)^(model.zDimension)*detZU);
    KU = simplifiedMultiprod(KUxz, ZUInv, model.xDimension, model.zDimension, model.zDimension, intensitySize); % Kalman gain
    SUTemp = rightMultiprod(KU, model.C, model.xDimension, model.zDimension, model.xDimension, intensitySize);
    SUUpdated = simplifiedMultiprod((I - SUTemp), SU, model.xDimension, model.xDimension, model.xDimension, intensitySize); % Update covariance matrices
    % Create update state components
    muZURep = repmat(muZU, [1 numberOfMeasurements]);
    zRep = reshape(repmat(measurements{i},[intensitySize 1]),[model.zDimension intensitySize*numberOfMeasurements]);
    % Determine likelihoods
    intensityDifference = reshape(zRep - muZURep, [model.zDimension 1 intensitySize numberOfMeasurements]);
    ZUInvStacked = repmat(reshape(ZUInv, [model.zDimension model.zDimension*intensitySize]), [1 numberOfMeasurements]);
    ZUInvShaped = reshape(ZUInvStacked, [model.zDimension model.zDimension intensitySize numberOfMeasurements]);
    leftProduct = sum(bsxfun(@times, ZUInvShaped, intensityDifference), 2);
    rightProduct = reshape(sum(intensityDifference.*leftProduct, 1), [intensitySize numberOfMeasurements]);
    clutterLikelihoodMatrix = repmat(normalisingConstants, [1 numberOfMeasurements]).*exp(-0.5*rightProduct);
    intensityLikelihoods = sum(clutterLikelihoodMatrix, 1) + model.clutterPerUnitVolume;
    %% Update the PPP states
    % Determine means
    intensityLikelihoods =  clutterLikelihoodMatrix./intensityLikelihoods;
    intensityProbabilities = reshape(intensityLikelihoods, [1 intensitySize numberOfMeasurements]);
    xUUpdated = reshape(xU, [model.xDimension 1 intensitySize]) + sum(bsxfun(@times, KU, permute(intensityDifference, [2 1 3 4])), 2);
    xUUpdated = permute(xUUpdated, [1 2 4 3]);
    if intensitySize == 1; xNew = sum(bsxfun(@times, squeeze(xUUpdated), intensityProbabilities), 2);
    else; xNew = sum(bsxfun(@times, permute(squeeze(xUUpdated), [1 3 2]), intensityProbabilities), 2); end
    % Determine covariance matrices
    xUShift = xUUpdated - xNew;
    outerProduct = bsxfun(@times, xUShift, permute(xUShift, [2 1 3 4]));
    intensityProbabilitiesReshaped = reshape(intensityProbabilities, [1 1 numberOfMeasurements intensitySize]);
    xNew = reshape(xNew, [model.xDimension numberOfMeasurements]);
    SNew = SUUpdated + squeeze(sum(bsxfun(@times, intensityProbabilitiesReshaped, outerProduct), 3));
    %% Thin the PPP
    lambdaU = (1-model.detectionProbability)*lambdaU;
    significantPoissonIndices = lambdaU > model.lambdaThreshold;
    lambdaU = lambdaU(significantPoissonIndices);
    xU = xU(:, significantPoissonIndices);
    SU = SU(:, :, significantPoissonIndices);
    %% Resolve data association by way of loopy belief propagation
    if (targetNumber ~= 0)
        [updatedAssociationMatrix, intensityLikelihoods] = loopyBeliefPropagation(associationMatrix, intensityLikelihoods, 10e-6, 5);
        %% Update states and compute weak marginals
        % Means
        associationProbabilties =  permute(reshape(updatedAssociationMatrix, [1 targetNumber numberOfMeasurements+1]), [1 3 2]);
        missedDetectionProbabilities = associationProbabilties(1, 1, :);
        muUpdated = zeros([model.xDimension 1 targetNumber numberOfMeasurements+1]);
        muUpdated(:, :, :, 1) = reshape(muPred, [model.xDimension 1 targetNumber]);
        muUpdated(:, :, :, 2:end) = muUpdated(:, :, :, 1) + sum(bsxfun(@times, K, permute(difference, [2 1 3 4])), 2);
        if targetNumber == 1; mu = sum(bsxfun(@times, squeeze(muUpdated), associationProbabilties), 2);
        else; mu = sum(bsxfun(@times, permute(squeeze(muUpdated), [1 3 2]), associationProbabilties), 2); end
        % Covariance matrices
        muShift = permute(muUpdated - mu, [1 2 4 3]);
        outerProduct = bsxfun(@times, muShift, permute(muShift, [2 1 3 4]));
        associationProbabiltiesReshaped = reshape(associationProbabilties, [1 1 numberOfMeasurements+1 targetNumber]);
        rankOneMatrices = squeeze(sum(bsxfun(@times, associationProbabiltiesReshaped, outerProduct), 3));
        positiveDefiniteComponents = bsxfun(@times, missedDetectionProbabilities, SPred) + bsxfun(@times, 1- missedDetectionProbabilities, SUpdated);
        % Reshape and reassign
        mu = reshape(mu, [model.xDimension targetNumber]);
        S = rankOneMatrices + positiveDefiniteComponents;
    end
    %% Add in new tentative tracks
    %% Reassign
    mu = reshape(mu, [model.xDimension targetNumber]);
    S = rankOneMatrices + positiveDefiniteComponents;
    %% Save the data
    stateEstimates.means{i} = mu;
    stateEstimates.covariances{i} = S;
    stateEstimates.cardinality(i) = targetNumber;
end
end