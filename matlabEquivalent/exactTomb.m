function stateEstimates = exactTomb(model, measurements)
% EXACTTOMB -- Estimates the target's states using approximate JPDAF
%   stateEstimates = exactTOMB(model, targetPriors, measurements)
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
%% State estimates
stateEstimates.means = cell(1, simulationLength);
stateEstimates.covariances = cell(1, simulationLength);
stateEstimates.cardinality = zeros(1, simulationLength);
%% Declare state variables
r = zeros(0, 0); mu = zeros(model.xDimension, 0); S = zeros(model.xDimension, model.xDimension, 0); % Established targets' states
lambdaU = zeros(0); xU = zeros(model.xDimension, 0); SU = zeros(model.xDimension, model.xDimension, 0); % Poisson point process parameters
targetNumber = 0;
%% Run the filter
for i = 2:3%simulationLength
    %% Predict the established targets' states
    r = model.survivalProbability*r;
    mu = model.A*mu + model.u;
    for j = 1:targetNumber
        S(:, :, j) = model.A*S(:, :, j)*model.Atranspose + model.R;
    end
    %% Predict the state of the PPP
    lambdaU = model.poissonSurvivalProbability*lambdaU;
    survivingPoissonIndices = lambdaU > model.lambdaThreshold;
    lambdaU = lambdaU(survivingPoissonIndices);
    xU = model.A*xU(:, survivingPoissonIndices) + model.u;
    SU = SU(:, :, survivingPoissonIndices);
    intensitySize = size(lambdaU, 2);
    for j = 1:intensitySize
        SU(:, :, j) = model.A*SU(:, :, j)*model.Atranspose + model.R;
    end
    %% Add in newly spawned targets
    intensitySize = intensitySize + model.numberOfSpawningLocations;
    lambdaU(end+1:intensitySize) = model.newTargetProbability;
    xU(:, end+1:intensitySize) = model.spawnMeans;
    SU(:, :, end+1:intensitySize) = model.spawnCovariances;
    %% Create update components for the established targets
    numberOfMeasurements = size(measurements{i}, 2);
    zMeas = measurements{i};
    associationMatrix = zeros(targetNumber, numberOfMeasurements+1);
    rUpdated = zeros(targetNumber, numberOfMeasurements+1);
    muUpdated = zeros(model.xDimension, targetNumber, numberOfMeasurements+1);
    SUpdated = zeros(model.xDimension, model.xDimension, targetNumber, numberOfMeasurements+1);
    for j = 1:targetNumber
        associationMatrix(j, 1) = 1 - r(j) + r(j)*(1-model.detectionProbability);
        rUpdated(j, 1) = r(j)*(1-model.detectionProbability)/associationMatrix(j, 1);
        muUpdated(:, j, 1) = mu(:, j);
        SUpdated(:, :, j, 1) = S(:, :, j);
        
        Z = model.C*S(:, :, j)*model.Ctranspose + model.Q;
        detNorm = sqrt(det(2*pi*Z));
        K = S(:, :, j)*model.Ctranspose/Z;
        SFinal = S(:, :, j) - K*model.C*S(:, :, j);
        for k = 1:numberOfMeasurements
            v =  zMeas(:, k) - model.C*mu(:, j);
            associationMatrix(j, k+1) = r(j)*model.detectionProbability*exp(-0.5*(v'/Z)*v)/detNorm;
            rUpdated(j, k+1) = 1;
            muUpdated(:, j, k+1) = mu(:, j) + K*v;
            SUpdated(:, :, j, k+1) = SFinal;
        end
    end
    %% Create tentative track for each measurement
    intensityLikelihood = zeros(1, numberOfMeasurements);
    rNew = zeros(numberOfMeasurements, 1);
    xNew = zeros(model.xDimension, numberOfMeasurements);
    SNew = zeros(model.xDimension, model.xDimension, numberOfMeasurements);
    ZU = zeros(model.zDimension, model.zDimension, intensitySize);
    KU = zeros(model.xDimension, model.zDimension, intensitySize);
    SUTemp = zeros(model.xDimension, model.xDimension, intensitySize);
    detZU = zeros(intensitySize, 1);
    marginalLikelihoods = zeros(intensitySize, 1);
    muU = zeros(model.xDimension, intensitySize);
    % Precompute components
    for j = 1:intensitySize
        ZU(:, :, j) = model.C*SU(:, :, j)*model.Ctranspose + model.Q;
        detZU(j) = sqrt(det(2*pi*ZU(:, :, j)));
        KU(:, :, j) = SU(:, :, j)*model.Ctranspose/ZU(:, :, j);
        SUTemp(:, :, j) = SU(:, :, j) - KU(:, :, j)*model.C*SU(:, :, j);
    end
    % Create new tracks
    for j = 1:numberOfMeasurements
        for k = 1:intensitySize
            v = zMeas(:, j) - model.C*xU(:, k);
            marginalLikelihoods(k) = lambdaU(k)*model.detectionProbability*exp(-0.5*(v'/ZU(:, :, k))*v)/detZU(k);
            muU(:, k) = xU(:, k) + KU(:, :, k)*v;
        end
        totalMass = sum(marginalLikelihoods);
        intensityLikelihood(j) = totalMass + model.clutterPerUnitVolume;
        rNew(j) = totalMass/intensityLikelihood(j);
        marginalLikelihoods = marginalLikelihoods/totalMass;
        marginalLikelihoods(isnan(marginalLikelihoods)) = 0;
        xNew(:, j) = muU*marginalLikelihoods;
        for k = 1:intensitySize
           v = xNew(:, j) - muU(:, k);
           SNew(:, :, j) = SNew(:, :, j) + marginalLikelihoods(k)*(SUTemp(:, :, k) + v*v');
        end
    end
    %% Thin the PPP
    lambdaU = (1-model.detectionProbability)*lambdaU;
    survivingPoissonIndices = lambdaU > model.lambdaThreshold;
    lambdaU = lambdaU(survivingPoissonIndices);
    xU = model.A*xU(:, survivingPoissonIndices) + model.u;
    SU = SU(:, :, survivingPoissonIndices);
    %% Loopy belief propagation
    [updatedAssociationMatrix, updatedIntensityLikelihoods] = loopyBeliefPropagation(associationMatrix, intensityLikelihood, 1e-6, 200);
    %% Update the established targets' states
    numberOfTentativeTracks = targetNumber + numberOfMeasurements;
    r = zeros(numberOfTentativeTracks, 1);
    mu = zeros(model.xDimension, numberOfTentativeTracks);
    S = zeros(model.xDimension, model.xDimension, numberOfTentativeTracks);
    
    for j = 1:targetNumber
       rProbability = updatedAssociationMatrix(j, :).*rUpdated(j, :); 
       r(j) = sum(rProbability);
       rProbability = rProbability'/r(j);
       mu(:, j) = squeeze(muUpdated(:, j, :))*rProbability;
       for k = 1:(numberOfMeasurements+1)
           v = mu(:, j) - muUpdated(:, j, k);
           S(:, :, j) = S(:, :, j) + rProbability(k)*(SUpdated(:, :, j, k) + v*v'); 
       end
    end
    %% Add in new tenetative tracks
    r(targetNumber+1:numberOfTentativeTracks) = updatedIntensityLikelihoods.*intensityLikelihood;
    mu(:, targetNumber+1:numberOfTentativeTracks) = xNew;
    S(:, :, targetNumber+1:numberOfTentativeTracks) = SNew;
    %% Gate the existing tracks
    significantIndices = r > model.existenceThreshold;
    r = r(significantIndices);
    mu = mu(:, significantIndices);
    S = S(:, :, significantIndices);
    targetNumber = size(r, 1);
    %% Save the data
    stateEstimates.means{i} = mu;
    stateEstimates.covariances{i} = S;
    stateEstimates.cardinality(i) = targetNumber;
end
end