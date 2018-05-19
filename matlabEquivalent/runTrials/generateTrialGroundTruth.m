function [targetPriors, groundTruth, measurements] = generateTrialGroundTruth(model, simulationLength, numberOfTargets)
% GENERATEGROUNDTRUTH -- Generates the ground truth trajectories, beliefs
% and measurements.
%   [groundTruthTrajectories, groundTruthMeans, groundTruthCovariances,
%   measurements, cardinality] = generateGroundTruth(model);
%
%   Generates the measurements and the ground truth trajectories and
%   beliefs. For each trajectory, the Kalman filter is used to determine
%   the optimal belief a filter could hope to achieve.
%
%   See also generateModel.
%
%   Inputs
%       model - struct. The struct created by the function generateModel.
%
%   Outputs
%       targetPriors - struct. A struct with the following fields:
%           birthsTimes - (1, n) array. The targets' birth times.
%           deathTimes - (1, n) array. The targets' death times.
%           means - (d, n) array. The targets' prior means.
%           covariance - (d, d, n) array. The targets' prior 
%       groundTruth - struct. A structure with the following fields:
%           trajectories - (1, n) cell. A cell containing an array
%               for each ground truth trajectory.
%           rfsTrajectories - (1, n). A cell containing the ground truth as
%           an RFS.
%           means - (1, m) cell. Contains the Kalman filter belief
%               mean for each live trajectory at every time-step.
%           covariances - (1, m) cell. Contains the Kalman filter belief
%               covariance for each live trajectory at every time-step. 
%           cardinality - (1, m) array. The number of targets present at each
%               time-step of the simulation.
%       measurements - (1, m) cell. Contains the measurements --
%           target-generated and clutter generated measurements -- for each
%           time-step of the simulation.
%% Target birth and death times
targetPriors.birthTimes = ones(1, numberOfTargets);
targetPriors.deathTimes = simulationLength*ones(1, numberOfTargets);
%% Noise parameters
noiseMean = zeros(model.zDimension, 1);
noiseCovariance = (1^2)*eye(model.zDimension);
%% Target Priors
% Means
positionIndex = 1:model.numberOfSpawningLocations; %randi([1 model.numberOfSpawningLocations], [1 numberOfTargets]);
positions = model.spawnMeans(1:2, positionIndex);
velocities = ones(2, numberOfTargets);
targetPriors.means = [positions; velocities + randn([2 numberOfTargets])];
% Covariance
targetPriors.covariances = reshape(repmat(model.spawnCovariance, [1 numberOfTargets]), [model.xDimension model.xDimension numberOfTargets]);
%% Add in initial targets
targetLabels = find(targetPriors.birthTimes == 1);
x = targetPriors.means(:, targetLabels);
mu = targetPriors.means(:, targetLabels);
S = targetPriors.covariances(:, :, targetLabels);
currentTargetNumber = size(targetLabels, 2);
%% Preallocate variables
I = eye(model.xDimension);
measurements = cell(1, simulationLength);
measurements{1} = model.observationSpaceLimits(:, 1) + 2*model.observationSpaceLimits(:, 2).*rand(model.zDimension, poissrnd(model.clutterRate)); 
groundTruth.means = cell(1, simulationLength); groundTruth.means{1} = mu;
groundTruth.covariances = cell(1, simulationLength); groundTruth.covariances{1} = S;
groundTruth.cardinality = zeros(1, simulationLength); groundTruth.cardinality(1) = currentTargetNumber;
groundTruth.rfsTrajectory = cell(1, simulationLength); groundTruth.rfsTrajectory{1} = mu;
groundTruth.trajectories{numberOfTargets} = []; 
for i = targetLabels; groundTruth.trajectories{i} = [0; x(:, i)]; end
%% Determine the ground truth
for i = 2:simulationLength
    time = model.T*(i-1);
    %% Add in measurements
    missedDetection = rand([1 currentTargetNumber]) > model.detectionProbability;
    numberOfTargetGeneratedMeasurements = sum(~missedDetection);
    numberOfClutterReturns = poissrnd(model.clutterRate);
    measurements{i} = zeros(model.zDimension, numberOfTargetGeneratedMeasurements + numberOfClutterReturns);
    measurements{i}(:, 1:numberOfClutterReturns) = model.observationSpaceLimits(:, 1) + 2*model.observationSpaceLimits(:, 2).*rand(model.zDimension, numberOfClutterReturns); 
    %% Predict the targets' states
    x = model.A*x; % Exact ground truth
    muPred = model.A*mu; 
    SPred = quadraticMultiprod(S, model.A, model.Atranspose, model.xDimension, currentTargetNumber) + model.R; 
    %% Create the update components
    targetGeneratedMeasurements = model.C*x + mvnrnd(noiseMean, noiseCovariance, currentTargetNumber)'; % Exact value of ground truth
    muZ = model.C*muPred; % Predicted measurement mean
    Kxz = rightMultiprod(SPred, model.Ctranspose, model.xDimension, model.xDimension, model.zDimension, currentTargetNumber); % Cross covariance
    Z = leftMultiprod(model.C, Kxz, model.zDimension, model.xDimension, model.zDimension, currentTargetNumber) + model.Q; % Predicted measurement covariance
    ZInv = simplifiedMultinv(Z, model.zDimension, currentTargetNumber); % Inverse of Z, used to determine likelihoods
    K = simplifiedMultiprod(Kxz, ZInv, model.xDimension, model.zDimension, model.zDimension, currentTargetNumber); % Kalman gain
    STemp = rightMultiprod(K, model.C, model.xDimension, model.zDimension, model.xDimension, currentTargetNumber);
    SUpdated = simplifiedMultiprod((I - STemp), SPred, model.xDimension, model.xDimension, model.xDimension, currentTargetNumber); % Update covariance matrices
    %% Update the targets' states
    difference = reshape(targetGeneratedMeasurements - muZ, [1 model.zDimension currentTargetNumber]);
    muPredReshaped = reshape(muPred, [model.xDimension 1 currentTargetNumber]);
    muUpdated = reshape(muPredReshaped + sum(bsxfun(@times, K, difference), 2), [model.xDimension currentTargetNumber]);
    %% Decide whether a target missed a detection
    mu = muUpdated; mu(:, missedDetection) = muPred(:, missedDetection);
    S = SUpdated; S(:, :, missedDetection) = SPred(:, :, missedDetection);
    measurements{i}(:, numberOfClutterReturns+1:end) = targetGeneratedMeasurements(:, ~missedDetection);
    %% Add in new targets
    newTargetLabels = find(targetPriors.birthTimes == i);
    numberOfNewTargets = size(newTargetLabels, 2);
    if numberOfNewTargets ~= 0
        newIndices = (currentTargetNumber+1):(currentTargetNumber + numberOfNewTargets);
        targetLabels(newIndices) = newTargetLabels;
        x(:, newIndices) = targetPriors.means(:, newTargetLabels);
        mu(:, newIndices) = targetPriors.means(:, newTargetLabels);
        S(:, :, newIndices) = targetPriors.covariances(:, :, newTargetLabels);
        currentTargetNumber = currentTargetNumber + numberOfNewTargets;
    end
    %% Remove dead targets
    aliveTargetIndices = targetPriors.deathTimes(targetLabels) >= i; 
    if currentTargetNumber ~= 0
        targetLabels = targetLabels(aliveTargetIndices);
        x = x(:, aliveTargetIndices);
        mu = mu(:, aliveTargetIndices);
        S = S(:, :, aliveTargetIndices);
        currentTargetNumber = size(targetLabels, 2);
    end
    %% Update ground truth trajectories
    groundTruth.rfsTrajectory{i} = x;
    for j = 1:currentTargetNumber; groundTruth.trajectories{targetLabels(j)} = [groundTruth.trajectories{targetLabels(j)} [time; x(:, j)]]; end
    %% Update ground truth beliefs
    groundTruth.means{i} = mu;
    groundTruth.covariances{i} = S;
    groundTruth.cardinality(i) = currentTargetNumber;
end
end