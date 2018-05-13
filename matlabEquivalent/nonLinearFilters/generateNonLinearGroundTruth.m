function [targetPriors, groundTruth, measurements] = generateGroundTruth(model)
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
%% Simulation length
simulationLength = 1000;
%% Target birth and death times
numberOfTargets = 5;
targetPriors.birthTimes = ones(1, 5); %sort(randi([1 simulationLength], [1 numberOfTargets]));
targetPriors.deathTimes = randi([1 simulationLength], [1 numberOfTargets]) + targetPriors.birthTimes;
targetPriors.deathTimes(targetPriors.deathTimes > simulationLength) = simulationLength;
%% Noise parameters
noiseMean = zeros(model.zDimension, 1);
noiseCovariance = (1^2)*eye(model.zDimension);
%% Target Priors
% Means
positionIndex = randi([1 model.numberOfSpawningLocations], [1 numberOfTargets]);
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
ONE = ones(model.xDimension, model.xDimension);
numberOfSigmaPoints = 2*model.xDimension + 1;
measurements = cell(1, simulationLength);
groundTruth.means = cell(1, simulationLength); groundTruth.means{1} = mu;
groundTruth.covariances = cell(1, simulationLength); groundTruth.covariances{1} = S;
groundTruth.cardinality = zeros(1, simulationLength); groundTruth.cardinality(1) = currentTargetNumber;
groundTruth.rfsTrajectory = cell(1, simulationLength); groundTruth.rfsTrajectory{1} = mu;
groundTruth.trajectories{numberOfTargets} = [];
for i = targetLabels; groundTruth.trajectories{i} = [0; x(:, i)]; end
%% Add in initial measurements
for i = 1:model.numberOfSensors
    numberOfClutterReturns = poissrnd(model.clutterRate);
    rangeComponent = model.maximumSensorRange*rand(1, numberOfClutterReturns);
    dopplerComponent = -model.maximumDopplerVelocity + 2*model.maximumDopplerVelocity*rand(1, numberOfClutterReturns);
    measurements{1}(i).z = [rangeComponent; dopplerComponent];
end
%% Determine the ground truth
for i = 2:simulationLength
    time = model.T*(i-1);
    %% Predict the targets' states
    x = model.A*x; % Exact ground truth
    muPred = model.A*mu;
    SPred = quadraticMultiprod(S, model.A, model.Atranspose, model.xDimension, currentTargetNumber) + model.R;
    for j = 1:model.numberOfSensors
        %% Add in measurements
        missedDetection = rand([1 currentTargetNumber]) > model.detectionProbability;
        numberOfTargetGeneratedMeasurements = sum(~missedDetection);
        numberOfClutterReturns = poissrnd(model.clutterRate);
        rangeComponent = model.maximumSensorRange*rand(1, numberOfClutterReturns);
        dopplerComponent = -model.maximumDopplerVelocity + 2*model.maximumDopplerVelocity*rand(1, numberOfClutterReturns);
        measurements{i}(j).z = zeros(model.zDimension, numberOfTargetGeneratedMeasurements + numberOfClutterReturns);
        measurements{i}(j).z(:, 1:numberOfClutterReturns) = [rangeComponent; dopplerComponent];
        %% Determine target-generated measurements
        groundTruthRelativePosition = x(1:2, :) - model.sensorPosition(:, j);
        groundTruthRange = sqrt(groundTruthRelativePosition(1, :).^2 + groundTruthRelativePosition(2, :).^2);
        groundTruthRadialVelocity = (1./groundTruthRange).*sum(real(groundTruthRelativePosition).*x(3:4,:));
        targetGeneratedMeasurements = [groundTruthRange; groundTruthRadialVelocity] + mvnrnd(noiseMean, noiseCovariance, currentTargetNumber)';
        %% Do Cholesky decompostion
        blockMatrix = kron(speye(currentTargetNumber), ONE);
        blockIndices = logical(blockMatrix);
        blockMatrix(blockIndices)= model.utGamma*SPred(:);
        lowerCholesky = chol(blockMatrix, 'lower');
        lowerCholesky = full(lowerCholesky(blockIndices));
        %% Determine the targets' sigma points
        muShaped = reshape(muPred, [model.xDimension 1 1 currentTargetNumber]);
        choleskyColumns = reshape(lowerCholesky, [model.xDimension model.xDimension 1 currentTargetNumber]);
        sigmaPointsPacked = zeros(model.xDimension, 1, numberOfSigmaPoints, currentTargetNumber);
        sigmaPointsPacked(:, :, 2:end, :) = [choleskyColumns -choleskyColumns];
        sigmaPoints = reshape(sigmaPointsPacked + muShaped, [model.xDimension numberOfSigmaPoints*currentTargetNumber]);
        %% Put the sigma points through the transform
        relativePosition = sigmaPoints(1:2, :) - model.sensorPosition(:, j);
        range = sqrt(relativePosition(1, :).^2 + relativePosition(2, :).^2);
        radialVelocity = (1./range).*sum(real(relativePosition).*sigmaPoints(3:4,:));
        transformedSigmaPoints = reshape([range; radialVelocity], [model.zDimension 1 numberOfSigmaPoints currentTargetNumber]);
        %% Calculate weights
        meanWeights = reshape(repmat([model.meanWeight model.dimensionWeights], [1 1 currentTargetNumber]), [1 1 numberOfSigmaPoints currentTargetNumber]);
        covWeights = reshape(repmat([model.covarianceWeight model.dimensionWeights], [1 1 currentTargetNumber]), [1 1 numberOfSigmaPoints currentTargetNumber]);
        %% Calculate the mean
        weightedPoints = sum(bsxfun(@times, meanWeights, transformedSigmaPoints), 2);
        muZ = sum(weightedPoints, 3);
        %% Calculate the covariance
        shiftPoints = transformedSigmaPoints - muZ;
        shiftPointsPermuted = permute(shiftPoints, [2 1 3 4]);
        covarianceOuterProduct = bsxfun(@times, shiftPoints, shiftPointsPermuted);
        Z = sum(bsxfun(@times, covWeights, covarianceOuterProduct), 3) + model.Q;
        %% Calculate the cross-covariance
        crossOuterProduct = bsxfun(@times, sigmaPointsPacked, shiftPointsPermuted);
        Kxz = sum(bsxfun(@times, covWeights, crossOuterProduct), 3);
        %% Reshape the mean, calculate inverse and determinants
        muZ = reshape(muZ, [model.zDimension currentTargetNumber]);
        ZInv = simplifiedMultinv(Z, model.zDimension, currentTargetNumber);
        K = simplifiedMultiprod(Kxz, ZInv, model.xDimension, model.zDimension, model.zDimension, currentTargetNumber);
        STemp = simplifiedMultiprod(Kxz, permute(K, [2 1 3]), model.xDimension, model.zDimension, model.xDimension, currentTargetNumber);
        SUpdated = SPred - STemp;
        %% Update the mean
        difference = reshape(targetGeneratedMeasurements - muZ, [1 model.zDimension currentTargetNumber]);
        muPredReshaped = reshape(muPred, [model.xDimension 1 currentTargetNumber]);
        muUpdated = reshape(muPredReshaped + sum(bsxfun(@times, K, difference), 2), [model.xDimension currentTargetNumber]);
        %% Decide whether a target missed a detection
        mu = muUpdated; mu(:, missedDetection) = muPred(:, missedDetection);
        S = SUpdated; S(:, :, missedDetection) = SPred(:, :, missedDetection);
        measurements{i}(j).z(:, numberOfClutterReturns+1:end) = targetGeneratedMeasurements(:, ~missedDetection);
        %% Reallocate for next iteration
        muPred = mu;
        SPred = S;
    end
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