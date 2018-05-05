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
simulationLength = 250;
%% Target birth and death times
numberOfTargets = 200;
%targetPriors.birthTimes = ones(1, numberOfTargets);  
%targetPriors.deathTimes = simulationLength*ones(1, numberOfTargets);
targetPriors.birthTimes = sort(randi([1 simulationLength], [1 numberOfTargets]));
targetPriors.deathTimes = randi([1 simulationLength], [1 numberOfTargets]) + targetPriors.birthTimes;
targetPriors.deathTimes(targetPriors.deathTimes > simulationLength) = simulationLength;
%% Noise parameters
noiseMean = zeros(model.zDimension, 1);
noiseCovariance = (1^2)*eye(model.zDimension);
%% Target Priors
% Means
positions = model.observationSpaceLimits(:, 1) + 2*model.observationSpaceLimits(:, 2).*rand(model.xDimension/2, numberOfTargets);
velocities = randi([-2 2], [2 numberOfTargets]) + randn([2 numberOfTargets]);
targetPriors.means = [positions; velocities];
% Covariance
covarianceMatrix = (2^2)*eye(model.xDimension);
targetPriors.covariances = reshape(repmat(covarianceMatrix, [1 numberOfTargets]), [model.xDimension model.xDimension numberOfTargets]);
%% Preallocate variables
measurements = cell(1, simulationLength);
groundTruth.means = cell(1, simulationLength);
groundTruth.covariances = cell(1, simulationLength);
groundTruth.trajectories = cell(1, numberOfTargets);
groundTruth.cardinality = zeros(1, simulationLength);
%% Add in clutter and initialise beliefs
for i = 1:simulationLength
   measurements{i} = model.observationSpaceLimits(:, 1) + ...
       2*model.observationSpaceLimits(:, 2).*rand(model.zDimension, poissrnd(model.clutterRate)); 
   groundTruth.cardinality(i) = sum( (i >= targetPriors.birthTimes).*(i <= targetPriors.deathTimes), 2);
end
%% Determine the ground truth trajectories and beliefs
for i = 1:numberOfTargets
    lifeSpan = targetPriors.deathTimes(i) - (targetPriors.birthTimes(i)-1);
    groundTruth.trajectories{i} = zeros(model.xDimension+1, lifeSpan); 
    groundTruth.trajectories{i}(1, :) = model.T*(targetPriors.birthTimes(i)-1:targetPriors.deathTimes(i)-1);
    groundTruth.trajectories{i}(2:end, 1) = targetPriors.means(:, i);
    %% Initialise the track
    absoluteTime = targetPriors.birthTimes(i);
    mu = targetPriors.means(:, i);
    S = targetPriors.covariances(:, :, i);
    groundTruth.means{absoluteTime}(:, end+1) = mu;
    groundTruth.covariances{absoluteTime}(:, :, size(groundTruth.means{absoluteTime}, 2)) = S;
    %% Determine the Kalman filter belief for trajectory lifespan
    for j = 2:lifeSpan
        absoluteTime = absoluteTime + 1;
        groundTruth.trajectories{i}(2:end, j) = model.A*groundTruth.trajectories{i}(2:end, j-1) + model.u;
        z = model.C*groundTruth.trajectories{i}(2:end, j) + mvnrnd(noiseMean, noiseCovariance, 1)';
        %% Predicted state of the target
        muPred = model.A*mu + model.u;
        SPred = model.A*S*model.Atranspose + model.R;
        %% Measurement update
        if (rand > model.detectionProbability) % If a detection is missed
            mu = muPred;
            S = SPred;
        else
            %% Kalman filter update
            K = SPred*(model.Ctranspose)/(model.C*SPred*model.Ctranspose + model.Q);
            mu = muPred + K*(z - model.C*muPred);
            S = SPred - K*model.C*SPred;
            %% Add measurement to list
            measurements{absoluteTime}(:, end+1) = z; 
        end
        %% Append to ground truth beliefs
        groundTruth.means{absoluteTime}(:, end+1) = mu;
        groundTruth.covariances{absoluteTime}(:, :, size(groundTruth.means{absoluteTime}, 2)) = S;
    end
end
end