function [groundTruthTrajectories, groundTruthMeans, groundTruthCovariances, cardinality] = generateGroundTruth(model)
%% Simulation length
simulationLength = 50;
%% Target birth and death times
numberOfTargets = 2;
birthTimes = [1 10];
deathTimes = [50 50];
%% Target Priors
% Means
targetPriorMean = zeros(model.xDimension, numberOfTargets);
targetPriorMean(1, :) = [-40; 40; 1; -2];
targetPriorMean(2, :) = [-40; -40; 2; 2];
% Covariance
targetPriorCovariance = zeros(model.xDimension, model.xDimension, numberOfTargets);
targetPriorCovariance(1, :)  = (3^2)*eye(model.xDimension);
targetPriorCovariance(2, :)  = (3^2)*eye(model.xDimension);
%% Generate the ground truth
groundTruthTrajectories = cell(1, numberOfTargets);
groundTruthMeans = cell(1, simulationLength);
groundTruthCovariances = cell(1, simulationLength);
cardinality = zeros(1, simulationLength);


% Add in initial targets
for i = 1:numberOfTargets
   if (bi 
end

for i = 1:simulationLength
    
    
    
    
    
end

end