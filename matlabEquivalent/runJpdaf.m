% Run the approximate Jpdaf
%% Load the model
model = generateModel(20, 0.95);
%% Generate the ground truth
[groundTruthTrajectories, groundTruthMeans, groundTruthCovariances, cardinality] = generateGroundTruth(model);