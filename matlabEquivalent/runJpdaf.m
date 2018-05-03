% Run the approximate Jpdaf
close all; clc;
addpath common statisticsToolBox;
%% Load the model
model = generateModel(2, 0.95);
%% Generate the ground truth
[targetPriors, groundTruth, measurements] = generateGroundTruth(model);
%% Run the filter
[stateEstimates, cardinalityEstimates] = runFilter(model, targetPriors, measurements);
%% Plot results
%plotResults(model, groundTruth, measurements, cardinalityEstimate)