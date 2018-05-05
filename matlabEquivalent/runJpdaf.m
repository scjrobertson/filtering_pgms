% Run the approximate Jpdaf
close all; clc;
addpath common statisticsToolBox;
warning off; %% Right division is iffy sometimes.
%% Load the model
model = generateModel(2, 0.95);
%% Generate the ground truth
[targetPriors, groundTruth, measurements] = generateGroundTruth(model);
%% Run the filter
profile on;
tic; stateEstimates = runFilter(model, targetPriors, measurements); toc;
profile off;
profile viewer;
%% Plot results
plotResults(model, groundTruth, stateEstimates, measurements, false)