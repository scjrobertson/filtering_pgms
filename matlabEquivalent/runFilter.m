% Run the approximate Jpdaf
close all; clc;
addpath common statisticsToolBox;
warning off; %% Right division is iffy sometimes.
%% Load the model
model = generateModel(1, 0.95);
%% Generate the ground truth
[targetPriors, groundTruth, measurements] = generateGroundTruth(model);
%% Run the filter
%tic; stateEstimates = jpdaf(model, targetPriors, measurements); toc;
tic; stateEstimates = jipdaf(model, measurements); toc;
%% Plot results
%plotResults(model, groundTruth, stateEstimates, measurements, true);