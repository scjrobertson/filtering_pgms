% Run the approximate Jpdaf
close all; clc;
addpath common statisticsToolBox;
warning off; %% Right division is iffy - sometimes
%% Load the model
model = generateModel(6, 0.98);
%% Generate the ground truth
[targetPriors, groundTruth, measurements] = generateGroundTruth(model);
%% Run the filter
tic;
stateEstimates = runFilter(model, targetPriors, measurements);
toc;
%% Plot results
plotResults(model, groundTruth, stateEstimates, measurements, false)