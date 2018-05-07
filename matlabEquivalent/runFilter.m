% Run the approximate Jpdaf
close all; clc;
addpath common statisticsToolBox;
warning off; %% Right division is iffy sometimes.
%% Load the model
model = generateModel(4, 0.95);
%% Generate the ground truth
%[targetPriors, groundTruth, measurements] = generateGroundTruth(model);
%% Run the TOMB
%tic; stateEstimates = jpdaf(model, targetPriors, measurements); toc;
tic; stateEstimatesTomb = exactTomb(model, measurements); toc;
%tic; stateEstimates = jipdaf(model, measurements); toc;
%% Run the JIPDAF
plotResults(model, groundTruth, stateEstimatesTomb, measurements, true);
%plotResults(model, groundTruth, stateEstimates, measurements, true);