% Run the approximate Jpdaf
close all; clc;
addpath ../common ../statisticsToolBox;
warning off; %% Right division is iffy sometimes.
%% Load the model
model = generateNonLinearModel(20, 0.95);
%% Generate the ground truth
%[targetPriors, groundTruth, measurements] = generateNonLinearGroundTruth(model);
%% Run the JPDAF/JIPDAF
tic; stateEstimates = ukJipdaf(model, measurements); toc;
%% Run the JIPDAF
plotResults(model, groundTruth, stateEstimates, measurements, false);