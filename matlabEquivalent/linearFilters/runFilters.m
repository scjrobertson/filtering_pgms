% Run the approximate Jpdaf
close all; clc;
warning off; %% Right division is iffy sometimes.
%% Load the model
model = generateModel(20, 0.75);
%% Generate the ground truth
[targetPriors, groundTruth, measurements] = generateGroundTruth(model);
%% Run the JPDAF/JIPDAF
tic; stateEstimates = jpdaf(model, targetPriors, measurements); toc;
%tic; stateEstimates = jipdafWithLoops(model, measurements); toc;
%tic; stateEstimates = jipdaf(model, measurements); toc;
%% Run the JIPDAF
%plotResults(model, groundTruth, stateEstimates, measurements, true);
resultsToTikz(model, groundTruth, stateEstimates, measurements, 'resources/pdafExample');