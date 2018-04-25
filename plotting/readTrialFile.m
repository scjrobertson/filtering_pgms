function [ospa] = readTrialFile(fileName)
%% Load simulation information
information = csvread(fileName, 0, 0, [0 0 3 0]);
numberOfTrials = information(1);
numberOfSimulationsPerTrial = information(2);
simulationLength = information(3)+1;
ospaC = information(4);
%% Load raw data
rawOspa = csvread(fileName, 4, 0);
ospa = reshape(rawOspa, [numberOfTrials, numberOfSimulationsPerTrial, simulationLength]);
ospa = ospa(:, :, 1:end-1);
end