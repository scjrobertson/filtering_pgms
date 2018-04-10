addpath data;
close all; clc;
%% Load the data
fileName = 'detectionTrials.csv';
ospa = readTrialFile(fileName);
[numberOfTrials, numberOfSimulationsPerTrial, simulationLength] = size(ospa);
%% Calculate mean and standard deviation for each trial
meanOspa = zeros(1, numberOfTrials);
stdDevOspa = zeros(1, numberOfTrials);

for i = 1:numberOfTrials
    trialMean = zeros(1, numberOfSimulationsPerTrial);
    trialStd = zeros(1, numberOfSimulationsPerTrial);
    for j = 1:numberOfSimulationsPerTrial
        trialMean(j) = mean(ospa(i, j, :));
        trialStd(j) = std(ospa(i, j, :));
    end
    meanOspa(i) = mean(trialMean);
    stdDevOspa(i) = mean(trialStd);
end
%% Plot a bar chart
indices = find(meanOspa-stdDevOspa < 0);
stdDevOspa(indices) = meanOspa(indices);

xAxis = 0.05*(0:numberOfTrials-1);

figure(); hold on; grid on; box on;
bar(xAxis, meanOspa);
errorbar(xAxis, meanOspa, stdDevOspa, 'r.');
set(gca,'xtick', xAxis);
xlabel('\lambda_{c}'); ylabel('OSPA (m)');