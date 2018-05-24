close all; clc;
warning off; %% Right division is iffy sometimes.
%% Trial properties
numberOfTrials = 50;
numberOfSimulationsPerTrial = 5000;
simulationLength = 250;
numberOfTargets = 2;
eOspaMean = zeros(numberOfSimulationsPerTrial, numberOfTrials);
hOspaMean = zeros(numberOfSimulationsPerTrial, numberOfTrials);
%% Parameter
variableParameter = linspace(5, 250, numberOfTrials);
%% Run the trials
for i = 1:numberOfTrials
    fprintf('Trial %d\n', i);
    %% Generate the model
    model = generateModel(variableParameter(i), 0.75);
    %% Iterate through the simulations
    eOspa = zeros(simulationLength, numberOfSimulationsPerTrial);
    hOspa = zeros(simulationLength, numberOfSimulationsPerTrial);
    for j = 1:numberOfSimulationsPerTrial
        % Load the ground truth
        [targetPriors, groundTruth, measurements] = generateTrialGroundTruth(model, simulationLength, numberOfTargets);
        % Run the filter
        stateEstimates = jpdaf(model, targetPriors, measurements);
        % Calculate the Euclidean and Hellinger OSPA
        for k = 1:simulationLength
            eOspa(k, j) = euclideanOspa(groundTruth.rfsTrajectory{k}, stateEstimates.means{k}, model.eOspaC, model.ospaP);
            hOspaTemp = ospaSpecific(groundTruth.means{k}, groundTruth.covariances{k}, stateEstimates.means{k}, stateEstimates.covariances{k}, model.hOspaC, model.ospaP);
            hOspa(k, j) = hOspaTemp(1, :);
        end
        % Calculate mean OSPA for the simulation
        eOspaMean(j, i) = mean(eOspa(:, j), 1);
        hOspaMean(j, i) = mean(hOspa(:, j), 1);
    end
    save('resources/cache/simulationResults.mat', 'eOspaMean', 'hOspa');
end
%% Calculate the mean Euclidean and Hellinger OSPA
meanEOspa = mean(eOspaMean, 1);
stdDevEOspa = std(eOspaMean, 0, 1);
meanHOspa = mean(hOspaMean, 1);
stdDevHOspa = std(hOspaMean, 0, 1);
%% Plot a bar chart for the E-OSPA
indices = find(meanEOspa-stdDevEOspa < 0);
stdDevEOspa(indices) = meanEOspa(indices);

figure(); hold on; grid on; box on;
bar(variableParameter, meanEOspa);
errorbar(variableParameter, meanEOspa, stdDevEOspa, 'r.');
set(gca,'xtick', variableParameter);
xlabel('\lambda_{c}'); ylabel('E-OSPA (m)');
matlab2tikz('resources/multiplePdaf/multiplePdafClutterTrialsEOspa.tex');
%% Plot a bar chart for the E-OSPA
indices = find(meanHOspa-stdDevHOspa < 0);
stdDevHOspa(indices) = meanHOspa(indices);

figure(); hold on; grid on; box on;
bar(variableParameter, meanHOspa);
errorbar(variableParameter, meanHOspa, stdDevHOspa, 'r.');
set(gca,'xtick', variableParameter);
xlabel('\lambda_{c}'); ylabel('H-OSPA (m)');
matlab2tikz('resources/multiplePdaf/multiplePdafClutterTrialsHOspa.tex');