addpath data;
close all; clc;
%% Load the data
[groundTruth, C, measurements, stateEstimates, ospa, cardinality] = readIniFiles('simple');
colours = {'red', 'blue'};
%% Admin
simulationLength = size(measurements, 2);
numberOfGroundTruthTrajectories = size(groundTruth, 2);
xLimits = [-100 100];
yLimits = [-100 100];
cardinalityEstimate = zeros(1, simulationLength);
%% Plot the observation space
figure; fullObservationSpace = gcf; box on; hold on; grid on;
% Plot measurements
for i = 1:simulationLength
    if ~isempty(measurements{i})
       measurementLine = line(measurements{i}(1, :), measurements{i}(2, :), 'LineStyle','none','Marker','+','Markersize',8,'Color',0.7*ones(1,3)); 
    end
end
% Plot ground truth
for i = 1:numberOfGroundTruthTrajectories
    projectedGroundTruth = C*groundTruth{i}(2:end, :);
    
    birthPoint = line(projectedGroundTruth(1, 1), projectedGroundTruth(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 8,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(projectedGroundTruth(1, :), projectedGroundTruth(2, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',0*ones(1,3));
    deathPoint = line(projectedGroundTruth(1, end), projectedGroundTruth(2, end), 'LineStyle','none','Marker','^', 'Markersize', 8,'LineWidth',1,'Color',0*ones(1,3));
end
% State Estimates
for i = 1:simulationLength
    numberOfTargets = size(stateEstimates{i}, 2);
    cardinalityEstimate(i) = numberOfTargets;
    for j = 1:numberOfTargets
        projectedStateEstimate = C*stateEstimates{i}{j};
        estimateLine = line(projectedStateEstimate(1, :), projectedStateEstimate(2, :), 'LineStyle','none','Marker','*','Markersize',12,'Color', colours{j});
    end
end
% Limits
xlabel('x (m)'); ylabel('y (m)');
set(gca, 'XLim', xLimits); set(gca, 'YLim', yLimits);
legend([measurementLine, truthLine, estimateLine], 'Measurements', 'Ground Truth', 'State Estimates');
%% Plot observation dimensions vs time
figure; observationVersusTime = gcf; 
%% x-happenings vs. time
subplot(211); box on; hold on; grid on;
% Plot measurements
for i = 1:simulationLength
    if ~isempty(measurements{i})
       measurementLine = line(i-1, measurements{i}(1, :), 'LineStyle','none','Marker','+','Markersize',8,'Color',0.7*ones(1,3)); 
    end
end
% Plot ground truth
for i = 1:numberOfGroundTruthTrajectories
    groundTruthTime = groundTruth{i}(1, :);
    projectedGroundTruth = C*groundTruth{i}(2:end, :);
    
    birthPoint = line(groundTruthTime(1), projectedGroundTruth(1, 1), 'LineStyle','none','Marker','o', 'Markersize', 8, 'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(groundTruthTime, projectedGroundTruth(1, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',0*ones(1,3));
    deathPoint = line(groundTruthTime(end), projectedGroundTruth(1, end), 'LineStyle','none','Marker','^', 'Markersize', 8, 'LineWidth',1,'Color',0*ones(1,3));
end
% State Estimates
for i = 1:simulationLength
    numberOfTargets = size(stateEstimates{i}, 2);
    for j = 1:numberOfTargets
        projectedStateEstimate = C*stateEstimates{i}{j};
        estimateLine = line(i-1, projectedStateEstimate(1, :), 'LineStyle','none','Marker','*','Markersize',12,'Color', colours{j});
    end
end
% Limits
xlabel('t (s)'); ylabel('y (m)');
set(gca, 'XLim', [0 simulationLength-1]); set(gca, 'YLim', xLimits);
%legend([measurementLine, truthLine], 'Measurements', 'Ground Truth');
%% y-happenings vs. time
subplot(212); box on; hold on; grid on;
% Plot measurements
for i = 1:simulationLength
    if ~isempty(measurements{i})
       measurementLine = line(i-1, measurements{i}(2, :), 'LineStyle','none','Marker','+','Markersize',8,'Color',0.7*ones(1,3)); 
    end
end
% Plot ground truth
for i = 1:numberOfGroundTruthTrajectories
    groundTruthTime = groundTruth{i}(1, :);
    projectedGroundTruth = C*groundTruth{i}(2:end, :);
    
    birthPoint = line(groundTruthTime(1), projectedGroundTruth(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 8, 'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(groundTruthTime, projectedGroundTruth(2, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',0*ones(1,3));
    deathPoint = line(groundTruthTime(end), projectedGroundTruth(2, end), 'LineStyle','none','Marker','^', 'Markersize', 8, 'LineWidth',1,'Color',0*ones(1,3));
end
% State Estimates
for i = 1:simulationLength
    numberOfTargets = size(stateEstimates{i}, 2);
    for j = 1:numberOfTargets
        projectedStateEstimate = C*stateEstimates{i}{j};
        estimateLine = line(i-1, projectedStateEstimate(2, :), 'LineStyle','none','Marker','*','Markersize',12,'Color', colours{j});
    end
end
% Limits
xlabel('t (s)'); ylabel('y (m)');
set(gca, 'XLim', [0 simulationLength-1]); set(gca, 'YLim', xLimits);
%legend([measurementLine, truthLine], 'Measurements', 'Ground Truth');
%% Plot OSPA and cardinality vs time
figure; ospaVersusTime = gcf; 
%% Plot OSPA vs time
subplot(211); box on; hold on; grid on;
% OSPA
ospaLine = line(1:simulationLength, ospa(1, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[1 0 0]);
localisationLine = line(1:simulationLength, ospa(2, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[0 1 0]);
cardinalityLine = line(1:simulationLength, ospa(3, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[0 0 1]);
% Limits
xlabel('t (s)'); ylabel('OSPA (m)');
set(gca, 'XLim', [0 simulationLength-1]); set(gca, 'YLim', [0 1.25]);
legend([ospaLine, localisationLine, cardinalityLine], 'OSPA', 'Localisation', 'Cardinality');
%% Plot cardinality vs time
subplot(212); box on; hold on; grid on;

stairs(0:simulationLength-1, cardinality, '-k');
stairs(0:simulationLength-1, cardinalityEstimate, '-.r');
maximumCardinality = max(max(cardinality), max(cardinalityEstimate));

set(gca,'ytick',0:maximumCardinality+1);
xlabel('t(s)'); ylabel('Cardinality');
set(gca, 'XLim',[0 simulationLength-1]); set(gca, 'YLim',[0 (maximumCardinality + 1)]);
legend(gca,'True Cardinality','Estimated Cardinality');