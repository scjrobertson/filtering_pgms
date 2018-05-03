function plotResults(model, groundTruth, stateEstimates, measurements, plotMeasurements)
% PLOTRESULTS - Plots the filter output and the ground truth
%
%   Plots the filter output, ground truth and performance metrics.
%
%   Inputs
%       groundTruth - struct. A structure with the following fields.
%           trajectories - (1, n) cell. A cell containing an array
%               for each ground truth trajectory.
%           means - (1, m) cell. Contains the Kalman filter belief
%               mean for each live trajectory at every time-step.
%           covariances - (1, m) cell. Contains the Kalman filter belief
%               covariance for each live trajectory at every time-step.
%           cardinality - (1, m) array. The number of targets present at each
%               time-step of the simulation.
%       stateEstimates - struct. A structure with the following fields.
%           means - (1, m) cell. Contains the state estimates means.
%           covariances - (1, m) cell. Contains the state estimates
%               coviariance matrices.
%           cardinality - (1, m) array. The estimated number of targets present at each
%               time-step of the simulation.
%       measurements - (1, m) cell. Contains the measurements --
%           target-generated and clutter generated measurements -- for each
%           time-step of the simulation.
%       plotMeasurements - bool. Plot graphs which include measurements.
%           This can be very slow for even a reasonable clutter rate.
%% Plot the observation space
simulationLength = size(measurements, 2);
numberOfGroundTruthTrajectories = size(groundTruth.trajectories, 2);
xLimits = model.observationSpaceLimits(1, :)';
yLimits = model.observationSpaceLimits(2, :)';
%% Calculate the error metrics
% Euclidean OSPA
eOspa = zeros(3, simulationLength);
for i = 1:simulationLength
    [eOspa(1, i), eOspa(2, i), eOspa(3, i)] = euclideanOspa(groundTruth.means{i}, stateEstimates.means{i}, model.eOspaC, model.ospaP);
end
% Hellinger OSPA
hOspa = zeros(3, simulationLength);
for i = 1:simulationLength
    [hOspa(1, i), hOspa(2, i), hOspa(3, i)] = hellingerOspa(groundTruth.means{i}, groundTruth.covariances{i}, ...
        stateEstimates.means{i}, stateEstimates.covariances{i}, model.hOspaC, model.ospaP);
end
%% Plot OSPA and cardinality vs time
figure; ospaVersusTime = gcf;
%% Plot E-OSPA vs time
subplot(311); box on; hold on; grid on;
ospaTime = model.T*(0:simulationLength-1);
% OSPA
eOspaLine = line(ospaTime, eOspa(1, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[1 0 0]);
eLocalisationLine = line(ospaTime, eOspa(2, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[0 1 0]);
eCardinalityLine = line(ospaTime, eOspa(3, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[0 0 1]);
% Limits
xlabel('t (s)'); ylabel('E-OSPA (m)');
set(gca, 'XLim', [0 model.T*(simulationLength-1)]); set(gca, 'YLim', [0 max(eOspa(1, :))+0.1]);
legend([eOspaLine, eLocalisationLine, eCardinalityLine], 'OSPA', 'Localisation', 'Cardinality');
%% Plot H-OSPA vs time
subplot(312); box on; hold on; grid on;
% OSPA
hOspaLine = line(ospaTime, hOspa(1, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[1 0 0]);
hLocalisationLine = line(ospaTime, hOspa(2, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[0 1 0]);
hCardinalityLine = line(ospaTime, hOspa(3, :), 'LineStyle','-','Marker','none','LineWidth',2,'Color',[0 0 1]);
% Limits
xlabel('t (s)'); ylabel('H-OSPA (m)');
set(gca, 'XLim', [0 model.T*(simulationLength-1)]); set(gca, 'YLim', [0 max(hOspa(1, :))+0.1]);
legend([hOspaLine, hLocalisationLine, hCardinalityLine], 'H-OSPA', 'Localisation', 'Cardinality');
%% Plot cardinality vs time
subplot(313); box on; hold on; grid on;

stairs(0:simulationLength-1, groundTruth.cardinality, '-k');
stairs(0:simulationLength-1, stateEstimates.cardinality, '-.r');
maximumCardinality = max(max(groundTruth.cardinality), max(stateEstimates.cardinality));

set(gca,'ytick',0:maximumCardinality+1);
xlabel('t(s)'); ylabel('Cardinality');
set(gca, 'XLim',[0 model.T*(simulationLength-1)]); set(gca, 'YLim',[0 (maximumCardinality + 1)]);
legend(gca,'True Cardinality','Estimated Cardinality');
%% Risk plotting measurements
if plotMeasurements == true
    %% Project the ground truth and state estimates into the observation space
    % Project ground truth
    projectedGroundTruth = cell(1, numberOfGroundTruthTrajectories);
    for i = 1:numberOfGroundTruthTrajectories
        projectedGroundTruth{i} = model.C*groundTruth.trajectories{i}(2:end, :);
    end
    % Project state estimate means
    projectedStateEstimate = cell(1, simulationLength);
    for i = 1:simulationLength
        projectedStateEstimate{i} = model.C*stateEstimates.means{i};
    end
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
        birthPoint = line(projectedGroundTruth{i}(1, 1), projectedGroundTruth{i}(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
        truthLine = line(projectedGroundTruth{i}(1, :), projectedGroundTruth{i}(2, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
        deathPoint = line(projectedGroundTruth{i}(1, end), projectedGroundTruth{i}(2, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    end
    % State Estimates
    for i = 1:simulationLength
        estimateLine = line(projectedStateEstimate{i}(1, :), projectedStateEstimate{i}(2, :), 'LineStyle','none','Marker','*','Markersize',12,'Color', 'red');
    end
    % Limits
    xlabel('x (m)'); ylabel('y (m)');
    set(gca, 'XLim', xLimits); set(gca, 'YLim', yLimits);
    legend([measurementLine, truthLine], 'Measurements', 'Ground Truth');
    %% Plot observation dimensions vs time
    figure; observationVersusTime = gcf;
    %% x-happenings vs. time
    subplot(211); box on; hold on; grid on;
    % Plot measurements
    for i = 1:simulationLength
        if ~isempty(measurements{i})
            measurementLine = line(model.T*(i-1), measurements{i}(1, :), 'LineStyle','none','Marker','+','Markersize',8,'Color',0.7*ones(1,3));
        end
    end
    % Plot ground truth
    for i = 1:numberOfGroundTruthTrajectories
        groundTruthTime = groundTruth.trajectories{i}(1, :);
        birthPoint = line(groundTruthTime(1), projectedGroundTruth{i}(1, 1), 'LineStyle','none','Marker','o', 'Markersize', 12, 'LineWidth',1,'Color',0*ones(1,3));
        truthLine = line(groundTruthTime, projectedGroundTruth{i}(1, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
        deathPoint = line(groundTruthTime(end), projectedGroundTruth{i}(1, end), 'LineStyle','none','Marker','^', 'Markersize', 12, 'LineWidth',1,'Color',0*ones(1,3));
    end
    % State Estimates
    for i = 1:simulationLength
        if ~isempty(projectedStateEstimate{i})
            estimateLine = line(model.T*(i-1), projectedStateEstimate{i}(1, :), 'LineStyle','none','Marker','*','Markersize',12,'Color', 'red');
        end
    end
    % Limits
    xlabel('t (s)'); ylabel('y (m)');
    set(gca, 'XLim', [0 model.T*(simulationLength-1)]); set(gca, 'YLim', xLimits);
    %legend([measurementLine, truthLine], 'Measurements', 'Ground Truth');
    %% y-happenings vs. time
    subplot(212); box on; hold on; grid on;
    % Plot measurements
    for i = 1:simulationLength
        if ~isempty(measurements{i})
            measurementLine = line(model.T*(i-1), measurements{i}(2, :), 'LineStyle','none','Marker','+','Markersize',8,'Color',0.7*ones(1,3));
        end
    end
    % Plot ground truth
    for i = 1:numberOfGroundTruthTrajectories
        groundTruthTime = groundTruth.trajectories{i}(1, :);
        birthPoint = line(groundTruthTime(1), projectedGroundTruth{i}(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 12, 'LineWidth',1,'Color',0*ones(1,3));
        truthLine = line(groundTruthTime, projectedGroundTruth{i}(2, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
        deathPoint = line(groundTruthTime(end), projectedGroundTruth{i}(2, end), 'LineStyle','none','Marker','^', 'Markersize', 12, 'LineWidth',1,'Color',0*ones(1,3));
    end
    % State Estimates
    for i = 1:simulationLength
        if ~isempty(projectedStateEstimate{i})
            estimateLine = line(model.T*(i-1), projectedStateEstimate{i}(2, :), 'LineStyle','none','Marker','*','Markersize',12,'Color', 'red');
        end
    end
    % Limits
    xlabel('t (s)'); ylabel('y (m)');
    set(gca, 'XLim', [0 model.T*(simulationLength-1)]); set(gca, 'YLim', xLimits);
    %legend([measurementLine, truthLine], 'Measurements', 'Ground Truth');
end
end