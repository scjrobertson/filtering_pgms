function [groundTruth, C, measurements, stateEstimates, ospa, cardinality] = readIniFiles(formatName)
%% Load simulation information
groundTruthFileName = [formatName, 'GroundTruth.ini'];
xDimension = inifile(groundTruthFileName, 'read', {'SIMULATION INFO', '', 'xDimension'}); xDimension = str2num(xDimension{1});
zDimension = inifile(groundTruthFileName, 'read', {'SIMULATION INFO', '', 'zDimension'}); zDimension = str2num(zDimension{1});
simulationLength = inifile(groundTruthFileName, 'read', {'SIMULATION INFO', '', 'simulationLength'}); simulationLength = str2num(simulationLength{1});
numberOfTrajectories = inifile(groundTruthFileName, 'read', {'SIMULATION INFO', '', 'numberOfGroundTruthTrajectories'}); numberOfTrajectories = str2num(numberOfTrajectories{1});
observationModel = inifile(groundTruthFileName, 'read', {'SIMULATION INFO', '', 'observationModel'}); observationModel = str2num(observationModel{1});
C = [observationModel(1:xDimension); observationModel(xDimension+1:end)];
%% Read and process ground truth
groundTruth = cell(1, numberOfTrajectories);
for i = 1:numberOfTrajectories
    trajectoryName = sprintf('TRAJECTORY %d', i);
    result = inifile(groundTruthFileName, 'read', {trajectoryName, '', 'trajectory'});
    result = str2num(result{1});
    trajectoryLength = size(result, 2)/(xDimension+1);
    groundTruth{i} = reshape(result, [xDimension+1 trajectoryLength]);
end
%% Read and process measurements
measurementFileName = [formatName, 'Measurements.ini'];
measurements = cell(1, simulationLength);

for i = 1:simulationLength
    measurementName = sprintf('MEASUREMENTS %d', i);
    result = inifile(measurementFileName, 'read', {measurementName, '', 'measurments'});
    result = str2num(result{1});
    numberOfMeasurements = size(result, 2)/(zDimension);
    measurements{i} = reshape(result, [zDimension numberOfMeasurements]);
end
%% Read and process state estimates
stateEstimateFileName = [formatName, 'StateEstimates.ini'];
stateEstimates = cell(1, simulationLength);

for i = 1:simulationLength
    stateEstimateNames = sprintf('STATE ESTIMATES %d', i);
    targetNumber = inifile(stateEstimateFileName, 'read', {stateEstimateNames, '', 'numberOfTargets'}); targetNumber = str2num(targetNumber{1});
    stateEstimates{i} = cell(1, targetNumber);
    for j = 1:targetNumber
       targetName = sprintf('target%d', j);
       result = inifile(stateEstimateFileName, 'read', {stateEstimateNames, '', targetName});
       result = str2num(result{1});
       numberOfMixtureComponents = size(result, 2)/(xDimension)
       stateEstimates{i}{j} = reshape(result, [xDimension numberOfMixtureComponents]);
    end
end
%% Read and process OSPA
ospaResultFileName = [formatName, 'OspaResults.ini'];

ospaResult = inifile(ospaResultFileName, 'read', {'OSPA', '', 'ospa'});
ospaResult = str2num(ospaResult{1});
ospaLength = size(ospaResult, 2)/(3);
ospa = reshape(ospaResult, [3 ospaLength]);

cardinalityResult = inifile(ospaResultFileName, 'read', {'OSPA', '', 'cardinality'});
cardinality = str2num(cardinalityResult{1});
end