function [ospa] = readTrialFile(fileName)

trialFile = [fileName];
numberOfTrials = inifile(trialFile, 'read', {'SIMULATION INFORMATION', '', 'numberOfTrials'}); numberOfTrials = str2num(numberOfTrials{1});
numberOfSimulationsPerTrial = inifile(trialFile, 'read', {'SIMULATION INFORMATION', '', 'numberOfSimulationsPerTrial'}); numberOfSimulationsPerTrial = str2num(numberOfSimulationsPerTrial{1});
simulationLength = inifile(trialFile, 'read', {'SIMULATION INFORMATION', '', 'simulationLength'}); simulationLength = str2num(simulationLength{1});
ospaC = inifile(trialFile, 'read', {'SIMULATION INFORMATION', '', 'ospaC'}); ospaC = str2num(ospaC{1});

ospa = zeros(numberOfTrials, numberOfSimulationsPerTrial, simulationLength);

for i = 1:numberOfTrials
    trialName = sprintf('TRIAL %d', i);
    for j = 1:numberOfSimulationsPerTrial
        simulationName = sprintf('simulation%d', j);
        ospaResult = inifile(trialFile, 'read', {trialName, '', simulationName})
        ospa(i, j, :) = str2num(ospaResult{1});
    end
end

end