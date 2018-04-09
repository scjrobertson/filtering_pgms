addpath drawTools/;
addpath uHellingerGMM/;
%% Admin
numberOfTrials = 10;
numberOfSimulationsPerTrial = 50000;

lambda = 1;
error = zeros(numberOfTrials, numberOfSimulationsPerTrial);
%% Run the simulation
for n = 1:numberOfTrials
   for j = 1:numberOfSimulationsPerTrial 
        % Gaussian 1
        mu1 = rand(n, 1);
        S1 = generateCovarianceMatrix(n, 2);        
        % Gaussian 2
        mu2 = rand(n, 1);
        S2 = generateCovarianceMatrix(n, 2);
        % GMM for uHellinger
        gmm1.w = [1]; gmm1.Mu = mu1; gmm1.Cov{1} = S1;
        gmm2.w = [1]; gmm2.Mu = mu2; gmm2.Cov{1} = S2;
        % Distance
        exactDistance = hellingerDistance(mu1, mu2, S1, S2);
        approxDistance = uHellingerJointSupport2_ND(gmm1, gmm2);
        % Calculate error
        error(n, j) = abs(exactDistance-approxDistance);
   end
end
%% Calculate statistics
meanError = mean(error, 2);
stdDevError = std(error, 0, 2);
%% Plot
xAxis = 1:numberOfTrials;

figure(); hold on; grid on; box on;
bar(xAxis, meanError);
errorbar(xAxis, meanError, stdDevError, '.');
set(gca,'xtick', xAxis);
xlabel('\lambda_{c}'); ylabel('Error');


function A = generateCovarianceMatrix(n, eigenMean)
Q = randn(n,n);
A = Q'*diag(abs(eigenMean+randn(n,1)))*Q;
end