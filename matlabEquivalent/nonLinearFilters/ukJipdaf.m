function stateEstimates = ukJipdaf(model, measurements)
% UKJIPDAF -- Estimates the target's states using approximate JPDAF
%   stateEstimates = ukJipdaf(model, targetPriors, measurements)
%
%   Runs an approximate JIPDAF - it resolves the data association by way of
%   loopy belief propagation. This implementation avoid for loops and is
%   not very readable. Using matrices instead of 3D arrays throughout would
%   probably be more efficient.
%
%   See also generateModel, generateGroundTruth, jpdaf, jipdaf, loopyBeliefPropagation
%   and plotResults.
%
%   Inputs
%       model - struct. The struct declared in generateModel, it has fields
%           describing the linear motion and measurement models, etc.
%       measurements - (1, m) cell. Contains the measurements --
%           target-generated and clutter generated measurements -- for each
%           time-step of the simulation.
%   Output
%       stateEstimates - struct. A structure with the following fields.
%           means - (1, m) cell. Contains the state estimates means.
%           covariances - (1, m) cell. Contains the state estimates
%               covariance matrices.
%           cardinality - (1, m) array. The estimated number of targets present at each
%               time-step of the simulation.
%% Admin
simulationLength = size(measurements, 2);
I = eye(model.xDimension);
%% State estimates
stateEstimates.means = cell(1, simulationLength);
stateEstimates.covariances = cell(1, simulationLength);
stateEstimates.cardinality = zeros(1, simulationLength);
%% Declare some constants
ONE = ones(model.xDimension, model.xDimension);
numberOfSigmaPoints = 2*model.xDimension + 1;
%% Declare state variables
r = zeros(0, 0); mu = zeros(model.xDimension, 0); S = zeros(model.xDimension, model.xDimension, 0); % Established targets' states
lambdaU = zeros(0); xU = zeros(model.xDimension, 0); SU = zeros(model.xDimension, model.xDimension, 0); % Poisson point process parameters
targetNumber = 0;
%% Run the filter
for i = 1:simulationLength
    %% Thin, then predict the Poisson Point Process
    lambdaU = model.poissonSurvivalProbability*lambdaU;
    survivingPoissonIndices = lambdaU > model.lambdaThreshold;
    lambdaU = lambdaU(survivingPoissonIndices);
    intensitySize = size(lambdaU, 2);
    xU = model.A*xU(:, survivingPoissonIndices) + model.u;
    SU = quadraticMultiprod(SU(:, :, survivingPoissonIndices), model.A, model.Atranspose, model.xDimension, intensitySize) + model.R;
    %% Add tentative spawned targets to the PPP
    intensitySize = intensitySize + model.numberOfSpawningLocations;
    lambdaU(end+1:intensitySize) = model.newTargetProbability;
    xU(:, end+1:intensitySize) = model.spawnMeans;
    SU(:, :, end+1:intensitySize) = model.spawnCovariances;
    %% Predict the targets' states
    if (targetNumber ~= 0)
        rPred = model.survivalProbability*r;
        muPred = model.A*mu + model.u;
        SPred =  quadraticMultiprod(S, model.A, model.Atranspose, model.xDimension, targetNumber) + model.R;
    end
    %% Start update
    for j = 1:model.numberOfSensors
        numberOfMeasurements = size(measurements{i}(j).z, 2);
        if (targetNumber ~= 0)
            %% Do Cholesky decompostion
            blockMatrix = kron(speye(targetNumber), ONE);
            blockIndices = logical(blockMatrix);
            blockMatrix(blockIndices)= model.utGamma*SPred(:);
            lowerCholesky = chol(blockMatrix, 'lower');
            lowerCholesky = full(lowerCholesky(blockIndices));
            %% Determine the targets' sigma points
            XUShaped = reshape(muPred, [model.xDimension 1 1 targetNumber]);
            choleskyColumns = reshape(lowerCholesky, [model.xDimension model.xDimension 1 targetNumber]);
            sigmaPointsPacked = zeros(model.xDimension, 1, numberOfSigmaPoints, targetNumber);
            sigmaPointsPacked(:, :, 2:end, :) = [choleskyColumns -choleskyColumns];
            sigmaPoints = reshape(sigmaPointsPacked + XUShaped, [model.xDimension numberOfSigmaPoints*targetNumber]);
            %% Put the sigma points through the transform
            relativePosition = sigmaPoints(1:2, :) - model.sensorPosition(:, j);
            range = sqrt(relativePosition(1, :).^2 + relativePosition(2, :).^2);
            radialVelocity = (1./range).*sum(real(relativePosition).*sigmaPoints(3:4,:));
            transformedSigmaPoints = reshape([range; radialVelocity], [model.zDimension 1 numberOfSigmaPoints targetNumber]);
            %% Calculate weights
            meanWeights = reshape(repmat([model.meanWeight model.dimensionWeights], [1 1 targetNumber]), [1 1 numberOfSigmaPoints targetNumber]);
            covWeights = reshape(repmat([model.covarianceWeight model.dimensionWeights], [1 1 targetNumber]), [1 1 numberOfSigmaPoints targetNumber]);
            %% Caluculate the targets' means
            weightedPoints = sum(bsxfun(@times, meanWeights, transformedSigmaPoints), 2);
            muZ = sum(weightedPoints, 3);
            %% Calculate the covariance
            shiftPoints = transformedSigmaPoints - muZ;
            shiftPointsPermuted = permute(shiftPoints, [2 1 3 4]);
            covarianceOuterProduct = bsxfun(@times, shiftPoints, shiftPointsPermuted);
            Z = sum(bsxfun(@times, covWeights, covarianceOuterProduct), 3) + model.Q;
            %% Calculate the cross-covariance
            crossOuterProduct = bsxfun(@times, sigmaPointsPacked, shiftPointsPermuted);
            Kxz = sum(bsxfun(@times, covWeights, crossOuterProduct), 3);
            %% Reshape the mean, calculate inverse and determinants
            muZ = reshape(muZ, [model.zDimension targetNumber]);
            ZInv = simplifiedMultinv(Z, model.zDimension, targetNumber);
            detZ = reshape(Z(1, 1, :).*Z(2, 2, :) - Z(1, 2, :).*Z(2, 1, :), [targetNumber 1]);
            normalisingConstants = model.detectionProbability./sqrt( ((2*pi)^model.zDimension)*detZ);
            K = simplifiedMultiprod(Kxz, ZInv, model.xDimension, model.zDimension, model.zDimension, targetNumber);
            STemp = simplifiedMultiprod(Kxz, permute(K, [2 1 3]), model.xDimension, model.zDimension, model.xDimension, targetNumber);
            SUpdated = SPred - STemp;
            %% Create updated state components
            associationMatrix = repmat(1 - model.detectionProbability*rPred, [1 numberOfMeasurements+1]);
            rUpdated = ones(targetNumber, numberOfMeasurements+1); rUpdated(:, 1) = (1-model.detectionProbability)*rPred./associationMatrix(:, 1);
            muZRep = repmat(muZ, [1 numberOfMeasurements]);
            zRep = reshape(repmat(measurements{i}(j).z, [targetNumber 1]),[model.zDimension targetNumber*numberOfMeasurements]);
            %% Calculate the likelihoods
            difference = reshape(zRep - muZRep, [1 model.zDimension targetNumber numberOfMeasurements]);
            differencePermuted = permute(difference, [2 1 3 4]);
            ZInvShaped = repmat(ZInv, [1 1 1 numberOfMeasurements]);
            leftProduct = sum(bsxfun(@times, ZInvShaped, difference), 2);
            rightProduct = reshape(sum(bsxfun(@times, differencePermuted, leftProduct), 1), [targetNumber numberOfMeasurements]);
            associationMatrix(:, 2:end) = repmat(rPred, [1 numberOfMeasurements]).*repmat(normalisingConstants, [1 numberOfMeasurements]).*exp(-0.5*rightProduct);
        end
        %% Do Cholesky decompostion for the PPP
        blockMatrix = kron(speye(intensitySize), ONE);
        blockIndices = logical(blockMatrix);
        blockMatrix(blockIndices)= model.utGamma*SU(:);
        lowerCholesky = chol(blockMatrix, 'lower');
        lowerCholesky = full(lowerCholesky(blockIndices));
        %% Determine the targets' sigma points
        XUShaped = reshape(xU, [model.xDimension 1 1 intensitySize]);
        choleskyColumns = reshape(lowerCholesky, [model.xDimension model.xDimension 1 intensitySize]);
        sigmaPointsPacked = zeros(model.xDimension, 1, numberOfSigmaPoints, intensitySize);
        sigmaPointsPacked(:, :, 2:end, :) = [choleskyColumns -choleskyColumns];
        sigmaPoints = reshape(sigmaPointsPacked + XUShaped, [model.xDimension numberOfSigmaPoints*intensitySize]);
        %% Put the sigma points through the transform
        relativePosition = sigmaPoints(1:2, :) - model.sensorPosition(:, j);
        range = sqrt(relativePosition(1, :).^2 + relativePosition(2, :).^2);
        radialVelocity = (1./range).*sum(real(relativePosition).*sigmaPoints(3:4,:));
        transformedSigmaPoints = reshape([range; radialVelocity], [model.zDimension 1 numberOfSigmaPoints intensitySize]);
        %% Calculate weights
        meanWeights = reshape(repmat([model.meanWeight model.dimensionWeights], [1 1 intensitySize]), [1 1 numberOfSigmaPoints intensitySize]);
        covWeights = reshape(repmat([model.covarianceWeight model.dimensionWeights], [1 1 intensitySize]), [1 1 numberOfSigmaPoints intensitySize]);
        %% Caluculate the targets' means
        weightedPoints = sum(bsxfun(@times, meanWeights, transformedSigmaPoints), 2);
        muZU = sum(weightedPoints, 3);
        %% Calculate the covariance
        shiftPoints = transformedSigmaPoints - muZU;
        shiftPointsPermuted = permute(shiftPoints, [2 1 3 4]);
        covarianceOuterProduct = bsxfun(@times, shiftPoints, shiftPointsPermuted);
        ZU = sum(bsxfun(@times, covWeights, covarianceOuterProduct), 3) + model.Q;
        %% Calculate the cross-covariance
        crossOuterProduct = bsxfun(@times, sigmaPointsPacked, shiftPointsPermuted);
        KUxz = sum(bsxfun(@times, covWeights, crossOuterProduct), 3);    
        %% Create update components for the PPP
        ZUInv = simplifiedMultinv(ZU, model.zDimension, intensitySize); % Inverse of Z, used to determine likelihoods
        detZU = reshape(ZU(1, 1, :).*ZU(2, 2, :) - ZU(2, 1, :).*ZU(1, 2, :), [intensitySize 1]);
        normalisingConstants = (lambdaU').*(model.detectionProbability)./sqrt((2*pi)^(model.zDimension)*detZU);
        KU = simplifiedMultiprod(KUxz, ZUInv, model.xDimension, model.zDimension, model.zDimension, intensitySize); % Kalman gain
        SUTemp = simplifiedMultiprod(KUxz, permute(KU, [2 1 3]), model.xDimension, model.zDimension, model.xDimension, intensitySize);
        SUUpdated = SU - SUTemp; % Update covariance matrices
        %% Create update state components
        muZURep = repmat(reshape(muZU, [model.zDimension intensitySize]), [1 numberOfMeasurements]);
        zRep = reshape(repmat(measurements{i}(j).z,[intensitySize 1]),[model.zDimension intensitySize*numberOfMeasurements]);
        %% Determine likelihoods
        intensityDifference = reshape(zRep - muZURep, [1 model.zDimension intensitySize numberOfMeasurements]);
        intensityDifferencePermuted = permute(intensityDifference, [2 1 3 4]);
        ZUInvShaped = repmat(ZUInv, [1 1 1 numberOfMeasurements]);
        leftProduct = sum(bsxfun(@times, ZUInvShaped, intensityDifference), 2);
        rightProduct = reshape(sum(bsxfun(@times, intensityDifferencePermuted, leftProduct), 1), [intensitySize numberOfMeasurements]);
        clutterLikelihoodMatrix = repmat(normalisingConstants, [1 numberOfMeasurements]).*exp(-0.5*rightProduct);
        totalMass = sum(clutterLikelihoodMatrix, 1);
        totalIntensityLikelihoods = totalMass + model.clutterPerUnitVolume;
        rNew = totalMass./totalIntensityLikelihoods;
        %% Update the PPP states
        % Determine means
        intensityLikelihoods =  clutterLikelihoodMatrix./totalMass;
        intensityLikelihoods(isnan(intensityLikelihoods)) = 0;
        intensityProbabilities = reshape(intensityLikelihoods, [1 1 intensitySize numberOfMeasurements]);
        xUReshaped = reshape(xU, [model.xDimension 1 intensitySize]);
        xURepeated = repmat(xUReshaped, [1 1 1 numberOfMeasurements]);
        xUUpdated = xURepeated + sum(bsxfun(@times, KU, permute(intensityDifferencePermuted, [2 1 3 4])), 2);
        if intensitySize == 1; xNew = sum(bsxfun(@times, squeeze(xUUpdated), intensityProbabilities), 2);
        else; xNew = sum(bsxfun(@times, xUUpdated, intensityProbabilities), 3); end
        % Determine covariance matrices
        xUShift = xUUpdated - xNew;
        outerProduct = bsxfun(@times, xUShift, permute(xUShift, [2 1 3 4]));
        xNew = reshape(xNew, [model.xDimension numberOfMeasurements]);
        SUReshaped = repmat(SUUpdated, [1 1 1 numberOfMeasurements]);
        SNew = squeeze(sum(bsxfun(@times, intensityProbabilities, SUReshaped), 3)) + squeeze(sum(bsxfun(@times, intensityProbabilities, outerProduct), 3));
        %% Thin the PPP
        lambdaU = (1-model.detectionProbability)*lambdaU;
        significantPoissonIndices = lambdaU > model.lambdaThreshold;
        lambdaU = lambdaU(significantPoissonIndices);
        xU = xU(:, significantPoissonIndices);
        SU = SU(:, :, significantPoissonIndices);
        intensitySize = size(lambdaU, 2);
        %% Resolve data association by way of loopy belief propagation
        if (targetNumber ~= 0)
            [updatedAssociationMatrix, updatedIntensityLikelihoods] = loopyBeliefPropagation(associationMatrix, totalIntensityLikelihoods, 1e-6, 200);
            %% Update states and compute weak marginals
            % Existence probability
            rProbability = rUpdated.*updatedAssociationMatrix;
            r = sum(rProbability, 2);
            rProbability = rProbability./r;
            % Means
            associationProbabilties =  permute(reshape(rProbability, [1 targetNumber numberOfMeasurements+1]), [1 3 2]);
            missedDetectionProbabilities = associationProbabilties(1, 1, :);
            muUpdated = zeros([model.xDimension 1 targetNumber numberOfMeasurements+1]);
            muUpdated(:, :, :, 1) = reshape(muPred, [model.xDimension 1 targetNumber]);
            muUpdated(:, :, :, 2:end) = muUpdated(:, :, :, 1) + sum(bsxfun(@times, K, difference), 2);
            if targetNumber == 1; mu = sum(bsxfun(@times, squeeze(muUpdated), associationProbabilties), 2);
            else; mu = sum(bsxfun(@times, permute(squeeze(muUpdated), [1 3 2]), associationProbabilties), 2); end
            % Covariance matrices
            muShift = permute(muUpdated - mu, [1 2 4 3]);
            outerProduct = bsxfun(@times, muShift, permute(muShift, [2 1 3 4]));
            associationProbabiltiesReshaped = reshape(associationProbabilties, [1 1 numberOfMeasurements+1 targetNumber]);
            rankOneMatrices = squeeze(sum(bsxfun(@times, associationProbabiltiesReshaped, outerProduct), 3));
            positiveDefiniteComponents = bsxfun(@times, missedDetectionProbabilities, SPred) + bsxfun(@times, 1- missedDetectionProbabilities, SUpdated);
            % Reshape and reassign
            mu = reshape(mu, [model.xDimension targetNumber]);
            S = rankOneMatrices + positiveDefiniteComponents;
        else
            updatedIntensityLikelihoods = totalIntensityLikelihoods./sum(totalIntensityLikelihoods, 2);
        end
        %% Add in new tentative tracks
        targetNumber = targetNumber + numberOfMeasurements;
        r(end+1:targetNumber, :) = (updatedIntensityLikelihoods.*rNew)';
        mu(:, end+1:targetNumber) = xNew;
        S(:, :, end+1:targetNumber) = SNew;
        %% Gate the tracks
        significantIndices = r > model.existenceThreshold;
        r = r(significantIndices);
        mu = mu(:, significantIndices);
        S = S(:, :, significantIndices);
        targetNumber = size(r, 1);
        %% Reallocate
        rPred = r;
        muPred = mu;
        SPred = S;
    end
    %% Save the data
    stateEstimates.means{i} = mu;
    stateEstimates.covariances{i} = S;
    stateEstimates.cardinality(i) = targetNumber;
end
end