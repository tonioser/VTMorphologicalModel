function [basisVector, meanCnts, cntsRes, errRes] = regress_Data_Scores_2_BasisVector(cnts, scores)
% 
% Multiple linear regressions to get a basis vector from data and scores.
% 
% It build a model where the points can be estimated thanks as:
%        cntsEst = scores * basisVector + meanCnts
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Articulation contours
%                                    Typically of size 41 x 1036 x 2
%     scores(nbSubjects)           : Column vector of the centred scores
% 
% Outputs
%     basisVector(nbPts,nbDim)        : Basis vector
%     meanCnts(nbPts,nbDim)           : Mean contour of the articulations
%     cntsRes(nbSubjects,nbPts,nbDim) : Residual contours
%     errRes(1)                       : Residual error
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Sizes
nbObs = size(cnts,1);
nbPts = size(cnts,2);
nbDim = size(cnts,3);

% Mean observation
meanCnts = squeeze(mean(cnts));

% Centred data
cntsC = cnts - permute(repmat(meanCnts,[1,1,nbObs]),[3,1,2]);

% Indices of the non-NaN points in the data
indAnA = find(~isnan(cnts(1,:,1)));
nbPtsAnA = length(indAnA);

% Non-NaN data points
ptsAnAC = cntsC(:,indAnA,:);

% Reshape the data in a 2D matrix
matPtsAnAC = reshape(ptsAnAC, nbObs, nbPtsAnA*nbDim);

% Restriction of the data to the linearly independant observations to avoid a
% rank deficiency in the linear regression
[~, indObsIndep]=licols(matPtsAnAC'); % Linearly independant observations
matPtsAnACIndep = matPtsAnAC(indObsIndep,:); % Restriction of the data to the linearly independant observations 
scoresIndep = scores(indObsIndep); % % Restriction of the scores to the linearly independant observations

% Multiple linear regressions
matBasisVector = scoresIndep \ matPtsAnACIndep;

% Prediction of the centred data with the obtained basis vector
matPtsAnAC_Est = scores * matBasisVector;

% Residue
matResAnA = matPtsAnAC - matPtsAnAC_Est;

% Residual error
errRes = norm(matResAnA(:));

% Reshape for the output
basisVector = NaN(nbPts,nbDim);
basisVector(indAnA,:) = reshape(matBasisVector, nbPtsAnA, nbDim);
cntsRes = NaN(size(cnts));
cntsRes(:,indAnA,:) = reshape(matResAnA, [nbObs, nbPtsAnA, nbDim]);

end