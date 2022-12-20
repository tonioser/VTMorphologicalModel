function [scores, cntsRes, errRes] = regress_Data_BasisVector_2_Scores(cnts, basisVector, meanCnts)
% 
% Linear regression of data on a basis vector to estimate the control parameter/score
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Data contours
%                                    Typically of size 1 x 108 x 1
%     basisVector(nbPts,nbDim)     : Basis vector
%                                    Typically of size 108 x 1
%     meanCnts(nbPts,nbDim)        : Mean contour of the articulations
%                                    Typically of size 108 x 1
% 
% Outputs
%     scores(nbSubjects)              : Column vector of the scores
%     cntsRes(nbSubjects,nbPts,nbDim) : Residual contours
%     errRes(1)                       : Residual error
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Sizes
nbObs = size(cnts,1);
nbPts = size(cnts,2);
nbDim = size(cnts,3);

% Centred data points
cntsC = cnts - permute(repmat(meanCnts,[1,1,nbObs]),[3,1,2]);

% Indices of the non-NaN points in the data
indAnA = find(~isnan(cnts(1,:,1)));
nbPtsAnA = length(indAnA);

% Non-NaN data points
ptsAnAC = cntsC(:,indAnA,:);
basisVectorAnA = basisVector(indAnA,:);

% Reshape the data in a 2D matrix
matPtsAnAC = reshape(ptsAnAC, nbObs, nbPtsAnA*nbDim);

% Reshape the basisVector
matBasisVectorAnA = reshape(basisVectorAnA, 1, nbPtsAnA*nbDim);

% Linear regression on the basis vector
scores = (matBasisVectorAnA' \ matPtsAnAC')';

% Prediction of the centred data with the obtained scores
matPtsAnAC_Est = scores * matBasisVectorAnA;

% Residue
matResAnA = matPtsAnAC - matPtsAnAC_Est;

% Residual error
errRes = norm(matResAnA(:));

% Reshape for the output
cntsRes = NaN(size(cnts));
cntsRes(:,indAnA,:) = reshape(matResAnA, [nbObs, nbPtsAnA, nbDim]);

end