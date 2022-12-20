function cnts = predict_Scores_BasisVectors_2_Data(scores, basisVectors, meanCnts)
% 
% Prediction of data from scores and basis vectors.
% 
% It assumes the following model:
%        cntsEst = scores * basisVector + meanCnts
% 
% Inputs
%     scores(nbSubjects,nbComp)       : Column vector of the centred scores
%                                       Typically of size 41 x 5
%     basisVector(nbComp,nbPts,nbDim) : Basis vectors
%                                       Typically of size 41 x 1692 x 2
%     meanCnts(nbPts,nbDim)           : Mean contour of the articulations
%                                       Typically of size 1692 x 2
% 
% Outputs
%     cnts(nbSubjects,nbPts,nbDim) : Prediction of the articulation contours
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Sizes
nbObs = size(scores,1);
nbComp = size(scores,2);
nbPts = size(meanCnts,1);
nbDim = size(meanCnts,2);

% Reshape
matBasisVectors = reshape(basisVectors, nbComp, nbPts*nbDim);
matMeanCnts = reshape(meanCnts, 1, nbPts*nbDim);

% Prediction of the centred data
matCntsC = scores * matBasisVectors;

% Add the mean
matCnts = matCntsC + repmat(matMeanCnts, nbObs, 1);

% Reshape
cnts = reshape(matCnts, [nbObs, nbPts, nbDim]);

end