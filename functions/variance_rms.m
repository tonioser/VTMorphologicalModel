function [varexTot, RMSTot, varex, RMS, varPts] = variance_rms(pts, ptsEst);
% 
% Percentage of variance explanation and RMS reconstruction error between
% reference data and data estimated by means of a (guided) PCA model
% 
% Inputs
%     pts(nbSubjects,nbPts,nbDim) : Reference contour points
%                                   Typically of size 41 x 1036 x 2
%     pts(nbSubjects,nbPts,nbDim) : Estimated contour points, typically by means of a PCA model
%                                   Typically of size 41 x 1036 x 2
% 
% Outputs
%     varexTot(1)         : Overall percentage of variance explanation (between 0 and 1)
%     RMSTot(1)           : Overall RMS reconstruction error
%     varex(nbPts,nbDim)  : Percentage of variance explanation, per point (between 0 and 1)
%     RMS(nbPts,nbDim)    : RMS reconstruction error, per point
%     varPts(nbPts,nbDim) : Variance of the reference data
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Sizes
[nbObs, nbPts, nbDdim] = size(pts);

% Reshape
matPts = reshape(pts, nbObs, nbPts*nbDdim);
matPtsEst = reshape(ptsEst, nbObs, nbPts*nbDdim);

% Data variance, per point
matVar = var(matPts,1);

% Percentage of variance explanation, per point
% = variance of the estimated data in reference to the variance of the data
matVarex = var(matPtsEst,1) ./ matVar;

% RMS error, per point
matRMS = std(matPtsEst - matPts, 1);

% Overall percentage of variance explanation
varexTot = sum(var(matPtsEst,1)) / sum(matVar);

% Overal RMS error
RMSTot = sqrt(mean(var(matPtsEst - matPts, 1)));

% Reshape
varPts = squeeze(reshape(matVar, 1, nbPts, nbDdim));
varex = squeeze(reshape(matVarex, 1, nbPts, nbDdim));
RMS = squeeze(reshape(matRMS, 1, nbPts, nbDdim));

end