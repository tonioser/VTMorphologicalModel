function [MPC, MPCC, coefsPalMPC, meanPalMPC] = gPCA_getMPC(cnts, indPal)
% 
% Calculate the MPC control parameter for a set of articulation contours.
%
% Returns also the coeffcicients and mean of the raw PCA used for
% calculating the control parameter.
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Articulation contours
%                                    Typically of size 41 x 1692 x 2
%     indPal(nbPtsPal)             : Indices of the vocal tract points corresponding to the hard palate for an articulation contour
%                                    Typically of length 108
% 
% Outputs
%     MPC(nbSubjects)                 : Column vector of the control parameter MPC
%     MPCC(nbSubjects)                : Column vector of the centred control parameter MPC (=MPC)
%     coefsPalMPC(1, nbPtsPal, nbDim) : Coefficients of the raw PCA performed on the coordinates of the palate points
%     meanPalMPC(nbPtsPal, nbDim)     : Mean of the coordinates of the palate points over the subjects 
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Sizes
nbObs = size(cnts,1);
nbPtsPal = length(indPal);
nbDim = size(cnts,3);

% PCA on the hard palate contours (from which the anteriority has been removed)
ptsPal = cnts(:,indPal,:);
matPtsPal = reshape(ptsPal, nbObs, nbPtsPal*nbDim);
[coeff, MPCC] = pca(matPtsPal, 'NumComponents', 1);

% Morphological parameter
MPC = MPCC;
MPCC;

% Return also the palate PCA model
coefsPalMPC = reshape(coeff, [1, nbPtsPal, nbDim]);
meanPalMPC = squeeze(mean(ptsPal));

end