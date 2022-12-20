function [MPA, MPAC, coefsPalMPA, meanPalMPA] = gPCA_getMPA(cnts, indPal)
% 
% Calculate the MPA control parameter for a set of articulation contours.
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
%     MPA(nbSubjects)          : Column vector of the control parameter MPA
%     MPAC(nbSubjects)         : Column vector of the centred control parameter MPA (=MPA)
%     coefsPalMPA(1, nbPtsPal) : Coefficients of the raw PCA performed on the X-values of the palate points
%     meanPalMPA(nbPtsPal, 1)  : Mean of the X-values of the palate points over the subjects 
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% PCA on the X-coordinates of the hard palate contours
ptsPal = cnts(:,indPal,1);
[coeff, MPAC] = pca(ptsPal, 'NumComponents', 1);

% Morphological parameter
MPA = MPAC;
MPAC;

% Return also the palate PCA model
coefsPalMPA = reshape(coeff, [1, length(coeff), 1]);
meanPalMPA = mean(ptsPal)';

end