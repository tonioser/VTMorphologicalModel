function [MY, MYC] = gPCA_getMY(cnts, iGF, iGB)
% 
% Calculate the MY control parameter for a set of articulation contours.
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Articulation contours
%                                    Typically of size 41 x 1692 x 2
%     iGF(1)                       : Index of the point corresponding to the anterior of the glottis for an articulation contour
%                                    Typically of value 1631
%     iGB(1)                       : Index of the point corresponding to the posterior of the glottis for an articulation contour
%                                    Typically of value 1650
% 
% Outputs
%     MY(nbSubjects)  : Column vector of the control parameter MY
%     MYC(nbSubjects) : Column vector of the centred control parameter MY
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Middle of the anterior and posterior glottis points (larynx height)
LH = (cnts(:,iGF,2) + cnts(:,iGB,2)) / 2;

% Morphological parameter
MY = LH;
MYC = MY - mean(MY);

end