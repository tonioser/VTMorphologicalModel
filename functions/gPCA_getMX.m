function [MX, MXC] = gPCA_getMX(cnts, UTy, indPhaVT)
% 
% Calculate the MX control parameter for a set of articulation contours.
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Articulation contours
%                                    Typically of size 41 x 1692 x 2
%     UTy(1)                       : Y value of the Upper Teeth point (UT)
%                                    Typically equal to 10
%     indPhaVT(nbPtsVT)            : Indices of the vocal tract points corresponding to the pharynx for an articulation contour
%                                    Typically of length 200
% 
% Outputs
%     MX(nbSubjects)  : Column vector of the control parameter MX
%     MXC(nbSubjects) : Column vector of the centred control parameter MX
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Index of the closest pharynx point to Y = UTy on the average contour
cnt = squeeze(mean(cnts,1));
iPtPhaVT = argmin(abs(cnt(indPhaVT,2) - UTy));
iMX = indPhaVT(iPtPhaVT);

% Morphological parameter
MX = cnts(:,iMX,1);
MXC = MX - mean(MX);

end