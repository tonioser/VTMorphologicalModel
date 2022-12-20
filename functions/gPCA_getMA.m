function [MA, MAC] = gPCA_getMA(cnts, iPhL, iPhU)
% 
% Calculate the MA control parameter for a set of articulation contours.
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim) : Articulation contours
%                                    Typically of size 41 x 1692 x 2
%     iPhL(1)                      : Index of the lower point of the pharynx for an articulation contour
%                                    Typically of value 328
%     iPhU(1)                      : Index of the upper point of the pharynx for an articulation contour
%                                    Typically of value 527
% 
% Outputs
%     MA(nbSubjects)  : Column vector of the control parameter MA
%     MAC(nbSubjects) : Column vector of the centred control parameter MA
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Sizes
nbObs = size(cnts,1);

% Lower part of the pharynx
percentPha = 0.75; % = 3/4 - Empirical

% Average contour
cnt = squeeze(mean(cnts));

% Pharynx points from bottom to top
indPhaLow2Up = iPhL:((iPhL<iPhU)*2-1):iPhU;

% Restriction to the lower part of the pharynx
indPhaLower = indPhaLow2Up(1:round(percentPha*length(indPhaLow2Up)));

% Initialisation of the morphological parameter
MA = NaN(nbObs,1);

% Loop on the articulations
for iObs = 1:nbObs
    
    % Upper and lower points of the pharynx restriction
    ptPhaLowL = squeeze(cnts(iObs, indPhaLower(1),:))';
    ptPhaLowU = squeeze(cnts(iObs, indPhaLower(end),:))';
    
    % Inversion of X and Y to calculate the angle vertically and not horizontally
    Ypseudo = ptPhaLowU(1) - ptPhaLowL(1);
    Xpseudo = ptPhaLowU(2) - ptPhaLowL(2);
    
    % Angle
    % = 0 for a perfectly vertical line
    % > 0 for an obtuse angle of the vocal tract
    % < 0 for an acute angle of the vocal tract
    MA(iObs) = - atan2(Ypseudo, Xpseudo);
    
end  % for iObs = 1:nbObs

% Morphological parameter
MA;
MAC = MA - mean(MA);

end