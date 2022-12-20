function cntsEst = gPCA_predict_Data_BasisVectors_2_Data(cnts, basisVectors, meanCnts, meanScores,...
    coefPalMPA, meanPalMPA, coefPalMPC, meanPalMPC, namesComp, indPhaVT, iGF, iGB, iPhL, iPhU, indPal)
% 
% Reconstruct the input data by means of the morphological model provided in input.
% It calculates iteratively for the input data the control parameters and
% multiply the control parameters with their corresponding basis vectors.
% 
% Inputs
%     cnts(nbSubjects,nbPts,nbDim)     : Articulation contours to approach with the model
%                                        Typically of size 41 x 1692 x 2
%     basisVectors(nbComp,nbPts,nbDim) : Morphological basis vectors
%                                        Typically of size 5 x 1692 x 2
%     meanCnts(nbPts,nbDim)            : Mean average-articulation of the morphological model
%                                        Typically of size 1692 x 2
%     meanScores(1,nbComp)             : Mean of the non-centred control parameters of the morphological model
%                                        Typically of size 1 x 5
%     coefsPalMPA(1, nbPtsPal)         : Coefficients of the raw PCA performed on the X-values of the palate points for the calculation of MPA
%                                        Typically of size 1 x 108
%     meanPalMPA(nbPtsPal, 1)          : Mean of the X-values of the palate points over the subjects for the calculation of MPA
%                                        Typically of size 108 x 1
%     coefsPalMPC(1, nbPtsPal, nbDim)  : Coefficients of the raw PCA performed on the coordinates of the palate points  for the calculation of MPC
%                                        Typically of size 1 x 108 x 2
%     meanPalMPC(nbPtsPal, nbDim)      : Mean of the coordinates of the palate points over the subjects for the calculation of MPC
%                                        Typically of size 108 x 2
%     namesComp({nbComp})              : Vector of cells of the names of the morphologcal components of the model
%                                        Typically of size 1 x 5
%     indPhaVT(nbPtsVT)                : Indices of the vocal tract points corresponding to the pharynx for an articulation contour
%                                        Typically of length 200
%     iGF(1)                           : Index of the point corresponding to the anterior of the glottis for an articulation contour
%                                        Typically of value 1631
%     iGB(1)                           : Index of the point corresponding to the posterior of the glottis for an articulation contour
%                                        Typically of value 1650
%     iPhL(1)                          : Index of the lower point of the pharynx for an articulation contour
%                                        Typically of value 328
%     iPhU(1)                          : Index of the upper point of the pharynx for an articulation contour
%                                        Typically of value 527
%     indPal(nbPtsPal)                 : Indices of the vocal tract points corresponding to the hard palate for an articulation contour
%                                        Typically of length 108
% 
% Outputs
%     cntsEst(nbSubjects,nbPts,nbDim) : Articulation contours estimated with the model
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% UT point
UT = [5, 10];

% Sizes
nbObs = size(cnts,1);
nbPts = size(cnts,2);
nbDim = size(cnts,3);
nbComp = size(basisVectors,1);

% Initialise the scores (control parameters)
scores = NaN(nbObs, nbComp);

% Initialise the reconstructed data
cntsEst = permute(repmat(meanCnts, [1,1,nbObs]),[3,1,2]);

% Initialise the residue
res = cnts - cntsEst;

% Iterative reconstruction
for iComp = 1:nbComp
    
    % Residue plus mean
    res_noC = res + permute(repmat(meanCnts, [1,1,nbObs]),[3,1,2]);

    % Calculate the control parameter
    switch namesComp{iComp}
        case 'MX'  % switch namesComp{iComp}
            scoresRaw = gPCA_getMX(res_noC, UT(2), indPhaVT);
        case 'MY'  % switch namesComp{iComp}
            scoresRaw = gPCA_getMY(res_noC, iGF, iGB);
        case 'MA'  % switch namesComp{iComp}
            scoresRaw = gPCA_getMA(res_noC, iPhL, iPhU);
        case 'MPA'  % switch namesComp{iComp}
            scoresRaw = regress_Data_BasisVector_2_Scores(res_noC(:,indPal,1), coefPalMPA', meanPalMPA);
        case 'MPC'  % switch namesComp{iComp}
            scoresRaw = regress_Data_BasisVector_2_Scores(res_noC(:,indPal,:), squeeze(coefPalMPC), meanPalMPC);
    end  % switch namesComp{iComp}
    
    % Centred control parameter
    scores(:,iComp) = scoresRaw - meanScores(iComp);

    % Reconsctruction with the current component
    cntsEst_comp = predict_Scores_BasisVectors_2_Data(scores(:,iComp), basisVectors(iComp,:,:), zeros(nbPts,nbDim));

    % Update the reconstructed points
    cntsEst = cntsEst + cntsEst_comp;

    % Update the residue
    res = cnts - cntsEst;
end  % for iComp = 1:nbComp

end