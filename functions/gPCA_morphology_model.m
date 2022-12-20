function [scoresC, basisVectors, meanCnts, meanScores, varexTot, RMSTot, namesComp,...
    namesCompLong, nbComp, scores, coefsPalMPA, meanPalMPA, coefsPalMPC, meanPalMPC] =...
    gPCA_morphology_model(averageArticulations, iGF, iGB, iPhL, iPhU, indPhaVT, indPalVT, indFocus)
% 
% Construction of a morphological model of the vocal tract by means of guided PCA according to the following article:
% 
% Antoine Serrurier and Christiane Neuschaefer-Rube (2023, in review)
% Morphological and acoustic modelling of the vocal tract
% Journal of the Acoustical Society of America
% 
% Inputs
%     averageArticulations(nbSubjects,nbPts,nbDim) : Morphological average-articulations
%                                                    Typically of size 41 x 1692 x 2
%     iGF(1)                                       : Index of the point corresponding to the anterior of the glottis for an articulation contour
%                                                    Typically of value 1631
%     iGB(1)                                       : Index of the point corresponding to the posterior of the glottis for an articulation contour
%                                                    Typically of value 1650
%     iPhL(1)                                      : Index of the lower point of the pharynx for an articulation contour
%                                                    Typically of value 328
%     iPhU(1)                                      : Index of the upper point of the pharynx for an articulation contour
%                                                    Typically of value 527
%     indPhaVT(nbPtsVT)                            : Indices of the vocal tract points corresponding to the pharynx for an articulation contour
%                                                    Typically of length 200
%     indPalVT(nbPtsPalVT)                         : Indices of the vocal tract points corresponding to the hard palate for an articulation contour
%                                                    Typically of length 108
%     indFocus(nbPtsFocus)                         : Indices of the points for which the model if optimised for an articulation contour (=vocal tract points)
%                                                    Typically of length 1036
% 
% Outputs
%     scoresC(nbSubjects,nbComp)        : Centred control parameters
%     basisVectors(nbComp,nbPts,nbDim)  : Morphologcal basis vectors
%     meanCnts(nbPts,nbDim)             : Mean average-articulation
%     meanScores(1,nbComp)              : Mean of the non-centred control parameters
%     varexTot(nbComp)                  : Percentage of variance explanation per component (between 0 and 1)
%     RMSTot(nbComp)                    : Cumulated RMS reconstruction error per component (in cm)
%     namesComp({nbComp})               : Names of the morphologcal components of the model
%     namesCompLong({nbComp})           : Long names of the morphologcal components of the model
%     nbComp(1)                         : Number of morphological components
%     scores(nbSubjects,nbComp)         : Non-centred control parameters
%     coefsPalMPA(1, nbPtsPalVT)        : Coefficients of the raw PCA performed on the X-values of the palate points for the calculation of MPA
%     meanPalMPA(nbPtsPalVT, 1)         : Mean of the X-values of the palate points over the subjects for the calculation of MPA
%     coefsPalMPC(1, nbPtsPalVT, nbDim) : Coefficients of the raw PCA performed on the coordinates of the palate points  for the calculation of MPC
%     meanPalMPC(nbPtsPalVT, nbDim)     : Mean of the coordinates of the palate points over the subjects for the calculation of MPC
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

%% ======================================================================
% Constants / Global variables

% UT point
UT = [5, 10];

% Names components
namesComp = {'MX', 'MY', 'MA', 'MPA', 'MPC'};
namesCompLong = {'Horizontal Length', 'Vertical Height', 'Angle',...
    'Palate Anteriority', 'Palate Concavity'};

% Dimensions
nbObs = size(averageArticulations, 1);
nbPts = size(averageArticulations, 2);
nbPtsFocus = length(indFocus);
nbDim = size(averageArticulations, 3);

% Number of components
nbComp = length(namesComp);

% Mean of the data
meanCnts = squeeze(mean(averageArticulations));

%% ======================================================================
% Iterative guided PCA

% Initialisation of the output arguments
scores = NaN(nbObs, nbComp); % Scores
scoresC = NaN(nbObs, nbComp); % Centred scores
basisVectors = NaN([nbComp, nbPts, nbDim]); % Basis vectors

% Initialisation of the residue
cntsRes = averageArticulations - permute(repmat(meanCnts, [1,1,nbObs]),[3,1,2]);

% Loop on the components
for iComp = 1:nbComp

    % Considered component
    nameComp = namesComp{iComp};

    % Residue + mean
    cntsRes_noC = cntsRes + permute(repmat(meanCnts,[1,1,nbObs]),[3,1,2]);

    % Control parameter calculation
    switch nameComp
        case 'MX'  % switch nameComp
            [scoreComp, scoreCompC] = gPCA_getMX(cntsRes_noC, UT(2), indPhaVT);
        case 'MY'  % switch nameComp
            [scoreComp, scoreCompC] = gPCA_getMY(cntsRes_noC, iGF, iGB);
        case 'MA'  % switch nameComp
            [scoreComp, scoreCompC] = gPCA_getMA(cntsRes_noC, iPhL, iPhU);
        case 'MPA'  % switch nameComp
            [scoreComp, scoreCompC, coefsPalMPA, meanPalMPA] = gPCA_getMPA(cntsRes_noC, indPalVT);
        case 'MPC'  % switch nameComp
            [scoreComp, scoreCompC, coefsPalMPC, meanPalMPC] = gPCA_getMPC(cntsRes_noC, indPalVT);
    end  % switch nameComp

    % Linear regression to get the component basis vector of the focus points
    basisVectorFocus = regress_Data_Scores_2_BasisVector(cntsRes(:,indFocus,:), scoreCompC);

    % Linear regression to get the component basis vector of all points
    basisVector = regress_Data_Scores_2_BasisVector(cntsRes, scoreCompC);

    % Replace the coefficients of the focus points by their values (slightly more optimal)
    basisVector(indFocus,:) = basisVectorFocus;

    % Update the residue for all points
    cntsResEst = predict_Scores_BasisVectors_2_Data(scoreCompC, permute(basisVector, [3,1,2]), zeros(nbPts,nbDim));
    cntsRes = cntsRes - cntsResEst;
    
    % Update the output arguments
    scores(:,iComp) = scoreComp;
    scoresC(:,iComp) = scoreCompC;
    basisVectors(iComp,:,:) = basisVector;

end  % for iComp = 1:nbComp

% Provide the mean of the scores in output
meanScores = mean(scores);

%% ======================================================================
% Evaluation on the focus points

% Initialisation of the error variables
varexTot = NaN(nbComp, 1); % Percentage of variance explained per component
RMSTot = NaN(nbComp, 1); % Cumulated RMS reconstruction error per component

% Restriction of the data and model on the focus points
basisVectorsFocus = basisVectors(:,indFocus,:);
meanCntsFocus = meanCnts(indFocus,:);
cntsFocus = averageArticulations(:,indFocus,:);

% Initialisation of the predicted data
cntsFocusEst = permute(repmat(meanCntsFocus, [1,1,nbObs]),[3,1,2]);

% Reshape
matCntsFocus = reshape(cntsFocus, nbObs, nbPtsFocus*nbDim);
matCntsFocusEst = reshape(cntsFocusEst, nbObs, nbPtsFocus*nbDim);

% Initial variance of the data
% Per point
varCnts = var(matCntsFocus, 1);
% Total
varTotCnts = sum(varCnts);

% Loop on the components
for iComp = 1:nbComp

    % Data prediction of this comp only
    cntsFocusCEst_comp = predict_Scores_BasisVectors_2_Data(scoresC(:,iComp), basisVectorsFocus(iComp,:,:), zeros(nbPtsFocus,nbDim));

    % Reshape
    matCntsFocusCEst_comp = reshape(cntsFocusCEst_comp, nbObs, nbPtsFocus*nbDim);

    % Variance of this component
    % Per point
    varCntsEst_comp = var(matCntsFocusCEst_comp, 1);
    % Total
    varTotCntsEst_comp = sum(varCntsEst_comp);

    % Percentage of variance exlained for the component
    varexTot(iComp) = varTotCntsEst_comp / varTotCnts;

    % Update the data prediction
    matCntsFocusEst = matCntsFocusEst + matCntsFocusCEst_comp;

    % RMS error
	RMSTot(iComp) = sqrt(mean(var(matCntsFocusEst - matCntsFocus, 1)));
end  % for iComp = 1:nbComp

end

