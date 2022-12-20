function [varexTot, RMSTot] = eval_model_LOO(averageArticulations, iGF, iGB, iPhL, iPhU, indPhaVT, indPalVT, indFocus)
%
% Evaluation of vocal tract morphological modelling in a leave-one-subject-out scheme
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
%     varexTot(nbComp) : Percentage of variance explanation per component (between 0 and 1)
%     RMSTot(nbComp)   : Cumulated RMS reconstruction error per component (in cm)
% 
% Author : Antoine Serrurier
% Date: 19/12/2022

% Apply once the model on all data to get the number of components
[~, ~, ~, ~, ~, ~, ~, ~, nbComp] =...
    gPCA_morphology_model(averageArticulations, iGF, iGB, iPhL, iPhU, indPhaVT, indPalVT, indFocus);

% Number of subjects
nbSubjects = size(averageArticulations,1);

% Initialiation of the reconstructed average articulations, per subject out and per component 
averageArticulationsEst = NaN([nbComp, size(averageArticulations)]);

% Leave one subject out in each loop
for iSubOut = 1:nbSubjects
    
    % Remaining subjects
    indSubTrain = setdiff(1:nbSubjects, iSubOut);

    % Build model on the remaining subjects
    [~, basisMorph, meanMorph, meanScores, ~, ~, namesComp,...
        ~, ~, ~, coefsPalMPA, meanPalMPA, coefsPalMPC, meanPalMPC] =...
        gPCA_morphology_model(averageArticulations(indSubTrain,:,:), iGF, iGB, iPhL, iPhU, indPhaVT, indPalVT, indFocus);

    % Estimate the subject out from the model, component by component
    for iComp = 1:nbComp
        averageArticulationsEst(iComp,iSubOut,:,:) = gPCA_predict_Data_BasisVectors_2_Data(averageArticulations(iSubOut,:,:),...
            basisMorph(1:iComp,:,:), meanMorph, meanScores, coefsPalMPA, meanPalMPA, coefsPalMPC, meanPalMPC,...
            namesComp, indPhaVT, iGF, iGB, iPhL, iPhU, indPalVT);
    end  % for iPred = 1:size(ptsEstVTLoo,2)

end  % for iSubOut = 1:nbSubjectsAll

% Calculate the performance, overall and per component 
varexTot_cum = NaN(nbComp,1); %  Initialisation of the cum. variance explanation per component
RMSTot = NaN(nbComp,1); % Initialisation of the cum. RMS reconstruction error per component
for iComp = 1:nbComp
    [varexTot_cum(iComp), RMSTot(iComp)] = variance_rms(averageArticulations(:,indFocus,:), squeeze(averageArticulationsEst(iComp,:,indFocus,:)));
end  % for iPred = 1:size(ptsEstVTLoo,2)
varexTot = [varexTot_cum(1); diff(varexTot_cum)]; %  Variance explanation per component

end

