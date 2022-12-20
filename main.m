% Morphological model of the vocal tract
% 
% Replicate the morphological model described in the following article:
% Antoine Serrurier and Christiane Neuschaefer-Rube (2023, in review)
% Morphological and acoustic modelling of the vocal tract
% Journal of the Acoustical Society of America
% 
% The code does as follows:
%   - Set the path
%   - Load the data: the morphological average-articulations + required landmarks
%   - Run the morphological model
%   - Display the performance (percentage of variance explanation and RMS reconstruction error) of the model
%   - Evaluate the morphological modelling in a leave-one-subject-out scheme 
%   - Display the performance (percentage of variance explanation and RMS reconstruction error) of the model in the leave-one-subject-out scheme
% 
% Cite:
% Antoine Serrurier and Christiane Neuschaefer-Rube (2023, in review)
% Morphological and acoustic modelling of the vocal tract
% Journal of the Acoustical Society of America
% 
% Author: Antoine Serrurier
% Date: 19/12/2022
%

% Set path
addpath(genpath('./functions/'))

% Load data
load('./data/AverageArticulations')

% Morphological Model
[scoresC, basisMorph, meanMorph, meanScores, varexTot, RMSTot, namesComp] =...
    gPCA_morphology_model(averageArticulations, iGF, iGB, iPhL, iPhU, indPhaVT, indPalVT, indVT);

% Display model performance
disp(table(varexTot*100, cumsum(varexTot)*100, RMSTot,...
    'VariableNames', {'Var. Expl. (%)', 'Cum. Var. Expl. (%)', 'Cum. RMS (cm)'}, 'RowNames', namesComp))

% Morphological model evaluation in leave-one-subject-out
[varexTotLOO, RMSTotLOO] = eval_model_LOO(averageArticulations, iGF, iGB, iPhL, iPhU, indPhaVT, indPalVT, indVT);

% Display model performance in leave-one-subject-out
disp(table(varexTotLOO*100, cumsum(varexTotLOO)*100, RMSTotLOO,...
    'VariableNames', {'Var. Expl. (%)', 'Cum. Var. Expl. (%)', 'Cum. RMS (cm)'}, 'RowNames', namesComp))

