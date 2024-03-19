%% Outline of which scripts to run for formal NSD manuscript analysis

% subjectConfig.m configures subject-specific inf

% i: Preprocessing mef data
preprocNSDMef

%% Part 1: stationary maps

% Perform analyses using scalar broadband values
script01_analyzePreprocScalarBB

% % Merge across subjects
% script02_mniAcrossSubjects
% 
% %% Part 2: correlate NSD fMRI with broadbabnd
% 
% % Save beta weights for each of the 8 fMRI subjects for the shared1000 images
% saveFmriShared1000
% 
% % Save mean z-scored beta weights and save together for all subjects
% script03_meanShared1000Betas
% 
% % calculates the average vertex-wide activation for each subject and across subjects
% script03b_fmriExpected.m
% 
% % correlation maps between NSD subjects and fMRI
% script04_correlateBBFmri
% 
% % expected correlations between maximum fMRI vertex and rest of cortex, for each iEEG electrode
% script05_correlateFmriVert2Cort
% 
% %% Abstract for VSS
% 
% % Calculate stats for all electrodes BBsnr -- how many significant, anatomical locations
% scriptAbstract1_statsBBncsnrAllSubs.m
% 
% % For the different fMRI patterns, look at which NSD shared1000 stimuli are preferred
% scriptAbstract2_patternsAnalysis.m

