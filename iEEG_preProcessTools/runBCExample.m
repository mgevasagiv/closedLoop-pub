bc = BadChannels;

% Example - 
%before using any of the methods - the user should define the properties for macro source
%folder, macro montage source, spike results source (if there is a wish to
%check high spikes), and load the EXP DATA
bc.sourceFolderMacro = 'D:\data_p\p485\EXP8\Denoised_Downsampled_InMicroVolt\MACRO';
bc.macroMontageFileName = 'E:\Data_p\MACRO_MONTAGE\p545\EXP3\MM485';
bc.loadHeaderExpData('D:\data_p\p485\EXP8\p485_EXP8_dataset');
bc.spikeResultsFileName = 'D:\data_p\p485\EXP8\Denoised_Downsampled_InMicroVolt\MACRO\spikesResultsShdema\SpikesResults';
disp('p485');

%see more documentation for all of the following methods in the class and DOC
%file
%the method that returns the results of the various tests on the channels
mm485 = bc.findBadChannels;
%the method that changes segments in the channel to nan according to the
%input matrix nan485
bc.correctChanWithNaN(nan485);
%get a list of areas for the patient and a map from area to channel indices
[uArea, cellAreaChan] = bc.getChansByArea;
%find which channels had high spike rate
[si485,sr485] = bc.findHighSpikes;


%% 545
%before using any of the methods - the user should define the properties for macro source
%folder, macro montage source, spike results source (if there is a wish to
%check high spikes), and load the EXP DATA
subj = 545;
exp = 3;
processedDataPath = 'E:\Data_p\p545\EXP3\';
bc.sourceFolderMacro = fullfile(processedDataPath, '\Denoised_Downsampled_InMicroVolt\MACRO');

bc.macroMontageFileName = 'E:\Data_p\MACRO_MONTAGE\p545\EXP3\MacroMontage.mat';

datasetFilename = fullfile(processedDataPath,sprintf('p%d_EXP%d_dataset.mat',subj,exp));
bc.loadHeaderExpData(datasetFilename);

bc.spikeResultsFileName = fullfile(bc.sourceFolderMacro, 'spikesResults','SpikesResults');
disp(subj);

%see more documentation for all of the following methods in the class and DOC
%file
%the method that returns the results of the various tests on the channels
mm545 = bc.findBadChannels;

%the method that changes segments in the channel to nan according to the
%input matrix nan485
bc.correctChanWithNaN(nan485);
%get a list of areas for the patient and a map from area to channel indices
[uArea, cellAreaChan] = bc.getChansByArea;
%find which channels had high spike rate
bc.spikeFileName = sprintf('MacroInterictalSpikeTimesFor_p%d_EXP%d_',subj,exp);
[si545,sr545] = bc.findHighSpikes;

% save EDF for viewing all channels
bc.saveEDF(1:20, 'p545_EXP3_ch1_20')
bc.saveEDF(21:60, 'p545_EXP3_ch21_60')
bc.saveEDF(61:100, 'p545_EXP3_ch61_100')

%% -- 
bc = BadChannels;

subj = 541;
exp = 8;
processedDataPath = 'E:\Data_p\p541\EXP8\';
bc.sourceFolderMacro = fullfile(processedDataPath, '\Denoised_Downsampled_InMicroVolt\MACRO');
bc.macroMontageFileName = 'E:\Data_p\MACRO_MONTAGE\p541\EXP8\MacroMontage.mat';

datasetFilename = fullfile(processedDataPath,sprintf('p%d_EXP%d_dataset.mat',subj,exp));
bc.loadHeaderExpData(datasetFilename);

bc.spikeResultsFileName = fullfile(bc.sourceFolderMacro, 'spikesResults','SpikesResults');
disp(subj);
%find which channels had high spike rate
bc.spikeFileName = sprintf('MacroInterictalSpikeTimesFor_p%d_EXP%d_',subj,exp);
[si541,sr541] = bc.findHighSpikes;


% save EDF for viewing all channels
bc.saveEDF(1:30, 'p541_EXP8_ch1_30')
bc.saveEDF(31:60, 'p541_EXP8_ch31_60')
bc.saveEDF(61:99, 'p541_EXP8_ch61_99')

%%
%% -- 
bc = BadChannels;

subj = 544;
exp = 8;
processedDataPath = 'E:\Data_p\p544\EXP8\';
bc.sourceFolderMacro = fullfile(processedDataPath, '\Denoised_Downsampled_InMicroVolt\MACRO');
bc.macroMontageFileName = 'E:\Data_p\MACRO_MONTAGE\p544\EXP8\MacroMontage.mat';

datasetFilename = fullfile(processedDataPath,sprintf('p%d_EXP%d_dataset.mat',subj,exp));
bc.loadHeaderExpData(datasetFilename);

bc.spikeResultsFileName = fullfile(bc.sourceFolderMacro, 'spikesResults','SpikesResults');
disp(subj);
%find which channels had high spike rate
bc.spikeFileName = sprintf('MacroInterictalSpikeTimesFor_p%d_EXP%d_',subj,exp);
[si544,sr544] = bc.findHighSpikes;

% bad data epochs were nanned via P544_patient_data_preprocess_E8.m

mm544 = bc.findBadChannels;

% save EDF for viewing all channels
bc.saveEDF(1:30, 'p544_EXP8_ch1_30')
bc.saveEDF(31:60, 'p544_EXP8_ch31_60')
bc.saveEDF(61:95, 'p544_EXP8_ch61_95')

%%
%% -- 
bc = BadChannels;

subj = 538;
exp = 3;
processedDataPath = 'E:\Data_p\p538\EXP3\';
bc.sourceFolderMacro = fullfile(processedDataPath, '\Denoised_Downsampled_InMicroVolt\MACRO');
bc.macroMontageFileName = 'E:\Data_p\MACRO_MONTAGE\p538\EXP3\MacroMontage.mat';

datasetFilename = fullfile(processedDataPath,sprintf('p%d_EXP%d_dataset.mat',subj,exp));
bc.loadHeaderExpData(datasetFilename);

bc.spikeFileName = sprintf('MacroInterictalSpikeTimesFor_p%d_EXP%d_',subj,exp);
[si538,sr538] = bc.findHighSpikes;

% save EDF for viewing all channels
bc.saveEDF(1:30, 'p538_EXP8_ch1_30')
bc.saveEDF(31:60, 'p538_EXP8_ch31_60')
bc.saveEDF(61:81, 'p538_EXP8_ch61_81')


%% -- Prep EDF for sleep scoring
bc = BadChannels;
subj = 486;
exp = 8;
processedDataPath = 'E:\Data_p\p486\EXP8\';
bc.sourceFolderMacro = fullfile(processedDataPath, '\Denoised_Downsampled_InMicroVolt\MACRO');
bc.macroMontageFileName = 'E:\Data_p\MACRO_MONTAGE\p486\EXP8\MacroMontage.mat';
datasetFilename = fullfile(processedDataPath,sprintf('p%d_EXP%d_dataset.mat',subj,exp));
bc.loadHeaderExpData(datasetFilename);

% Get a sleep-scoring montage

% save EDF for viewing all channels
bc.sourceFolderMacro = fullfile(processedDataPath, '\Denoised_Downsampled_InMicroVolt\MACRO');
obj.BRBasedStim = 1;
bc.saveEDF([97:104], 'p538_EXP8_BR_channels')
obj.BRBasedStim = 0;
bc.saveEDF([3,35,93], 'p538_EXP8_NLX_sleepScoring_channels')


