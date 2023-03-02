classdef BadChannels < handle
    properties
        header;
        expData;
        
        sourceFolderMicro;
        sourceFolderMacro;
        montageFileName;
        macroMontageFileName;
        spikeResultsFileName;
        samplerate = 1000;
        EDFPath;
        generateMacroMicro = 'MACRO';
        uMontage;
        macroMontage;
        BRBasedStim = 0;
        stimVal = 100;
        spikeFileName = '\spikesResultsShdema\SpikesResults';
        
        %the dictionary is from different types of "outlier" behaviour
        %to their index in the output matrix
        badChannelsDict = containers.Map({'LowThreshAmpSTD','MidThreshAmpSTD','HighThreshAmpSTD','LowThreshAmpSkewness','HighThreshAmpSkewness','LowThreshAmpKurtosis','HighThreshAmpKurtosis','HighNaN','HighSpikes','HighCorrAdj','LowCorrAdj','NegCorrAdj'},[1:12]);
        NaNThreshPerc = 0.3;
        SpikesRateThresh = 5; %per min
        
        lowAmpSTDThresh = 1.5; %>thresh*median of electrode
        midAmpSTDThresh = 2; %>thresh*median of electrode
        highAmpSTDThresh = 4; %>thresh*median of electrode
        
        lowAmpSkewThresh = 2; %>thresh*median of electrode
        highAmpSkewThresh = 4; %>thresh*median of electrode
        
        lowAmpKurtThresh = 2; %>thresh*median of electrode
        highAmpKurtThresh = 4; %>thresh*median of electrode
        
        highCorrThresh = 0.95;
        lowCorrThresh = 0.3;
        
        minContactsInArea = 6;
        timeWinForAmp = 1; %seconds
    end
    
    methods
        
        function loadHeaderExpData(obj, EXPDatafilename)
            %load exp_data and macro montage, should be run before calling saveEDF
            %or findBadChannels
            load(EXPDatafilename);
            obj.header = header;
            obj.header.processed_MICRO = obj.sourceFolderMicro;
            obj.header.processed_MACRO = obj.sourceFolderMacro;
            obj.header.montagePath = obj.montageFileName;
            obj.header.macroMontagePath = obj.macroMontageFileName;
            obj.header.samplerate = obj.samplerate;
            obj.expData = EXP_DATA;
            if ~isempty(obj.macroMontageFileName)
                obj.macroMontage = load(obj.macroMontageFileName,'MacroMontage'); obj.macroMontage = obj.macroMontage.MacroMontage;
            end
            if ~isempty(obj.montageFileName)
                obj.uMontage = load(obj.montageFileName,'Montage'); obj.uMontage = obj.uMontage.Montage;
            end
        end
        
        
        function saveEDF(obj,chansToSave, nameEDF)
            
            %A method that saves EDF file - 
            %chansToSave - an array of integers specifying what channels to save for the current patient and
            %parameters
            %nameEdf - file name of the EDF
            %The code is based on code from Maya, uses SaveEDF
            
            if nargin < 3
                nameEDF = '';
            end
            %saves EDF file - based on code from Maya, uses SaveEDF
            if strcmp(obj.generateMacroMicro,'MICRO')
                if isempty(obj.sourceFolderMicro)
                    disp 'No micro source folder provided';
                    return;
                end
                sourceFolderEDF = obj.sourceFolderMicro;
                if isempty(obj.montageFileName)
                    disp 'No micro monatage file name provided';
                    return;
                else
                    microMontage = load(obj.montageFileName); microMontage = microMontage.Montage;                 
                end
            else
                if isempty(obj.sourceFolderMacro)
                    disp 'No macro source folder provided';
                    return;
                end
                sourceFolderEDF =  obj.sourceFolderMacro;
                if isempty(obj.macroMontageFileName)
                    disp 'No macro monatage file name provided';
                    return;
                end
                %                     EEG_SIG = [{'C3'},{'C4'},{'Pz'},{'Ez'}];
                %                     str = 'MACRO';
            end
            
            
            for iChannel = 1:length(chansToSave)
                
                nChannel = chansToSave(iChannel);
                disp(sprintf('loading channel #%d/%d', iChannel, length(chansToSave)))
                
                try
                    if strcmp(obj.generateMacroMicro,'MICRO')
                        channel = microMontage(nChannel).Channel;
                        % channel = obj.montage .NS5ChannelLabels(nChannel);
                    else
                        channel = obj.macroMontage(nChannel).Channel;
                    end
                catch
                    channel = nChannel;
                end
                
                filename = fullfile(sourceFolderEDF, sprintf('CSC%d.mat',channel));
                
                if isempty(dir(filename))
                    fprintf('channel %d missing\n',channel);
                    continue
                else
                    data = load(filename,'data'); % 30kHz downsampled channels
                    data = data.data;
                end
                if iChannel == 1
                    dataLength = length(data);
                else
                    if length(data) ~= dataLength
                        disp('%s length is uncompatible - dropped')
                        data = zeros(dataLength,1);
                    end
                end
                
                dataMatEDF(iChannel,:) = data(:)';
                try
                    if strcmp(obj.generateMacroMicro,'MICRO')
                        EDFHeader.labels{iChannel} = sprintf('%d-%s',channel,microMontage(nChannel).Area);
                    else
                        EDFHeader.labels{iChannel} = sprintf('MACRO-%d-%s',channel,obj.macroMontage(nChannel).Area);
                    end
                catch
                    EDFHeader.labels{iChannel} = sprintf('channel %d',nChannel);
                end
            end
            
            %                 % Add MACRO channels - bypassed from clinical system
            %                 channel_vec = 1:length(header.NS3ChannelLabels);
            %                 for ii_cc = 1:length(channel_vec)
            %
            %                     clear data
            %                     channel = header.NS3ChannelLabels(ii_cc);
            %                     filename = fullfile(Source_folder_EDF, sprintf('CSC%d.mat',channel));
            %
            %                     if isempty(dir(filename))
            %                         fprintf('channel %d missing\n',channel) % These were inputs to closed loop
            %                         continue
            %                     else
            %                         data = load(filename,'data'); % 30kHz downsampled channels
            %                         data = data.data;
            %                     end
            %
            %                     DATA_MAT_EDF(ii_c+ii_cc,:) = data;
            %                     EDFHeader.labels{ii_c+ii_cc} = sprintf('MACRO-%d-%s',channel,header.MontageCell{ii_c+ii_cc,2});
            %                 end
            %
            %                 % add spike sorted cells to list
            %                 if header.N_spikeSortedCells > 0
            %                     Nchannels =  size(DATA_MAT_EDF,1);
            %                     Nsamples = size(DATA_MAT_EDF,2);
            %                     source_folder = header.processed_AverageRef;
            %                     load(fullfile(source_folder,sprintf('%s_spike_timestamps_post_processing',header.id)))
            %                     cnt = 0;
            %                     for jj = [3 20 25 28 44 46]; % 1:length(micro_channels_spike_summary.spike_timestamps)
            %                         cnt = cnt + 1;
            %                         area_name = micro_channels_spike_summary.unit_list_xls(jj).Location;
            %                         Cluster = micro_channels_spike_summary.unit_list_xls(jj).Cluster;
            %                         Channel = micro_channels_spike_summary.unit_list_xls(jj).Channel;
            %
            %                         Unit_spike_timestamps_ms = micro_channels_spike_summary.spike_timestamps{jj};
            %                         DATA_MAT_EDF(Nchannels+cnt,:) = zeros(Nsamples,1);
            %                         DATA_MAT_EDF(Nchannels+cnt,int64(sort(unique([Unit_spike_timestamps_ms-1, Unit_spike_timestamps_ms,Unit_spike_timestamps_ms+1])))) = int8(100);
            %                         EDFHeader.labels{Nchannels+cnt} = sprintf('%s-SU-Ch%dU%d',area_name,Channel,Cluster);
            %                     end
            %
            %                 end
            
            %%% 10 Annotation
            % header.annotation.event     structure of cells whith event name
            % header.annotation.duration   vector with event duration (seconds)
            % header.annotation.starttime  vector with event startime  (seconds)
            nChannels =  size(dataMatEDF,1);
            NStimCh = nChannels + 1;
            Nsamples = size(dataMatEDF,2);
            A = zeros(Nsamples,1);
            EDFHeader.labels{NStimCh} = sprintf('stimulation');

            if isfield(obj.expData,'stimTiming') % adding a last channel with stimulation timing
                
                if obj.BRBasedStim
                    tstim = obj.expData.stimTiming.validatedTTL_BR_msec;
                    for jStim = 1:length(tstim)
                        A(int64(sort(unique([tstim-1, tstim,tstim+1])))) = obj.stimVal;
                    end
                    dataMatEDF(NStimCh,:) = A;
                    
                else
                    
                    tstim = obj.expData.stimTiming.validatedTTL_NLX;
                    %                         else
                    %                             TTL_struct = load(fullfile(header.processedDataPath,'STIM_TTLs_NLX.mat'),'TTL_struct'); TTL_struct = TTL_struct.TTL_struct;
                    %                             tstim = TTL_struct(1).t_stim_TRAIN_END_MICRO_based_sec *obj.samplerate; % msec, NLX time
                    %                             tstim(isnan(tstim)) = [];
                    %                         end
                    for jStim = 1:length(tstim)
                        A(int64(sort(unique([tstim-1, tstim,tstim+1])))) = obj.stimVal;
                    end
                    dataMatEDF(NStimCh,:) = A;
                end
                
            elseif isfield(obj.expData,'timestamps_us')
                tstim = obj.expData.timestamps_us.X/1000;
                A(int64(tstim),1) = obj.stimVal;
                dataMatEDF(NStimCh,:) = A;
            end
            

            EDFHeader.samplerate = obj.samplerate;
            EDFHeader.annotation = [];
            EDFHeader.annotation.event = [];
            EDFHeader.annotation.duration = [];
            EDFHeader.annotation.starttime = [];
            SaveEDF(fullfile(obj.sourceFolderMacro,sprintf('%s_EXP%d_EDF_%s_1kHz_%s.edf',obj.header.id,obj.header.experimentNum,obj.generateMacroMicro,nameEDF)), dataMatEDF', EDFHeader);
        end
        
        function badChannelsMat = findBadChannels(obj)
            
            %The method runs a series of tests on the data in order to
            %recognize outlier channels. The output is a matrix where each
            %row represents a channel and each column the results of on of
            %the test (1 if positive, 0 otherwise). 
            %The dictionary which maps test name to column number is the
            %property badChannelsDict.
            %
            %The tests that are run are:
            %A. Tests that look for outlier behaviour of the amplitude of
            %data segments.
            %B. Tests that look for outlier behaviour of correlation
            %between adjacent channels.
            %C. Too many NaNs in channel.
            %D. High spike rate (if spike file name property exists)
            
            timeWin = obj.timeWinForAmp*obj.samplerate;
            
            badChannelsMat = zeros(length(obj.macroMontage),length(obj.badChannelsDict));
            
            %build montage dictionary - uArea all the areas for this
            %patient, cellAreaChan - a cell from area index to channel
            %indices of that area
            areaNames = {obj.macroMontage.Area};
            chans = [obj.macroMontage.Channel];
            [uArea,~,indsArea] = unique(areaNames);
            cellAreaChan = cell(1,length(uArea));
            for iArea=1:length(uArea)
                currInds = find(indsArea==iArea);
                cellAreaChan{iArea} = [chans(currInds)' currInds];
            end
            
            %if spike results file name is provided as a property, also checks high
            %spikes
            if ~isempty(obj.spikeResultsFileName)
                checkSpikes = true;
            else
                checkSpikes = false;
            end
            
            %The tests are run within each area - outliers are recognized
            %compared to other channels within the area
            for iArea=1:length(uArea)
                disp(['Area #',num2str(iArea),'/',num2str(length(uArea))]);
                nChans = size(cellAreaChan{iArea},1);
                %skip areas that have small number of contacts (those are
                %irrelevant areas, such as scalp EEG or eye movement)
                if nChans>=obj.minContactsInArea
                    
                    currChan = cellAreaChan{iArea}(:,1);
                    currInds = cellAreaChan{iArea}(:,2);
                    datas = [];
                    if checkSpikes
                        spikes = zeros(1,nChans);
                    end
                    
                    stdsNorms = zeros(1,nChans);
                    skewNorms = zeros(1,nChans);
                    kurtNorms = zeros(1,nChans);
                    corrsAdjElecs = zeros(1,nChans-1);
                    
                    %go over electrode indices in area and load data and spikes for each electrode
                    for iChan=1:nChans
                        try
                            data = load([obj.sourceFolderMacro,'\CSC',num2str(currChan(iChan)),'.mat']);
                            data = data.data;
                            datas(iChan,:) = data;
                        catch
                            datas(iChan,:) = nan;
                        end
                        
                        %find number of spikes for each electrode
                        if checkSpikes
                            try
                                if currChan(iChan)>9
                                    peakTimes = load([obj.spikeResultsFileName,num2str(currChan(iChan)),'_threshold=5.mat']);
                                else
                                    peakTimes = load([obj.spikeResultsFileName,'0',num2str(currChan(iChan)),'_threshold=5.mat']);
                                end
                                peakTimes = peakTimes.peakTimes;
                                spikes(iChan) = length(peakTimes);
                            catch
                                spikes(iChan) = 0;
                            end
                        end
                    end
                    
                    %the data is divided to time windows of size timeWin,
                    %the norm of each time window is calculated and
                    %parameters of the resulting distribution (norms per
                    %chan) are calculated. Channels who have distribution with outlier parameters
                    %are marked in the matrix.
                    dataL = size(datas,2);
                    nTimeWins = floor(size(datas,2)/timeWin);
                    dataLmin = dataL/obj.samplerate/60;
                    
                    normsTimeWin = zeros(nChans,nTimeWins);
                    %do checks
                    for iChan=1:nChans
                        %check if there are too many NaNs in the channel
                        if sum(isnan(datas(iChan,:)))/dataL >= obj.NaNThreshPerc
                            %mark in the output matrix that this channel is
                            %high NaN
                            badChannelsMat(currInds(iChan),obj.badChannelsDict('HighNaN')) = 1;
                            stdsNorms(iChan) = NaN;
                            corrsAdjElecs(iChan) = NaN;
                            continue;
                        end
                        
                        %calculate norm of all tie windows
                        for iTimeWin=1:nTimeWins
                            normsTimeWin(iChan,iTimeWin) = sqrt(nansum(datas(iChan,(iTimeWin-1)*timeWin+1:iTimeWin*timeWin).^2));
                        end
                        %calculated parameters of the norms' distribution -
                        %std, skewness and kurtosis. This is aimed to
                        %detect the existence of data segments with large
                        %atypical fluctuations
                        stdsNorms(iChan) = std(normsTimeWin(iChan,:));
                        skewNorms(iChan) = skewness(normsTimeWin(iChan,:));
                        kurtNorms(iChan) = kurtosis(normsTimeWin(iChan,:));
                        
                        %Another test - calculate correlation with the adjacent contact
                        if iChan < nChans
                            tmpData1 = datas(iChan,:);
                            tmpData1(isnan(tmpData1)) = 0;
                            tmpData2 = datas(iChan+1,:);
                            tmpData2(isnan(tmpData2)) = 0;
                            
                            rTmp = corrcoef(tmpData1,tmpData2);
                            corrsAdjElecs(iChan) = rTmp(2,1);
                        end
                    end
                    
                    %find the median of the channels' parameters for this
                    %area (std, skewness, kurtosis)
                    medStd = nanmedian(stdsNorms);
                    medSkew = nanmedian(skewNorms);
                    medKurt = nanmedian(kurtNorms);
                    
                    %go over the channels and find which channels have
                    %parameters who are large compared to the median and mark them in the output matrix. 
                    %For STD there are three threshold - highAmpSTDThresh,
                    %midAmpSTDThresh, lowAmpSTDThresh, for skewness and
                    %kurtosis there are two thresholds (high/low)
                    for iChan = 1:nChans
                        if ~isnan(stdsNorms(iChan))
                            if stdsNorms(iChan)>obj.highAmpSTDThresh*medStd
                                badChannelsMat(currInds(iChan),obj.badChannelsDict('HighThreshAmpSTD')) = 1;
                            else if stdsNorms(iChan)>obj.midAmpSTDThresh*medStd
                                    badChannelsMat(currInds(iChan),obj.badChannelsDict('MidThreshAmpSTD')) = 1;
                                else if stdsNorms(iChan)>obj.lowAmpSTDThresh*medStd
                                        badChannelsMat(currInds(iChan),obj.badChannelsDict('LowThreshAmpSTD')) = 1;
                                    end
                                end
                            end
                            
                            if skewNorms(iChan)>obj.highAmpSkewThresh*medSkew
                                badChannelsMat(currInds(iChan),obj.badChannelsDict('HighThreshAmpSkewness')) = 1;
                            else if skewNorms(iChan)>obj.lowAmpSkewThresh*medSkew
                                    badChannelsMat(currInds(iChan),obj.badChannelsDict('LowThreshAmpSkewness')) = 1;
                                end
                            end
                            
                            if kurtNorms(iChan)>obj.highAmpKurtThresh*medKurt
                                badChannelsMat(currInds(iChan),obj.badChannelsDict('HighThreshAmpKurtosis')) = 1;
                            else if kurtNorms(iChan)>obj.lowAmpKurtThresh*medKurt
                                    badChannelsMat(currInds(iChan),obj.badChannelsDict('LowThreshAmpKurtosis')) = 1;
                                end
                            end
                        end
                        
                        %go over channels and mark which channels are
                        %especially high correlation with next channel
                        %(i.e. are suspicious of being a duplicate),
                        %which have especially low correlation (which 
                        %may point to a bad channel) and which have
                        %negative correlation - which may point to an
                        %opposite polarity channel
                        if iChan < nChans
                            if ~isnan(corrsAdjElecs(iChan))
                                if abs(corrsAdjElecs(iChan))>obj.highCorrThresh
                                    badChannelsMat(currInds(iChan),obj.badChannelsDict('HighCorrAdj')) = 1;
                                else if abs(corrsAdjElecs(iChan))<obj.lowCorrThresh
                                        badChannelsMat(currInds(iChan),obj.badChannelsDict('LowCorrAdj')) = 1;
                                    end
                                end
                                if corrsAdjElecs(iChan)<0
                                    badChannelsMat(currInds(iChan),obj.badChannelsDict('NegCorrAdj')) = 1;
                                end
                            end
                        end
                        
                        %mark high spikes - higher rate than
                        %SpikesRateThresh, if spikes were loaded
                        if checkSpikes
                            if spikes(iChan)/dataLmin>obj.SpikesRateThresh
                                badChannelsMat(currInds(iChan),obj.badChannelsDict('HighSpikes')) = 1;
                            end
                        end
                    end
                end
            end
            
        end
        
        function saveSpikesResults(obj)
            
            %Runs spikes detections for all channels and saves the results with the
            %filenames as defined by the property spikeFileName
            
            IIS_det = SpikeWaveDetector;
            IIS_det.samplingRate = obj.samplerate;
            
            nChans = length(obj.macroMontage);
            for iChan = 1:nChans
                %load data
                disp(['Channel ',num2str(iChan),'/',num2str(nChans)]);
                try
                    data = load([obj.sourceFolderMacro,'\CSC',num2str(iChan),'.mat']);
                catch
                    disp([obj.sourceFolderMacro,'\CSC',num2str(iChan),'.mat doesnt exist']);
                    continue;
                end
                data = data.data(:)';
                %detect spikes and save
                [peakTimes, peakStats] = IIS_det.detectTimes(data, true);
                save([obj.sourceFolderMacro,obj.spikeFileName,num2str(iChan),'.mat'],'peakTimes','peakStats','IIS_det');
            end
        end
        
        function saveSpikesResultsRevPolar(obj)
            %Detects and saves spike times only for channels which have opposite polarity (the badChannel field in 
            %the Macro Montage is 4). The method assumes the data wasn’t flipped and flips it before running the 
            %detection
            IIS_det = SpikeWaveDetector;
            
            nChans = length(obj.macroMontage);
            for iChan = 1:nChans
                %check if opposite polarity
                if obj.macroMontage(iChan).badChannel == 4
                    
                    disp(['Channel ',num2str(iChan),'/',num2str(nChans)]);
                    %load data
                    try
                        data = load([obj.sourceFolderMacro,'\CSC',num2str(iChan),'.mat']);
                    catch
                        disp([obj.sourceFolderMacro,'\CSC',num2str(iChan),'.mat doesnt exist']);
                        continue;
                    end
                    
                    %flip data
                    data = data.data;
                    data = -data;
                    
                    %detect spikes and save
                    [peakTimes, peakStats] = IIS_det.detectTimes(data, true);
                    save([obj.sourceFolderMacro,obj.spikeFileName,num2str(iChan),'_RevPol.mat'],'peakTimes','peakStats','IIS_det');
                end
            end
        end
        
        function [highSpikesInds, spikeRates] = findHighSpikes(obj)
            
            %A method that finds channels with high spike rate
            %The method receives no input and returns as output:
            %highSpikesInds – the inds of the channels which passed the threshold for this patient.
            %spikeRates – the spike rates of all the channels for this patient (i.e. not just those that passed the threshold).
            
            highSpikesInds = [];
            
            nChans = length(obj.macroMontage);
            spikeRates = nan(1,nChans);
            %find data's length in order to calculate duration of the night
            data = load([obj.sourceFolderMacro,'\CSC1.mat']);
            data = data.data;
            dataDurMin = length(data)/obj.samplerate/60;
            for iChan = 1:nChans
                %load spikes times data
                disp(['Channel ',num2str(iChan),'/',num2str(nChans)]);
                try
                    spikes = load(fullfile(obj.sourceFolderMacro,[obj.spikeFileName,num2str(iChan),'.mat']));
                catch
                    disp([obj.sourceFolderMacro,obj.spikeFileName,num2str(iChan),'.mat doesnt exist']);
                    continue;
                end
                
                %find spike rate and compare to threshold
                spikes = length(spikes.peakTimes);
                spikeRates(iChan) = spikes/dataDurMin;
                if spikeRates(iChan) > obj.SpikesRateThresh
                    highSpikesInds(end+1) = iChan;
                end
            end
        end
        
        function [highSpikesInds, spikeRates, chans] = findHighSpikesRevPol(obj)
            highSpikesInds = [];
            spikeRates = [];
            chans = [];
            
            nChans = length(obj.macroMontage);
            %find data length
            data = load([obj.sourceFolderMacro,'\CSC1.mat']);
            data = data.data;
            dataDurMin = length(data)/obj.samplerate/60;
            for iChan = 1:nChans
                %load data
                if obj.macroMontage(iChan).badChannel == 4
                    disp(['Channel ',num2str(iChan),'/',num2str(nChans)]);
                    chans(end+1) = iChan;
                    try
                        spikes = load([obj.sourceFolderMacro,obj.spikeFileName,num2str(iChan),'_RevPol.mat']);
                    catch
                        disp([obj.sourceFolderMacro,obj.spikeFileName,num2str(iChan),'_RevPol.mat doesnt exist']);
                        continue;
                    end
                    spikes = length(spikes.peakTimes);
                    spikeRates(end+1) = spikes/dataDurMin;
                    if spikeRates(end) > obj.SpikesRateThresh
                        highSpikesInds(end+1) = iChan;
                    end
                end
            end
        end
        
        function correctChanWithNaN(obj,sourcefolder, nanInds)
            %This method receives indices of the segments, loads the data, marks the segments as NaN and saves the data.
            %The input is a matrix (nanInds) with the following format:
            %Each row includes 3 numbers: <channel index> <start point of bad segment in seconds> <duration of bad segment in seconds)
            
            chansToCorrect = unique(nanInds(:,1));
            %go over the channels in nanInds
            for iChan = 1:length(chansToCorrect)
                disp(['Channel ',num2str(iChan),'/',num2str(length(chansToCorrect))]);
                currFileName = [sourcefolder,'\CSC',num2str(chansToCorrect(iChan)),'.mat'];
                try
                    mfile = matfile(currFileName,'Writable',true);
                    data = mfile.data;
                catch
                    disp([currFileName,' doesn''t exist']);
                    continue;
                end
                %mark bad segments as nan
                indsChan = nanInds(nanInds(:,1)==chansToCorrect(iChan),:);
                for iSeg = 1:size(indsChan,1)
                    data(indsChan(iSeg,2)*obj.samplerate+1:indsChan(iSeg,2)*obj.samplerate+indsChan(iSeg,3)*obj.samplerate) = NaN;
                end
                % save(currFileName,'data'); % this overrided the current
                                             % file structure
                mfile.data = data;
                mfile.nanned = 1;
            end
        end
        
        function [uArea, mapAreaChan] = getChansByArea(obj)
            % This method returns two variables:
            % uArea - a list of all the areas for the current patients
            % mapAreaChan – a map from each area to its list of channel indices. i.e. mapAreaChan(‘RAH’) will return a 
            % list of the RAH channel indices.

            areaNames = {obj.macroMontage.Area};
            chans = [obj.macroMontage.Channel];
            [uArea,~,indsArea] = unique(areaNames);
            mapAreaChan = containers.Map;
            for iArea=1:length(uArea)
                mapAreaChan(uArea{iArea}) = chans(indsArea==iArea);
            end
        end
        
        function [detSR, undetSR] = plotSR(obj, spikeRatesPerChannel, humanDecMap)
            [uArea, cellAreaChan] = obj.getChansByArea;
            areasDetected = humanDecMap.keys;
            
            detChans = [];
            undetChans = [];
            
            for iArea = 1:length(uArea)
                currChans = cellAreaChan{iArea};
                if length(currChans) < 6
                    continue;
                end
                if ismember(uArea{iArea},areasDetected)
                    detChans = [detChans currChans(humanDecMap(uArea{iArea}))];
                    undetChans = [undetChans setdiff(currChans, currChans(humanDecMap(uArea{iArea})))];
                else
                    undetChans = [undetChans currChans];
                end
            end
            
            detSR = spikeRatesPerChannel(detChans);
            undetSR = spikeRatesPerChannel(undetChans);
            
            plot(zeros(1,length(detSR)),detSR,'*');
            hold all;
            plot(ones(1,length(undetSR)),undetSR,'*');
            xlim([-1 2]);
            line([-1 2], [5 5], 'color', 'k');
            legend('SR of chans human-positive','SR of chans human-negative')
        end
    end
end
