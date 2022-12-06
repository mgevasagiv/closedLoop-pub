classdef AnalyzeStimulationEffect_pub < handle
    properties
        samplingRate = 1000;
        
        %filtering constants
        defaultFilterOrder = 1;
        nanWarning = 0.01;
        
        lowLimitSW = 0.5; %Hz
        highLimitSW = 2; % % MGS - This is different than our SW-detection scheme, but aligns to the SW detection in realtime
        minUpPhase = pi/4;
        maxUpPhase = pi;
        minDnPhase = -3*pi/4;
        maxDnPhase = 0;
        
        minNumUpDown = 50;
        
        % Defining a SHAM-control set
        desiredPercentOfSlowWaves = 40;
        
        %single events co-occurance constants, based on Maingret 2016
        minDistSWSpindle = 0; %ms ++
        maxDistSWSpindle = 1300; %ms ++
        minDistRippleSW = 50; %ms ++
        maxDistRippleSW = 500; %ms ++
        avgFigBinning = 100; %ms ++
        timeToCheckAfterSW = 3000; %ms ++
        timeToCheckAfterSp = 3000; %ms
        timeToCheckAfterRip = 200; %ms
        timeToCheckBeforeSW = 0; %ms ++
        timeToCheckBeforeSp = 0; %ms
        timeToCheckBeforeRip = 0; %ms
        controlDistRip = 1500; %ms, relevant to ripples only
        nControls = 1000;
        
        % indexing the all*Event 4-way matrices
        event_index_im = 1;
        event_index_im_control = 2;
        event_index_persistent_control = 3;
        event_index_persistent = 4;
        event_index_persistent_matched_length_control = 6;
        
        minSpikeRateToIncludeUnit = 0.1; %Hz ++
        windowSpikeRateAroundStim = 500; %ms ++
        windowSpikeRateForComparison = 500; %ms - for comparing between stim and control
        controlDistForStim = 1000; %ms ++
        firingRateWinSize = 10; %ms ++
        
        % coherence check
        slowWavesRange = [0.5 4];
        spindlesRange = [9 16];
        
        % phaseLock analysis
        minSpikesToIncludeUnitPhaseAnalysis = 30; % per epoch
        minRatetoIncludeUnitPhaseAnalysis = 0.1;
        spikeMinDistIIS = 500;
        minT_spike = 500; %ms
        binSizeRad = 20*(pi/180);
        minPhaseLockingPToIncludeInAnalysis = 0.1;
        phaseLockingBuffer_ms = 60e3; % ms
        pval_th = 0.05;
        rsq_th = 0.2;
        phaseLockbaselineInd = 1; % Pre-stimulation
        phaseLockStimSessionID = 2; % Using all stimulation epochs pulled together
        phaseLockconditions_str  = {'pre-session','all stim','60s post stim-blocks','first pause','all pause','60s pre stim-blocks'}
        
        % bootStrapping phase-lock
        phaseLocking_BS_N = 10;
        
        %cross-corr analysis
        crossCorrBuffer_ms = 300e3; % ms
        minNumSpikesCrossCor = 30;
        crossCorWindowSpikeRateAfterStim = 2500; %ms ++
        crossCorfsz = 10;
        ind_for_auc_comparison = 190:210; % assuming 10mS bins, [-100,100] msec around peak
        
        % cross cor plots
        XLIM_auc_comparison = [170,230];
        
        % TFR analysis
        timeBeforeAfterEventStim = 2.5;  % sec
        timeForBaselineStim = 1; % sec
        minNCyclesStim = 5;
        
        % ROITimeRegionSWSp = [-1 1]; %seconds
        ROIFreqRegionSWSp = [9 16]; %frequency
        ROIFreqRegionControl = [20 27]; %frequency
        
        % spectrum constants
        freqRangeSWSp = [5:30]; %Hz
        
        % plotting phase in [0-360]
        minUpPhase_360 = pi/2;
        maxUpPhase_360 = 3*pi/2;
        
        %final analysis constants
        shortTimeRangeAfterStim = 3; %seconds ++
        midTimeRangeAfterStim = 60; %seconds
        timeBeforeStim = 10; %ms
        nColsInFig = 4;
        
        pThresh = 0.05;
        dataFilePrefix = 'CSC';
        sleepEpochs = 1;
        
        %presentation constants
        nSTDofDataToShow = 2;
        
        numBinsPhase = 10;
    end
    
    methods
        
        % Collect stimulation stats for Supp. table S3
        function results = collectStimStats(obj, runData, fileNameResults)
            if nargin < 2
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName,' ',num2str(iPatient),'/',num2str(nPatients)]);
                
                ptVec(iPatient) = str2num(runData(iPatient).patientName(2:end));
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                header = expData.header;
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX; stimTimes = floor(stimTimes(:))';
                
                %load sleep scoring
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                
                %check only stim times during NREM
                if ~isempty(sleepScoring)
                    stim_times_ms = stimTimes(sleepScoring(round(stimTimes))==1);
                end
                
                
                EXP_num(iPatient) = header.experimentNum;
                Nstim(iPatient) = length(expData.stimTiming.validatedTTL_NLX);
                Nstim_NREM(iPatient) = length(stim_times_ms);
                NstimBlocks(iPatient)  = (size(expData.Session_start_end_msec,1)-1)/2;
                tSessionMin(iPatient)  = (expData.Session_start_end_msec(end-1,2)-expData.Session_start_end_msec(2,1))/60e3;
                
            end
            
            [ptVec, id] = sort(ptVec);
            EXP_num = EXP_num(id);
            NstimBlocks = NstimBlocks(id);
            Nstim = Nstim(id);
            Nstim_NREM = Nstim_NREM(id);
            tSessionMin(id) = tSessionMin;
            
            results.ptVec = ptVec;
            results.EXP_num = EXP_num;
            results.Nblocks = NstimBlocks;
            results.Nstim = Nstim;
            results.Nstim_NREM = Nstim_NREM;
            results.tSessionMin = tSessionMin;
            
        end
        
        % Generate an equivalent set of timestamps with a similar distribution as
        % the stim one - relative to up-states in pre-session
        % SHAM_stimDist - choose SHAM points preserving the distance from peaks
        function genControlTimestamps(obj, runData, saveFilename)
            
            nPatients = length(runData);
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName,' ',num2str(iPatient),'/',num2str(nPatients)]);
                filename_str = runData(iPatient).ExpDataFileName(end-20:end-12);
                filename = [saveFilename,sprintf('_%s_stim_control_sets.mat',filename_str)];
                
                %                 if ~isempty(dir(filename))
                %                     disp('timestamps file exists')
                %                     continue
                %                 end
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX; stimTimes = floor(stimTimes(:))';
                
                %probe chan
                probeChan = runData(iPatient).probeChan;
                
                %load probe data
                try
                    probeData = load([runData(iPatient).DataFolder,'\',obj.dataFilePrefix ,num2str(probeChan),'.mat']);
                    probeData = probeData.data;
                catch
                    disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(probeChan),'.mat doesnt exist']);
                    continue;
                end
                
                %load sleep scoring
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                
                %check only stim times during NREM
                if ~isempty(sleepScoring)
                    stim_times_ms = stimTimes(sleepScoring(round(stimTimes))==1);
                end
                results(iPatient).stimTimes_NREM = stim_times_ms;
                
                nStims = length(stim_times_ms);
                
                % get interictal spikes
                try
                    SpikesFileNames = [runData(iPatient).SpikesFileNames,num2str(probeChan),'.mat'];
                    ISIpeakTimes = load(SpikesFileNames);
                    ISIpeakTimes = ISIpeakTimes.peakTimes;
                catch
                    disp([SpikesFileNames,'.mat doesn''t exist']);
                    disp('running spikes');
                    IIS_det = SpikeWaveDetector;
                    [ISIpeakTimes, ~] = IIS_det.detectTimes(probeData, true);
                end
                
                ISIpeakTimes(ISIpeakTimes<obj.minT_spike) = obj.minT_spike+1;
                probeData(floor(ISIpeakTimes) - obj.minT_spike: floor(ISIpeakTimes) + obj.minT_spike)= 0;
                
                for ii = 1:length(stim_times_ms); probeData(stim_times_ms(ii)-100:stim_times_ms(ii)+100) = 0;end
                
                % low pass filtering
                dataFilteredSW = obj.bandpass(probeData, obj.lowLimitSW, obj.highLimitSW);
                nanInds = isnan(dataFilteredSW);
                dataFilteredSW(nanInds) = 0;
                dataFilteredSW(sleepScoring ~= 1) = NaN;
                
                
                
                % SHAM v0 - choose time-points based on ISI of stims
                t0 = floor(expData.Session_start_end_msec(1,1));
                tend = floor(expData.Session_start_end_msec(end,2));
                
                peaks = getSW_local(obj, dataFilteredSW(t0:tend)); % get the *local* peaks
                if isempty(peaks)
                    error('no peaks found')
                end
                events_ms = floor(t0 + peaks);
                allPeaks = dataFilteredSW(events_ms);
                
                stim_diff_event_ms = []; downState_peak_stim = [];
                for ii = 1:length(stim_times_ms)
                    tstim = stim_times_ms(ii);
                    if (isnan(probeData(floor(tstim)-1e3:floor(tstim)+1e3)))
                        continue
                    end
                    indices = events_ms < tstim;
                    [t_diff idx] = min(tstim - events_ms(indices));
                    downState_peak_stim = [downState_peak_stim, dataFilteredSW(events_ms(idx))];
                    stim_diff_event_ms = [stim_diff_event_ms, t_diff];
                end
                
                % get 'pause' sessions peaks
                pausePeaks = [];
                timeFromStimBuffer = 60e3;
                for i_ss = 1:2:length(expData.Session_start_end_msec)
                    t0 = expData.Session_start_end_msec(i_ss,1);
                    tend = expData.Session_start_end_msec(i_ss,2);
                    
                    indices = (events_ms > (t0 + timeFromStimBuffer)) &  (events_ms < tend);
                    pausePeaks = [pausePeaks, events_ms(indices)'];
                end
                indices = randi(length(pausePeaks),1,length(stim_diff_event_ms));
                shamTimestamps = pausePeaks(indices) + stim_diff_event_ms;
                downStatePausePeaks =  dataFilteredSW(pausePeaks(indices));
                

                A = [];
                for ii = 1:length(stim_diff_event_ms);
                    T = stim_times_ms(ii) - stim_diff_event_ms(ii);
                    A(ii,:)  = dataFilteredSW(T-1000:T+1000);
                end
                 B = [];
                for ii = 1:length(stim_times_ms);
                    T = stim_times_ms(ii);
                    B(ii,:)  = dataFilteredSW(T-1000:T+1000);
                end
                mean_filtered = nanmean(A);
                std_filtered = nanstd(A);
                stim_control_sets.stim_diff_event_ms = stim_diff_event_ms;
                stim_control_sets.shamTimestamps = shamTimestamps;
                stim_control_sets.allPeaks = allPeaks;
                stim_control_sets.upState_peak_stim = downState_peak_stim;
                stim_control_sets.pausePeaks = downStatePausePeaks;
                stim_control_sets.mean_filtered_MACRO_PROBE = nanmean(A);
                stim_control_sets.std_filtered_MACRO_PROBE = nanstd(A);
                save(filename,'stim_control_sets')
                
                
                str = runData(iPatient).ExpDataFileName(end-20:end-12);
                figName = sprintf('%s_stimTriggered_probeLFP',str);
                f0 = figure('Name', figName,'NumberTitle','off');
                % Some WYSIWYG options:
                set(gcf,'DefaultAxesFontSize',8);
                set(gcf,'DefaultAxesFontName','arial');
                set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 10 10]); % this size is the maximal to fit on an A4 paper when printing to PDF
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
                colormap('jet');
                axes('position',[0.1,0.1,0.6,0.85])
                
                subplot(2,2,1:2)
                A = [];
                for ii = 1:length(stim_diff_event_ms);
                    T = stim_times_ms(ii) - stim_diff_event_ms(ii);
                    A(ii,:)  = dataFilteredSW(T-1000:T+1000);
                end
                mean_filtered = nanmean(A);
                std_filtered = nanstd(A);
                shadedErrorBar([-1000:1000], mean_filtered, std_filtered/sqrt(nStims),'lineprops','-b');
                YLIM = get(gca,'ylim');
                title('Macro signal eliciting stim (aligned on detected peak)')
                ylabel('Filtered iEEG (uV, [0.5,2]Hz)')
                xlabel('ms')
                
                axis([-1000 1000,floor(YLIM)]);
                YLIM = get(gca,'ylim');
                set(gca,'xtick',[-1000,0,500,1000],'ytick',[floor(YLIM(1)),0,floor(YLIM(2))]);
                if iPatient == 7
                    axis([-1000 1000,-150 300]);
                    set(gca,'xtick',[-1000,0,500,1000],'ytick',[-150,0,300]);
                end
                
                subplot(2,2,3)
                histogram(stim_diff_event_ms,[0:25:500])
                YLIM = get(gca,'ylim');
                if YLIM(2) < 20
                    set(gca,'xtick',[0,250,500],'ytick',[0,10,20],'ylim',[0,20]);
                else
                    set(gca,'xtick',[0,250,500],'ytick',[0,20,40],'ylim',[0,40]);
                end
                xlabel('ms')
                ylabel('counts')
                title({'Time from closest','down-state peak'})
                
                subplot(2,2,4)
                bins = [0:10:floor(max([downStatePausePeaks,downState_peak_stim]))];
                histogram(downStatePausePeaks,bins,'facecolor',[0.7,0.7,0.7]);
                hold on
                histogram(downState_peak_stim,bins,'facecolor',[0.7,0,0]);
                [h,p] = ttest2(downState_peak_stim,downStatePausePeaks);
                if max(bins) > 250
                    set(gca,'xlim',[0,500],'xtick',[0,500],'ytick',[0,20,40]);
                else
                    set(gca,'xlim',[0,250],'xtick',[0,250],'ytick',[0,20,40]);
                end
                if iPatient == 7
                    axis([150 650,0 30]);
                    set(gca,'xtick',[150,650],'ytick',[0,30]);
                end
                title(sprintf('down-state peaks, p = %2.2f',p))
                xlabel('uV')
                ylabel('counts')
                legend('pause epochs','stim epochs')
                
                outputFigureFolder = 'E:\Data_p\ClosedLoopDataset\stimPhaseFigures\perPt';
                a  = get(gcf);
                res = 400;
                eval(['print ', [outputFigureFolder,'\',a.Name], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
                close(gcf)
            end
            
        end
        
        
        % Very simple peak detector - no amplitude threshold
        % (every zero crossing is a candidate)
        % based on - detect_SW_positive_dealWithSAW()
        function peaks = getSW_local(obj, data_vec)
            
            pos_index=zeros(length(data_vec),1);
            pos_index(find(data_vec>0))=1; %index of all positive points for EEG
            difference=diff(pos_index); poscross=find(difference==1) ; negcross=find(difference==-1); %find neg ZX and pos ZX
            EEGder = obj.meanfilt(diff(data_vec),5); %meanfilt is a function that uses a 5 sample moving window to smooth derivative
            pos_index=zeros(length(EEGder),1);
            pos_index(find(EEGder>0.1))=1; %index of all positive points above minimum threshold
            difference=diff(pos_index);
            peaks=find(difference==-1)+1; troughs=find(difference==1)+1; %find pos ZX and neg ZX of the derivative (the peaks & troughs)
            peaks( data_vec(peaks)<0 | isnan( data_vec( peaks)))=[]; % rejects peaks below zero and troughs above zero
            troughs(data_vec(troughs)>0 | isnan(data_vec( troughs)))=[]; % rejects peaks below zero and troughs above zero
            
            %% Now work separately on good or bad waves -
            %% First - GOOD WAVES
            %% select a subset of slow waves with highest amplitudes
            numOfSlowWaves = size(peaks,1);
            allAmplitudes = data_vec(peaks);
            [sortedAmplitudes,sortedIndices] = sort(allAmplitudes, 'descend');
            slowWavesSortedByAmplitude = peaks(sortedIndices, :);
            cutoffNumber = round((obj.desiredPercentOfSlowWaves/100) * numOfSlowWaves);
            peaks = sort(slowWavesSortedByAmplitude(1:cutoffNumber, :));
            
        end
        
        % same local meanfilt as used in SW-detector
        function [filtdata]=meanfilt(obj, datatofilt,pts)
            
            if length(datatofilt)>=pts
                filtdata=[];
                ptsaway=floor(pts/2); isEven = (ptsaway*2 == pts);
                filtdata([1:pts])=datatofilt([1:pts]);
                filtdata([length(datatofilt)-(pts-1):length(datatofilt)])=datatofilt([length(datatofilt)-(pts-1):length(datatofilt)]);
                for wndw=pts-ptsaway:length(datatofilt)-(pts-ptsaway)
                    filtdata(wndw)=nanmean(datatofilt([wndw-(ptsaway)+ isEven:wndw+(ptsaway)] ));
                end
            else filtdata=datatofilt;
            end
        end
        
        %% neuralPhaseLocking - wrapper to calculate phase locking per SU
        function results = neuralPhaseLocking(obj, runData, fileNameResults, whatTorun)
            
            if nargin < 3
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            
            if whatTorun.useStableUnits
                mm = matfile('E:\Data_p\ClosedLoopDataset\spikeStimResults\validation\su_info_summary');
                results = mm.results;
                % Create of list of unstable units
                rmvUnitInd = results.userVisualValidation == 3;
                usUnit_pt = results.ptVec(rmvUnitInd);
                usUnit_ch = results.chanVec(rmvUnitInd);
                usUnit_suId = results.suVec(rmvUnitInd);
                clear results
            end
            
            %go over the patients requested
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName]);
                results(iPatient).patientName = runData(iPatient).patientName;
                results(iPatient).EXP = runData.EXP;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX;
                nStims = length(stimTimes);
                dataDuration = max(stimTimes) + obj.shortTimeRangeAfterStim*obj.samplingRate;
                
                %load spike data
                try
                    spikeData = load(runData(iPatient).spikeData);
                catch
                    disp([runData(iPatient).spikeData ' doesn''t exist, continuing']);
                    continue;
                end
                
                if isempty(spikeData.micro_channels_spike_summary.unit_list_xls)
                    continue;
                end
                
                %create a list of all the areas for which there is multi
                %unit data
                allAreas = {spikeData.micro_channels_spike_summary.unit_list_xls.Location};
                areasList = unique(allAreas);
                results(iPatient).areasList = areasList;
                results(iPatient).nStims = nStims;
                nAreas = length(areasList);
                
                channelsList = cell(1,nAreas);
                stimTriggeredRates = cell(1,nAreas);
                contTriggeredRates = cell(1,nAreas);
                ratesPerChan = cell(1,nAreas);
                nChansInAvg = zeros(1,nAreas);
                ps = zeros(1,nAreas);
                
                %get phase of slow waves
                if whatTorun.useProbe
                    macroCh = runData(iPatient).probeChan;
                    try
                        probeData = load([runData(iPatient).DataFolder,'\','CSC' ,num2str(macroCh),'.mat']);
                        probeData= probeData.data;
                        probeData(isnan(probeData)) = 0;
                    catch
                        disp('Probe wasnt found, aborting run');
                        disp([fileDir,'\','CSC' ,num2str(macroCh),'.mat doesnt exist']);
                        return;
                    end
                    
                    % sleepScoreFile
                    % load sleep scoring
                    if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                        try
                            sleepScoring = load(runData(iPatient).sleepScoringFileName);
                            sleepScoring = sleepScoring.sleep_score_vec;
                        catch
                            disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                            sleepScoring = [];
                        end
                        isSleep = sleepScoring==obj.sleepEpochs;
                    end
                    
                    % get interictal spikes
                    try
                        SpikesFileNames = [runData(iPatient).SpikesFileNames,num2str(macroCh),'.mat'];
                        ISIpeakTimes = load(SpikesFileNames);
                        ISIpeakTimes = ISIpeakTimes.peakTimes;
                    catch
                        disp([SpikesFileNames,'.mat doesn''t exist']);
                        disp('running spikes');
                        IIS_det = SpikeWaveDetector;
                        [ISIpeakTimes, ~] = IIS_det.detectTimes(probeData, true);
                    end
                    ISIpeakTimes(ISIpeakTimes < obj.minT_spike) = obj.minT_spike + 1;
                    ISIpeakTimes(ISIpeakTimes > (length(probeData)-obj.minT_spike-1)) = length(probeData)-obj.minT_spike-1;
                    
                    dataFilteredSW = bandpass(probeData, [obj.lowLimitSW, obj.highLimitSW], obj.samplingRate);
                    nanInds = isnan(dataFilteredSW);
                    dataFilteredSW(nanInds) = 0;
                    dataFilteredSW(~isSleep) = 0;
                    dataFilteredSW(ISIpeakTimes - obj.minT_spike:ISIpeakTimes + obj.minT_spike) = 0;
                    currentPhases = angle(hilbert(dataFilteredSW-mean(dataFilteredSW))); % valuse range: [-pi:pi]
                    
                end
                
                results(iPatient).currentPhases = currentPhases;
                
                % Extract stimulation phase per epoch
                for ii_e = 1:length(expData.Session_start_end_msec)
                    t0 = expData.Session_start_end_msec(ii_e,1);
                    tend = expData.Session_start_end_msec(ii_e,2);
                    
                    stimTimes_e = stimTimes(stimTimes > t0 & stimTimes < tend);
                    if length(stimTimes_e) > obj.minSpikesToIncludeUnitPhaseAnalysis &&... % minimal N spikes
                            sum(dataFilteredSW(floor(t0):floor(tend))>0) % NREM sleep
                        [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, stimTimes_e );
                        
                        results(iPatient).stimCELL.phaseLockStructEpochs{ii_e} = phaseLockStruct;
                        results(iPatient).stimCELL.IndStruct = IndStruct;
                    else
                        disp('no spikes in epoch')
                        results(iPatient).stimCELL.phaseLockStructEpochs{ii_e} = {};
                        results(iPatient).stimCELL.IndStructEpochs{ii_e} = {};
                    end
                end
                
                
                
                %go over all the areas for the patient
                for iArea = 1:nAreas
                    
                    %find all relevant units
                    unitInds = find(strcmp(allAreas,areasList{iArea}));
                    %find all the channels for that area
                    currChannels = [spikeData.micro_channels_spike_summary.unit_list_xls(unitInds).Channel];
                    channelsList{iArea} = unique(currChannels);
                    
                    for iChannel = 1:length(channelsList{iArea})
                        
                        currUinds = unitInds(currChannels==channelsList{iArea}(iChannel));
                        
                        % Remove units that weren't stable throughout the
                        % session
                        if whatTorun.useStableUnits
                            rmvInd = logical(zeros(1,length(currUinds)));
                            for ii_su = 1:length(currUinds)
                                su_ind = currUinds(ii_su);
                                ptNum = str2num(results(iPatient).patientName(2:end));
                                currCh = channelsList{iArea}(iChannel);
                                if sum((ptNum == usUnit_pt) & (usUnit_ch == currCh) & (usUnit_suId == su_ind))
                                    rmvInd(ii_su) = 1;
                                else
                                    rmvInd(ii_su) = 0;
                                end
                            end
                            
                            if sum(rmvInd)
                                currUinds(rmvInd) = [];
                                disp(sprintf('removed %d units for instability',sum(rmvInd)));
                            end
                        end
                        
                        
                        % For each single unit - look at changes in phase
                        % locking
                        for ii_su = 1:length(currUinds)
                            
                            su_ind = currUinds(ii_su);
                            spikeTimes = spikeData.micro_channels_spike_summary.spike_timestamps{su_ind};
                            
                            spikeInfo.spikeTimes_ms =  spikeData.micro_channels_spike_summary.spike_timestamps{ii_su};
                            spikeInfo.spike_shapes_mean =  spikeData.micro_channels_spike_summary.spike_shapes_mean{ii_su};
                            spikeInfo.spike_shapes_std =  spikeData.micro_channels_spike_summary.spike_shapes_std{ii_su};
                            spikeInfo.su_ind = su_ind;
                            results(iPatient).spikeCELL{iArea, iChannel,ii_su}.spikeInfo = spikeInfo;
                            
                            % phaseLockStructEpochs
                            % ii_a = 1 - baseline session
                            % ii_a = 2 - immediate - all stimulation epochs pulled together
                            % ii_a = 3 - persistent - 60sec post stimulation block
                            % ii_a = 4 - first pause epoch
                            % ii_a = 5 - all pause (1-5min)
                            % ii_a = 6 - equal to persistent condition -
                            %            last 60 sec of pause blocks
                            for ii_a = 1:6
                                events_ms_e = [];
                                
                                if ii_a == 1
                                    ii_e = 1; % baseline
                                    t0 = expData.Session_start_end_msec(ii_e,1);
                                    tend = expData.Session_start_end_msec(ii_e,2);
                                    duration_ms = tend-t0; %ms
                                    events_ms_e = spikeTimes(spikeTimes > t0 & spikeTimes < tend)';
                                    
                                    % Pull together all stimulation epochs
                                elseif ii_a == 2
                                    baselineDuration = results(iPatient).spikeCELL{iArea, iChannel,ii_su}.duration_ms(1);
                                    duration_ms = -nStims*obj.spikeMinDistIIS*2; % these vectors will be removed later
                                    for ii_e = 2:2:length(expData.Session_start_end_msec)-2
                                        t0 = expData.Session_start_end_msec(ii_e,1);
                                        tend = expData.Session_start_end_msec(ii_e,2);
                                        events_ms_e = [events_ms_e, spikeTimes(spikeTimes > t0 & spikeTimes < tend)'];
                                        duration_ms = duration_ms + (tend-t0); %ms
                                    end
                                    
                                    % Pull together all first 60Sec of 'pause' epochs
                                elseif ii_a == 3
                                    duration_ms = 0;
                                    for ii_e = 3:2:length(expData.Session_start_end_msec)
                                        t0 = expData.Session_start_end_msec(ii_e,1);
                                        tend = expData.Session_start_end_msec(ii_e,1) + obj.phaseLockingBuffer_ms;
                                        events_ms_e = [events_ms_e, spikeTimes(spikeTimes > t0 & spikeTimes < tend)'];
                                        duration_ms = duration_ms + tend-t0; %ms
                                    end
                                    
                                    % First 'pause' epoch
                                elseif ii_a == 4
                                    duration_ms = 0;
                                    for ii_e = 3
                                        t0 = expData.Session_start_end_msec(ii_e,1);
                                        tend = expData.Session_start_end_msec(ii_e,2);
                                        events_ms_e = [events_ms_e, spikeTimes(spikeTimes > t0 & spikeTimes < tend)'];
                                        duration_ms = duration_ms + tend-t0; %ms
                                    end
                                    
                                    
                                    % All 'pause' epoch (after mid-period)
                                elseif ii_a == 5
                                    duration_ms = 0;
                                    for ii_e = 3:2:length(expData.Session_start_end_msec)
                                        t0 = expData.Session_start_end_msec(ii_e,1) + obj.phaseLockingBuffer_ms;
                                        tend = expData.Session_start_end_msec(ii_e,2);
                                        events_ms_e = [events_ms_e, spikeTimes(spikeTimes > t0 & spikeTimes < tend)'];
                                        duration_ms = duration_ms + tend-t0; %ms
                                    end
                                    
                                    % All ends of 'pause' epoch (60sec before start of next session)
                                elseif ii_a == 6
                                    duration_ms = 0;
                                    for ii_e = 3:2:length(expData.Session_start_end_msec)
                                        t0 = expData.Session_start_end_msec(ii_e,2) - obj.phaseLockingBuffer_ms;
                                        tend = expData.Session_start_end_msec(ii_e,2);
                                        events_ms_e = [events_ms_e, spikeTimes(spikeTimes > t0 & spikeTimes < tend)'];
                                        duration_ms = duration_ms + tend-t0; %ms
                                    end
                                    
                                end
                                
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.duration_ms(ii_a) = duration_ms;
                                
                                
                                % remove events that are too close to
                                % stimulation times
                                rmvSpikes = sum( abs(events_ms_e(:) - repmat(stimTimes, length(events_ms_e),1)) < obj.spikeMinDistIIS,2) > 0;
                                
                                % remove events that are not during nrem
                                % sleep
                                rmvSpikes = sum( abs(events_ms_e(:) - repmat(stimTimes, length(events_ms_e),1)) < obj.spikeMinDistIIS,2) > 0;
                                rmvSpikes2 = (~ismember(floor(events_ms_e), find(isSleep)));
                                
                                disp(sprintf('%d spikes close to IIS removed from analysis (%1.1f%%)',sum(rmvSpikes),100*sum(rmvSpikes)/length(events_ms_e)))
                                disp(sprintf('%d spikes outside NREM sleep removed from analysis (%1.1f%%)',sum(rmvSpikes2),100*sum(rmvSpikes2)/length(events_ms_e)))
                                
                                events_ms_e( rmvSpikes2 | rmvSpikes') = [];
                                
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.spikeTimeForPhaseCalc{ii_a} = events_ms_e;
                                
                                
                            end % Loop over the conditions
                            
                            
                            % Calc phase-locking values
                            for ii_a = 1:6
                                
                                events_ms_e = results(iPatient).spikeCELL{iArea, iChannel,ii_su}.spikeTimeForPhaseCalc{ii_a};
                                if length(events_ms_e) > obj.minSpikesToIncludeUnitPhaseAnalysis &&... % minimal N spikes
                                        sum(dataFilteredSW(floor(t0):floor(tend))>0) % NREM sleep
                                    [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, events_ms_e );
                                    
                                    results(iPatient).spikeCELL{iArea, iChannel,ii_su}.phaseLockStructEpochs{ii_a} = phaseLockStruct;
                                    results(iPatient).spikeCELL{iArea, iChannel,ii_su}.IndStructEpochs{ii_a} = IndStruct;
                                else
                                    disp('no spikes or not sufficient long NREM epoch')
                                    results(iPatient).spikeCELL{iArea, iChannel,ii_su}.phaseLockStructEpochs{ii_a} = NaN;
                                    results(iPatient).spikeCELL{iArea, iChannel,ii_su}.IndStructEpochs{ii_a} = NaN;
                                end
                                
                            end
                            
                            % Calculate joint phase-locking when pulling
                            % together all 1-min post stim vs 1-min before
                            events_ms_baseline = [];
                            for ii_e = 2:2:length(expData.Session_start_end_msec)-1
                                tend = expData.Session_start_end_msec(ii_e,1); % start of stimulation session
                                t0 = tend - obj.phaseLockingBuffer_ms;
                                events_ms_e = spikeTimes(spikeTimes > t0 & spikeTimes < tend);
                                events_ms_baseline = [events_ms_baseline events_ms_e(:)'];
                            end
                            
                            if length(events_ms_baseline) > obj.minSpikesToIncludeUnitPhaseAnalysis &&... % minimal N spikes
                                    sum(dataFilteredSW(floor(events_ms_baseline))>0) % NREM sleep
                                [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, events_ms_baseline );
                                
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.phaseLockStructAllEpochs{1} = phaseLockStruct;
                            else
                                disp('no spikes in epoch')
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.phaseLockStructAllEpochs{1} = NaN;
                            end
                            
                            events_ms_post_stim = [];
                            for ii_e = 2:2:length(expData.Session_start_end_msec)-1
                                t0 = expData.Session_start_end_msec(ii_e,2); % end of stimulation session
                                tend = t0 + obj.phaseLockingBuffer_ms;
                                events_ms_e = spikeTimes(spikeTimes > t0 & spikeTimes < tend);
                                events_ms_post_stim = [events_ms_post_stim events_ms_e(:)'];
                            end
                            
                            if length(events_ms_post_stim) > obj.minSpikesToIncludeUnitPhaseAnalysis &&... % minimal N spikes
                                    sum(dataFilteredSW(floor(events_ms_post_stim))>0) % NREM sleep
                                [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, events_ms_post_stim );
                                
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.phaseLockStructAllEpochs{2} = phaseLockStruct;
                            else
                                disp('no spikes in epoch')
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.phaseLockStructAllEpochs{2} = NaN;
                            end
                            
                            
                        end % su
                    end % channel
                    results(iPatient).channelsList = channelsList;
                    
                    if ~isempty(fileNameResults)
                        save(fileNameResults,'results');
                    end
                    
                end % area
            end % pt
        end % func
        
        % Fit a cosine function to spike-phase distribution
        function [phaseLockStruct, structIND]= genSWSTriggeredSpikeHist(obj, currentPhases, spikeTimes_ms)
            structIND.phasesWhenSpikesOccur = 1;
            structIND.pval = 2;
            structIND.z = 3;
            structIND.fitRes = 4;
            structIND.estPhase = 5;
            structIND.r = 6;
            structIND.mu = 7;
            structIND.circ_length = 8;
            structIND.firingRate = 9;
            structIND.peakEstPhase = 10;
            
            bins = -pi:obj.binSizeRad:pi;
            midBins = bins(1:end-1) +obj.binSizeRad/2;
            
            phasesWhenSpikesOccur = currentPhases(int64(spikeTimes_ms));
            [pval, z] = circ_rtest(phasesWhenSpikesOccur);
            
            [N, edges] = histcounts(phasesWhenSpikesOccur,bins);
            fitRes = createCosFitModel(midBins',N');
            phase_r = (2*fitRes.results.a)/(2*fitRes.results.a + fitRes.results.c);
            [M, I] = max(fitRes.TC_fit);
            estPhase =  midBins(I); % all the wrap options didn't work well, so I'm looking at the max of cosine function
            
            % I'm using the positive
            % shift in order to be
            % aligned with the results of
            % circ_mean
            peakEstPhase = -estPhase; % a*cos(x+b) peaks at -b.
            phase_mu_all = circ_mean(phasesWhenSpikesOccur');
            phase_r_all = circ_r(phasesWhenSpikesOccur' );
            
            if phase_mu_all == 0
                warning('phase is exactly 0?')
            end
            phaseLockStruct{structIND.phasesWhenSpikesOccur} = phasesWhenSpikesOccur;
            phaseLockStruct{structIND.fitRes} = fitRes;
            phaseLockStruct{structIND.estPhase} = estPhase; % cos model
            phaseLockStruct{structIND.peakEstPhase} = peakEstPhase; % cos model
            phaseLockStruct{structIND.r} = phase_r; % cos model
            phaseLockStruct{structIND.mu} = phase_mu_all; % circ_mean
            phaseLockStruct{structIND.circ_length} = phase_r_all; % circ_r
            phaseLockStruct{structIND.pval} = pval;
            phaseLockStruct{structIND.z} = z;
            
            if(0)
                figure; hist(phasesWhenSpikesOccur,20); title(sprintf('mu = %2.2f, est_mu = %2.2f, p = `%2.5e',phase_mu_all, estPhase, pval))
            end
        end
        
        % BootStraping for phase locking calc - based on pre-run of neuralPhaseLocking
        function results = neuralPhaseLocking_BS(obj, NeuralPhaseLockingOutputFile, outputFile, whatTorun)
            mm = matfile(NeuralPhaseLockingOutputFile);
            results = mm.results;
            ptV = []; areaV= []; chV = []; cnt = 0;
            isMTL = []; isFrontal = []; isProbeHem = [];
            BS_N = obj.phaseLocking_BS_N ;
            
            for iPatient = 1:length(results)
                
                pt = str2num(results(iPatient).patientName(2:end));
                ptResults = results(iPatient);
                spikeCell = ptResults.spikeCELL;
                currentPhases = results(iPatient).currentPhases; % PROBE phase
                
                % Bootstrap values of each phase-locking calculation
                % to make sure it's independant of specific spike selection
                % Here we bootstrap condition 3 and 6 -
                % Prolonged effect and decay of effect
                ses_1 = 3;
                ses_2 = 6;
                
                [areaN, Nchannel , Nsu] = size(ptResults.spikeCELL);
                for aa = 1:areaN
                    
                    area = classifyArea(results(iPatient).areasList{aa});
                    
                    for cc = 1:Nchannel
                        for nn = 1:Nsu
                            
                            CELL = spikeCell{aa,cc,nn};
                            if isempty(CELL); continue; end
                            if length(CELL.phaseLockStructEpochs{ses_1}) == 1; continue; end % spike # too low / out of NREM sleep
                            
                            cnt = cnt + 1;
                            
                            ptV = [ptV pt];
                            areaV = [areaV aa];
                            chV = [chV cc];
                            isMTL = [isMTL area.isMTL];
                            isFrontal = [isFrontal area.isFrontal];
                            
                            % first entry is ground-truth
                            [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, CELL.spikeTimeForPhaseCalc{ses_1} );
                            BS_phaseLock_CELL{1,cnt, 1} = phaseLockStruct;
                            BS_IndStruct_CELL{1,cnt, 1} = IndStruct;
                            pval_bs(1,cnt,1) = phaseLockStruct{IndStruct.pval};
                            
                            
                            [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, CELL.spikeTimeForPhaseCalc{ses_2} );
                            BS_phaseLock_CELL{2,cnt, 1} = phaseLockStruct;
                            BS_IndStruct_CELL{2,cnt, 1} = IndStruct;
                            pval_bs(2,cnt,1) = phaseLockStruct{IndStruct.pval};
                            
                            for bs_i = 1:BS_N
                                
                                
                                Nmatch = min(length( CELL.spikeTimeForPhaseCalc{ses_1}),...
                                    length(CELL.spikeTimeForPhaseCalc{ses_2}));
                                Nmatch_BS = floor(0.9 * Nmatch); % Shuffle 90% of spikes from the lower number of spikes
                                x = randi(length( CELL.spikeTimeForPhaseCalc{ses_1}), 1, Nmatch_BS);
                                CELL.BS_spikeTimeForPhaseCalc{1} = CELL.spikeTimeForPhaseCalc{ses_1}(x);
                                x = randi(length( CELL.spikeTimeForPhaseCalc{ses_2}), 1, Nmatch_BS);
                                CELL.BS_spikeTimeForPhaseCalc{2} = CELL.spikeTimeForPhaseCalc{ses_2}(x);
                                
                                [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, CELL.BS_spikeTimeForPhaseCalc{1} );
                                BS_phaseLock_CELL{1,cnt,bs_i+1} = phaseLockStruct;
                                pval_bs(1,cnt,bs_i+1) = phaseLockStruct{IndStruct.pval};
                                
                                [phaseLockStruct, IndStruct] = obj.genSWSTriggeredSpikeHist(currentPhases, CELL.BS_spikeTimeForPhaseCalc{2} );
                                BS_phaseLock_CELL{2,cnt,bs_i+1} = phaseLockStruct;
                                pval_bs(2,cnt,bs_i+1) = phaseLockStruct{IndStruct.pval};
                                
                                
                            end % BS
                            
                            
                            % sanity check
                            for ii = 1:2
                                gt(ii) = BS_phaseLock_CELL{ii,cnt,1}{IndStruct.r};
                            end
                            gt_diff(cnt) = (gt(1)-gt(2))./(gt(1)+gt(2));
                            for ii_bs = 1:BS_N
                                for ii = 1:2
                                    r_bs(ii_bs,ii) = BS_phaseLock_CELL{ii,cnt,ii_bs}{IndStruct.r};
                                end
                                diff_bs(cnt,ii_bs) = (r_bs(ii_bs,1)-r_bs(ii_bs,2))./(r_bs(ii_bs,1)+r_bs(ii_bs,2));
                                
                            end
                            
                            
                        end % Nsu
                        
                    end % Nch
                end % Na
                
                % save sub-files
                midFilename = sprintf('shuffle_phaseLock_%d',ptV(end));
                folderName = 'E:\Data_p\ClosedLoopDataset\SWTriggeredSpikeHistResults\SHUFFLE';
                save(fullfile(folderName,midFilename),'BS_phaseLock_CELL', 'diff_bs','gt_diff','ptV','chV','areaV','-v7.3')
                BS_phaseLock_CELL =[];
            end % pt
            
            
            bs_results.ptV = ptV;
            bs_results.areaV = areaV;
            bs_results.chV = chV;
            bs_results.isMTL = isMTL;
            bs_results.isFrontal = isFrontal;
            bs_results.IndStruct = IndStruct;
            bs_results.diff_bs = diff_bs;
            bs_results.gt_diff = gt_diff;
            bs_results.pval_bs = pval_bs;
            
            ind =  find(sum(pval_bs(1:2,:,1) < 0.05) == 2);
            B = gt_diff(ind);
            [h,p] = signrank(gt_diff(ind));
            A = diff_bs(ind,:);
            [p,h] = signrank(A(:));
            figure;
            histogram(gt_diff(ind),'Normalization','Probability');
            hold all;
            h2 = histogram(A(:),'Normalizatio','Probability');
            
            [h,p ] = kstest2(B(:),A(:));
            
            % full dataset is huge
            save(outputFile, 'bs_results','-v7.3')
            
            %              bs_results.BS_phaseLock_CELL = BS_phaseLock_CELL;
            %              outputFile =  'E:\Data_p\ClosedLoopDataset\SWTriggeredSpikeHistResults\SWTriggeredSpikeHist_allSes_AllPts_stableUnits_BS_1000_FULL.mat';
            %              save(outputFile, 'bs_results','-v7.3')
            
        end % function
        
        
        function neuralPhaseLocking_BS_PLOT(obj, NeuralPhaseLocking_BS_OutputFile,  whatTorun)
            
            mm = matfile(NeuralPhaseLocking_BS_OutputFile);
            bs_results = mm.bs_results;
            gt_diff = bs_results.gt_diff;
            diff_bs = bs_results.diff_bs;
            
            M1 = median(diff_bs(:));
            M2 = median(gt_diff(:));
            
            ind =  find(sum(pval_bs(1:2,:,1) < 0.05) == 2);
            B = gt_diff(ind);
            [h,p] = signrank(gt_diff(ind));
            A = diff_bs(ind,:);
            [p,h] = signrank(A(:));
            figure; histogram(gt_diff(ind))
            hold all;
            h2 = histogram(A(:));
            
            
        end
        
        function results = stimulationEffectCortical(obj,runData,fileNameResults, whatToRun)
            %
            % In cortical channels the events that are of interest are: slow waves, spindles, and couplets of slow waves
            % followed by spindles (between minDistSWSpindle to
            % maxDistSWSpindle ms after a slow wave).
            % The method calculates the rates of these events around stimulations and around control points. The controls
            % are random point in the quiet epochs between stimulations. The number of controls is set by the property
            % nControls. The rate function is created by using movsum on the events scatter, using a
            % window of avgFigBinning ms (by default 100 ms). The duration of the window around the stimulations is defined
            % by timeToCheckBeforeSW and timeToCheckAfterSW for slow waves, timeToCheckBeforeSp and timeToCheckAfterSp for
            % spindles.mid
            
            %
            % In addition to the rate function, the method calculates: a. the probability of an event post stimulation:
            % #stimulations after which the event occurred / #total number of stimulations, calculated similarly for control
            % points, b. the event rate post stimulation: number of times the event occurs in the post stimulation epoch /
            % total duration of post stimulation epoch (in minutes). The baseline (“control”) for this variable is the
            % number of times the event occurred in the quiet epoch divided by that epoch’s duration. The post stimulation
            % window for calculating both the probability and the rate is set by timeToCheckAfterSW for slow waves and
            % timeToCheckBeforeSp.
            %
            % The input runData is a struct in the length of number of patients (for which the analysis is required). In
            % addition it receives the input parameter fileNameResults which includes the file name into which the results
            % will be saved (optional).
            % Each element (=patient) in runData should include the fields:
            % patientName
            % corticalChans – list of channel indices for which to perform the analysis.
            % ExpDataFileName – name (including path) of the EXP_DATA for the patient.
            % macroMontageFileName - the file name (including path) of the macromontage.
            % DataFolder – The folder in which the raw data files are saved (the method assumes the prefix for the files is
            % CSC, can be changed by the property dataFilePrefix).
            % SpindlesFileNames - name (including path) of the spindle mat files in which the spindle times for the macro
            % channels are saved (the method assumes the name of the file is SpindlesFileNames <#channel index>).
            % SWStaresinaFileName - name (including path) of the slow wave mat files in which the slow waves times for the
            % macro channels are saved (the method assumes the name of the file is SWStaresinaFileName <#channel index>).
            % sleepScoringFileName – file name (including path) of the sleep scoring mat file.
            %
            % The output struct results includes all the results of the analysis, which can then be plotted using
            % plotStimEffectCortical. The output struct is a struct with the length of the number of patients (=the length
            % of runData), where each element includes:
            % patientName
            % resultsPerChan – a struct in the length of the number of channels required for the analysis per the patient.
            % Each element in resultsPerChan includes the fields:
            % channelNum
            % nStims – The number of stimulations.
            % baselineSWRate – rate of slow waves in the epoch before the stimulations (#sw/second).
            % baselineSpRate – rate of spindles in the epoch before the stimulations (#spindles/second).
            % baselineSWSpRate – rate of slow waves-spindles couples in the epoch before the stimulations
            % (#sw-spindles/second).
            % timesAroundStimSW – a cell array with two elements: 1. A matrix where each row is the sw rate function around
            % a stimulation, 2. A matrix where each row is the sw rate function around a control.
            % timesAroundStimSp – a cell array with two elements: 1. A matrix where each row is the spindle rate function
            % around a stimulation, 2. A matrix where each row is the spindle rate function around a control.
            % timesAroundStimSWSp - a cell array with two elements: 1. A matrix where each row is the SW-spindle couples
            % rate function around a stimulation, 2. A matrix where each row is the SW-spindle couples rate function around
            % a control.
            % dataAroundStims - a cell array with two elements: 1. A matrix where each row is the raw data around a
            % stimulation, 2. A matrix where each row is the raw data around a control.
            % probsStimSW – an array with two elements: 1. The probability of slow waves post stimulation, 2. The
            % probability of slow waves post control.
            % probsStimSp - an array with two elements: 1. The probability of spindles post stimulation, 2. The probability
            % of spindles post control.
            % probsStimSWSp - an array with two elements: 1. The probability of SW-spindles post stimulation, 2. The
            % probability of SW-spindles post control.
            % eventsPerMinSW – an array with two elements: 1. The event rate of slow waves post stimulation, 2. The
            % event rate of slow waves post control.
            % eventsPerMinSp - an array with two elements: 1. The event rate of spindles post stimulation, 2. The
            % event rate of spindles post control.
            % eventsPerMinSWSp - an array with two elements: 1. The event rate of SW-spindles post stimulation, 2.
            % The event rate of SW-spindles post control.
            
            
            if nargin < 3
                fileNameResults = '';
            end
            
            minDistSWSpindle = obj.minDistSWSpindle*obj.samplingRate/1000;
            maxDistSWSpindle = obj.maxDistSWSpindle*obj.samplingRate/1000;
            avgFigBinning = obj.avgFigBinning*obj.samplingRate/1000;
            timeToCheckAfterSW = obj.timeToCheckAfterSW*obj.samplingRate/1000;
            timeToCheckAfterSp = obj.timeToCheckAfterSp*obj.samplingRate/1000;
            timeToCheckMax = max(timeToCheckAfterSW,timeToCheckAfterSp);
            timeToCheckBeforeSW = obj.timeToCheckBeforeSW*obj.samplingRate/1000;
            timeToCheckBeforeSp = obj.timeToCheckBeforeSp*obj.samplingRate/1000;
            midTimeRangeAfterStim = obj.midTimeRangeAfterStim*obj.samplingRate;
            
            %go over the required patients
            nPatients = length(runData);
            for iPatient = 1:nPatients
                clear results
                disp(['Patient ',runData(iPatient).patientName,' ',num2str(iPatient),'/',num2str(nPatients)]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX; stimTimes = floor(stimTimes(:))';
                
                
                %load macro montage for area name
                macroMontage = load(runData(iPatient).macroMontageFileName);
                macroMontage = macroMontage.MacroMontage;
                
                %load sleep scoring
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                
                % dataDuration = stimTimes(end)+timeToCheckMax;
                % extend duration to include post-stim session
                dataDuration = floor(expData.Session_start_end_msec(end,2));
                disp( sprintf('removing %d stimulations detected after end of session',sum(stimTimes > dataDuration)))
                stimTimes(stimTimes > dataDuration) = [];
                
                %check only stim times during NREM
                if ~isempty(sleepScoring)
                    stimTimes = stimTimes(sleepScoring(round(stimTimes))==1);
                end
                results(iPatient).stimTimes_NREM = stimTimes;
                
                nStims = length(stimTimes);
                
                resultsPerChan = [];
                nChans = length(runData(iPatient).corticalChans);
                
                
                % find the quiet epochs times between stimulations (breaks
                % which are beyond the midTimeRangeAfterStim)
                quietEpochsBetween = zeros(1,dataDuration);
                stimDiffs = diff(stimTimes);
                stimIndsWithMidPauseAfter = [find(stimDiffs >= midTimeRangeAfterStim) length(stimTimes)];
                stimIndsWithMidPauseAfterTimes = stimTimes(stimIndsWithMidPauseAfter); % start of pause block
                
                stimIndsAfterPauseTimes = stimTimes([1,stimIndsWithMidPauseAfter(1:end-1)+1]); % end of pause block
              
                for iStim = 1:length(stimIndsWithMidPauseAfter)-1
                    quietEpochsBetween(stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim+1:stimTimes(stimIndsWithMidPauseAfter(iStim)+1)-timeToCheckMax-obj.timeBeforeStim) = 1;
                end
                quietEpochsBetweenTimes = find(quietEpochsBetween);
                
                % March 21, add midRange estimate - times that are
                % immediatly after end of stimulation block
                quietEpochsMidRange = zeros(1,dataDuration);
                for iStim = 1:length(stimIndsWithMidPauseAfter)-1
                    quietEpochsMidRange(stimIndsWithMidPauseAfterTimes(iStim)+timeToCheckMax+1:stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim) = 1;
                end
                quietEpochsMidRangeTimes = find(quietEpochsMidRange);

                firstStim = stimTimes(1)-obj.timeBeforeStim;
                durationBeforeStim = sum(sleepScoring(1:firstStim)==obj.sleepEpochs)/obj.samplingRate; %seconds
                
                % March 2022 - add an *equal* time for CONTROL of midRange
                % estimate - the last chunck of the quiet epochs.
                quietEpochsMidRange_CONTROL = zeros(1,dataDuration);
                for iStim = 1:length(stimIndsAfterPauseTimes)-1
                    quietEpochsMidRange_CONTROL(stimIndsAfterPauseTimes(iStim)-midTimeRangeAfterStim+1:stimIndsAfterPauseTimes(iStim)-timeToCheckMax) = 1;
                end
                quietEpochsMidRange_CONTROLTimes = find(quietEpochsMidRange_CONTROL);

                %go over the channels of interest
                for iChan = 1:nChans
                    
                    currChan = runData(iPatient).corticalChans(iChan);
                    resultsPerChan(iChan).channelNum = currChan;
                    resultsPerChan(iChan).nStims = nStims;
                    
                    currArea = macroMontage(currChan).Area;
                    resultsPerChan(iChan).area = currArea;
                    
                    %load data (to show average data around stim)
                    try
                        data = [runData(iPatient).DataFolder '\',obj.dataFilePrefix, num2str(currChan) '.mat'];
                        data = load(data);
                        data = data.data;
                    catch
                        disp([runData(iPatient).DataFolder '\',obj.dataFilePrefix, num2str(currChan) '.mat doesn''t exist']);
                        continue;
                    end
                    
                    
                    %load spindles
                    if isfield(runData(iPatient), 'SpindlesFileNames') && ~isempty(runData(iPatient).SpindlesFileNames)
                        try
                            spindlesTimes = [runData(iPatient).SpindlesFileNames num2str(currChan) '.mat'];
                            spindlesTimes = load(spindlesTimes);
                            spindlesTimes = spindlesTimes.spindlesTimes;
                        catch
                            disp([runData(iPatient).SpindlesFileNames num2str(currChan) '.mat doesn''t exist']);
                            spindlesTimes = [];
                        end
                    end
                    
                    %load slow waves
                    if isfield(runData(iPatient), 'SWStaresinaFileName') && ~isempty(runData(iPatient).SWStaresinaFileName)
                        try
                            slowWavesTimes = [runData(iPatient).SWStaresinaFileName num2str(currChan) '.mat'];
                            slowWavesTimes = load(slowWavesTimes);
                            slowWavesTimes = slowWavesTimes.slowWavesTimes;
                        catch
                            disp([runData(iPatient).SWStaresinaFileName num2str(currChan) '.mat doesn''t exist']);
                            slowWavesTimes = [];
                        end
                    end
                    
                    %load interictal spikes
                    if isfield(runData(iPatient), 'SpikesFileNames') && ~isempty(runData(iPatient).SpikesFileNames)
                        try
                            interictalSpikeTimes = [runData(iPatient).SpikesFileNames num2str(currChan) '.mat'];
                            interictalSpikeTimes = load(interictalSpikeTimes);
                            interictalSpikeTimes = interictalSpikeTimes.peakTimes;
                        catch
                            disp([runData(iPatient).SpikesFileNames num2str(currChan) '.mat doesn''t exist']);
                            interictalSpikeTimes = [];
                        end
                    end
                    
                    interictalSpikeTimes = interictalSpikeTimes(interictalSpikeTimes<dataDuration);
                    spindlesTimes = spindlesTimes(spindlesTimes<dataDuration);
                    slowWavesTimes = slowWavesTimes(slowWavesTimes<dataDuration);
                    
                    
                    spIndices = zeros(1,dataDuration+maxDistSWSpindle);
                    spIndices(spindlesTimes) = 1;
                    
                    swSpTimes = [];
                    
                    %first find coupled slow wave - spindles
                    for iSW = 1:length(slowWavesTimes)
                        if sum(spIndices(slowWavesTimes(iSW)+minDistSWSpindle:slowWavesTimes(iSW)+maxDistSWSpindle))>0
                            swSpTimes(end+1) = slowWavesTimes(iSW);
                        end
                    end
                    
                    %find baseline rates (rates during the epoch before the
                    %stimulations)
                    baselineIISRate = length(interictalSpikeTimes(interictalSpikeTimes<firstStim))/durationBeforeStim;
                    baselineSWRate = length(slowWavesTimes(slowWavesTimes<firstStim))/durationBeforeStim;
                    baselineSpRate = length(spindlesTimes(spindlesTimes<firstStim))/durationBeforeStim;
                    baselineSWSpRate = length(swSpTimes(swSpTimes<firstStim))/durationBeforeStim;
                    
                    resultsPerChan(iChan).baselineIISRate = baselineIISRate;
                    resultsPerChan(iChan).baselineSWRate = baselineSWRate;
                    resultsPerChan(iChan).baselineSpRate = baselineSpRate;
                    resultsPerChan(iChan).baselineSWSpRate = baselineSWSpRate;
                    
                    %build stimulation triggered rates for all three events
                    %types (sw, spindles, sw-sp)
                    
                    spRate = movsum(spIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    swIndices = zeros(1,dataDuration);
                    swIndices(slowWavesTimes) = 1;
                    swRate = movsum(swIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    IIS_Indices = zeros(1,dataDuration);
                    IIS_Indices(interictalSpikeTimes) = 1;
                    IIS_Rate = movsum(IIS_Indices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    swSpIndices = zeros(1,dataDuration);
                    swSpIndices(swSpTimes) = 1;
                    swSpRate = movsum(swSpIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    
                    
                    
                    timesAroundStimIIS = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    probsStimIIS = [0 0];
                    eventsNIIS = [0 0];
                    
                    timesAroundStimSW = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    probsStimSW = [0 0];
                    eventsNSW = [0 0];
                    
                    timesAroundStimSp = nan(nStims,timeToCheckBeforeSp+timeToCheckAfterSp+1);
                    probsStimSp = [0 0];
                    eventsNSp = [0 0];
                    
                    timesAroundStimSWSp = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    probsStimSWSp = [0 0];
                    eventsNSWSp = [0 0];
                    
                    dataAroundStims = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    
                    %go over the stimulations
                    for iStim = 1:nStims
                        
                        %get raw data around stims
                        dataAroundStims(iStim,:) = data(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        
                        %get the IIS rate function around stimulations and
                        %update the prob and rates
                        timesAroundStimIIS(iStim,:) = IIS_Rate(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        currNIISStim = sum(IIS_Indices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSW));
                        if currNIISStim>0
                            probsStimIIS(1) = probsStimIIS(1)+1;
                            eventsNIIS(1) = eventsNIIS(1)+currNIISStim;
                        end
                        
                        %get the sw rate function around stimulations and
                        %update the prob and rates
                        timesAroundStimSW(iStim,:) = swRate(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        currNSWStim = sum(swIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSW));
                        if currNSWStim>0
                            probsStimSW(1) = probsStimSW(1)+1;
                            eventsNSW(1) = eventsNSW(1)+currNSWStim;
                        end
                        
                        %get the spindle rate function around stimulations and
                        %update the prob and rates
                        timesAroundStimSp(iStim,:) = spRate(stimTimes(iStim)-timeToCheckBeforeSp:stimTimes(iStim)+timeToCheckAfterSp);
                        currNSpStim = sum(spIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSp));
                        if currNSpStim>0
                            probsStimSp(1) = probsStimSp(1)+1;
                            eventsNSp(1) = eventsNSp(1)+currNSpStim;
                        end
                        
                        %get the SW-spindle rate function around stimulations and
                        %update the prob and rates
                        timesAroundStimSWSp(iStim,:) = swSpRate(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        currNSWSpStim = sum(swSpIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSW));
                        if currNSWSpStim>0
                            probsStimSWSp(1) = probsStimSWSp(1)+1;
                            eventsNSWSp(1) = eventsNSWSp(1)+currNSWSpStim;
                        end
                    end
                    
                    %calculate the prob - number of times the event
                    %occurred post stim divided by the total number of
                    %stimulations
                    probsStimIIS(1) = probsStimIIS(1)/nStims;
                    probsStimSW(1) = probsStimSW(1)/nStims;
                    probsStimSp(1) = probsStimSp(1)/nStims;
                    probsStimSWSp(1) = probsStimSWSp(1)/nStims;
                    
                    %calculate the event rate - rate of occurrence of the
                    %event post stim divided by the duration of post stim
                    %epochs
                    eventsPerMinIIS(1) =  eventsNIIS(1)/(nStims*obj.timeToCheckAfterSW/60000); %events per min
                    eventsPerMinSW(1) = eventsNSW(1)/(nStims*obj.timeToCheckAfterSW/60000); %events per min
                    eventsPerMinSp(1) = eventsNSp(1)/(nStims*obj.timeToCheckAfterSp/60000);
                    eventsPerMinSWSp(1) = eventsNSWSp(1)/(nStims*obj.timeToCheckAfterSW/60000);
                    
                    %create controls
                    timesAroundStimIISControl = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    timesAroundStimSWControl = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    timesAroundStimSpControl = nan(obj.nControls,timeToCheckBeforeSp+timeToCheckAfterSp+1);
                    timesAroundStimSWSpControl = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    dataAroundConts = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    
                    
                    if whatToRun.controlTimestamps
                        try
                            stim_control_sets = load(runData(iPatient).controlTimestampsFilename,'stim_control_sets');
                            controlTimestamps = stim_control_sets.stim_control_sets.shamTimestamps;
                            controlTimestamps(controlTimestamps > (length(IIS_Rate) - timeToCheckAfterSW)) = [];
                            obj.nControls = length(controlTimestamps);
                        catch
                            warning('control timestamps not found')
                        end
                    end
                    
                    for iControl = 1:obj.nControls
                        
                        if whatToRun.controlTimestamps
                            % controls are matched to stim distance to
                            % slow-waves in the quiet epochs
                            controlInd = controlTimestamps(iControl);
                        else
                            %controls are random indices in the quiet epochs
                            controlInd = quietEpochsBetweenTimes(randi(length(quietEpochsBetweenTimes)));
                        end
                        controlTimestamps(iControl) = controlInd;
                        
                        %get raw data around controls
                        dataAroundConts(iControl,:) = data(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        
                        %get the IIS rate function around controls and
                        %update the prob and rates
                        timesAroundStimIISControl(iControl,:) = IIS_Rate(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        currNIISCont = sum(IIS_Indices(controlInd:controlInd+timeToCheckAfterSW));
                        if currNIISCont>0
                            probsStimIIS(2) = probsStimIIS(2)+1;
                            eventsNIIS(2) = eventsNIIS(2)+currNIISCont; % Mar 2021
                        end
                        
                        %get the sw rate function around controls and
                        %update the prob and rates
                        timesAroundStimSWControl(iControl,:) = swRate(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        currNSWCont = sum(swIndices(controlInd:controlInd+timeToCheckAfterSW));
                        if currNSWCont>0
                            probsStimSW(2) = probsStimSW(2)+1;
                            eventsNSW(2) = eventsNSW(2)+currNSWCont; % Mar 2021
                        end
                        
                        %get the spindle rate function around controls and
                        %update the prob and rates
                        timesAroundStimSpControl(iControl,:) = spRate(controlInd-timeToCheckBeforeSp:controlInd+timeToCheckAfterSp);
                        currNSpCont = sum(spIndices(controlInd:controlInd+timeToCheckAfterSp));
                        if currNSpCont>0
                            probsStimSp(2) = probsStimSp(2)+1;
                            eventsNSp(2) = eventsNSp(2)+currNSpCont;
                        end
                        
                        %get the SW-spindle rate function around controls and
                        %update the prob and rates
                        timesAroundStimSWSpControl(iControl,:) = swSpRate(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        currNSWSpCont = sum(swSpIndices(controlInd:controlInd+timeToCheckAfterSW));
                        if currNSWSpCont>0
                            probsStimSWSp(2) = probsStimSWSp(2)+1;
                            eventsNSWSp(2) = eventsNSWSp(1)+currNSWSpCont;
                        end
                    end
                    
                    %calculate the prob - number of times the event
                    %occurred post control divided by the total number of
                    %controls
                    probsStimIIS(2) = probsStimIIS(2)/obj.nControls;
                    probsStimSW(2) = probsStimSW(2)/obj.nControls;
                    probsStimSp(2) = probsStimSp(2)/obj.nControls;
                    probsStimSWSp(2) = probsStimSWSp(2)/obj.nControls;
                    
                    %calculate the event rate - rate of occurrence of the
                    %control event post stim divided by the duration of post stim
                    %epochs
                    eventsPerMinIIS(2) =  eventsNIIS(2)/(obj.nControls*obj.timeToCheckAfterSW/60000); %events per min
                    eventsPerMinSW(2) = eventsNSW(2)/(obj.nControls*obj.timeToCheckAfterSW/60000); %events per min
                    eventsPerMinSp(2) = eventsNSp(2)/(obj.nControls*obj.timeToCheckAfterSp/60000);
                    eventsPerMinSWSp(2) = eventsNSWSp(2)/(obj.nControls*obj.timeToCheckAfterSW/60000);
                    
                    
                    % Orig - "control" baseline for the event rate is the event rate over entire quiet epoch
                    % FIX - Nov 2022 - make the control shorter to match
                    % the intermidiate condition below
%                 
                    %event rate (per min) over entire quiet epoch blocks
                    %(obj.midTimeRangeAfterStim sec after each stimulation block)
%                     eventsPerMinIIS(3) =  sum(IIS_Indices&quietEpochsBetween)/(sum(quietEpochsBetween)/60000); %events per min
%                     eventsPerMinSW(3) = sum(swIndices&quietEpochsBetween)/(sum(quietEpochsBetween)/60000); %events per min
%                     eventsPerMinSp(3) = sum(spIndices(1:dataDuration)&quietEpochsBetween)/(sum(quietEpochsBetween)/60000);
%                     eventsPerMinSWSp(3) = sum(swSpIndices&quietEpochsBetween)/(sum(quietEpochsBetween)/60000);

                    eventsPerMinIIS(3) = sum(IIS_Indices&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000); %events per min
                    eventsPerMinSW(3) = sum(swIndices&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000); %events per min
                    eventsPerMinSp(3) = sum(spIndices(1:dataDuration)&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000);
                    eventsPerMinSWSp(3) = sum(swSpIndices&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000);
                  

                    % event rate (per min) over short period immediatly following stimulation block)
                    eventsPerMinIIS(4) =  sum(IIS_Indices&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000); %events per min
                    eventsPerMinSW(4) = sum(swIndices&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000); %events per min
                    eventsPerMinSp(4) = sum(spIndices(1:dataDuration)&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000);
                    eventsPerMinSWSp(4) = sum(swSpIndices&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000);
                    
                    %get std of data (for clipping the figures if necessary)
                    resultsPerChan(iChan).dataSTD = nanstd(data(sleepScoring==1));
                    
                    resultsPerChan(iChan).timesAroundStimSW = {timesAroundStimSW, timesAroundStimSWControl};
                    resultsPerChan(iChan).timesAroundStimSp = {timesAroundStimSp, timesAroundStimSpControl};
                    resultsPerChan(iChan).timesAroundStimSWSp = {timesAroundStimSWSp, timesAroundStimSWSpControl};
                    resultsPerChan(iChan).dataAroundStims = {dataAroundStims, dataAroundConts};
                    resultsPerChan(iChan).probsStimIIS = probsStimIIS;
                    resultsPerChan(iChan).probsStimSW = probsStimSW;
                    resultsPerChan(iChan).probsStimSp = probsStimSp;
                    resultsPerChan(iChan).probsStimSWSp = probsStimSWSp;
                    resultsPerChan(iChan).eventsPerMinIIS = eventsPerMinIIS;
                    resultsPerChan(iChan).eventsPerMinSW = eventsPerMinSW;
                    resultsPerChan(iChan).eventsPerMinSp = eventsPerMinSp;
                    resultsPerChan(iChan).eventsPerMinSWSp = eventsPerMinSWSp;
                    
                    
                    % Adding spectral analysis triggered by stimulation
                    % times
                    if isfield(whatToRun,'TFR') && whatToRun.TFR
                        
                        timeAfterEventStim = obj.timeBeforeAfterEventStim*obj.samplingRate;
                        TimeAfterEventROIStimSp = 2*obj.samplingRate;
                        indsForROIStimSpTime = 1:TimeAfterEventROIStimSp;
                        
                        pacCalc = PACCalculator;
                        pacCalc.timeBeforeAfterEvent = obj.timeBeforeAfterEventStim; %seconds
                        pacCalc.timeForBaseline = obj.timeForBaselineStim; %seconds
                        pacCalc.minNCycles = obj.minNCyclesStim;
                        pacCalc.freqRange = [5:30]; % Hz - to include slow wave range
                        
                        indsForROIStimSpFreq = find(pacCalc.freqRange==obj.ROIFreqRegionSWSp(1)):find(pacCalc.freqRange==obj.ROIFreqRegionSWSp(2));
                        controlFreq = [20,27];
                        indsForROIStimSpFreqControl = find(pacCalc.freqRange==controlFreq(1)):find(obj.freqRangeSWSp==controlFreq(2));
                        
                        % meanSpec - average of the normalized spectograms
                        % meanEpochs - average of the event epochs
                        % allSpec - returns all the single spectograms
                        % nEpochs - number of epochs in the calculation
                        toPlot = false;
                        stimV = true; % limits analysis to post-stim
                        specFileName = [];
                        [meanSpec, meanEpochs, stdEpochs, allSpecs, nEpochs] = pacCalc.plotAvgSpecDiff(data, stimTimes, toPlot, stimV, specFileName);
                        
                        % use control timestamps used for probability calc,
                        % to create SHAM TFRs
                        [SHAM_meanSpec, SHAM_meanEpochs, SHAM_stdEpochs, SHAM_allSpecs, SHAM_nEpochs] = pacCalc.plotAvgSpecDiff(data, controlTimestamps, toPlot, stimV, specFileName);
                        
                        
                        resultsPerChan(iChan).meanSpecs{1} = meanSpec;
                        resultsPerChan(iChan).stdEpochs{1} = stdEpochs;
                        resultsPerChan(iChan).meanEpochs{1} = meanEpochs;
                        % resultsPerChan(iChan).allSpecs{1} = allSpecs;
                        % MGS - this caused the files to be huge
                        resultsPerChan(iChan).nEpochs{1} = nEpochs;
                        resultsPerChan(iChan).meanSpecs{2} = SHAM_meanSpec;
                        resultsPerChan(iChan).stdEpochs{2} = SHAM_stdEpochs;
                        resultsPerChan(iChan).meanEpochs{2} = SHAM_meanEpochs;
                        % resultsPerChan(iChan).allSpecs{2} = SHAM_allSpecs;
                        resultsPerChan(iChan).nEpochs{2} = SHAM_nEpochs;
                        
                        if ~isempty(allSpecs)
                            resultsPerChan(iChan).ROIStimSp{1} = nanmean(nanmean(allSpecs(:,indsForROIStimSpFreq,indsForROIStimSpTime),2),3);
                            resultsPerChan(iChan).ROIcontrol{1} = nanmean(nanmean(allSpecs(:,indsForROIStimSpFreqControl,indsForROIStimSpTime),2),3);
                            [~, resultsPerChan(iChan).ROIdiff_p(1)] = ttest(resultsPerChan(iChan).ROIStimSp{1}, resultsPerChan(iChan).ROIcontrol{1});
                        else
                            resultsPerChan(iChan).ROIStimSp(1) = NaN;
                            resultsPerChan(iChan).ROIcontrol(1) = NaN;
                            resultsPerChan(iChan).ROIdiff_p(1) = NaN;
                        end
                        if ~isempty(SHAM_allSpecs)
                            resultsPerChan(iChan).ROIStimSp{2} = nanmean(nanmean(SHAM_allSpecs(:,indsForROIStimSpFreq,indsForROIStimSpTime),2),3);
                            resultsPerChan(iChan).ROIcontrol{2} = nanmean(nanmean(SHAM_allSpecs(:,indsForROIStimSpFreqControl,indsForROIStimSpTime),2),3);
                            [~, resultsPerChan(iChan).ROIdiff_p(2) ] = ttest(resultsPerChan(iChan).ROIStimSp{2}, resultsPerChan(iChan).ROIcontrol{2});
                            [~, resultsPerChan(iChan).spROI_diff_SHAM_STIM_p ] = ttest2(resultsPerChan(iChan).ROIStimSp{1}, resultsPerChan(iChan).ROIStimSp{2});
                        else
                            resultsPerChan(iChan).ROIStimSp{2} = NaN;
                            resultsPerChan(iChan).ROIcontrol{2} = NaN;
                            resultsPerChan(iChan).ROIdiff_p(2) = NaN;
                        end
                        
                    end
                    
                    
                end
                
                results(iPatient).resultsPerChan = resultsPerChan;
                
                if ~isempty(fileNameResults)
                    [filepath,name,ext] = fileparts(fileNameResults);
                    sep_filename = fullfile(filepath, 'perPt', sprintf('%s_%s',name,runData(iPatient).ExpDataFileName(end-20:end-12)));
                    disp(['saving ',sep_filename])
                    results_per_pt = results(iPatient);
                    save(sep_filename,'results_per_pt','-v7.3');
                end
                
                clear results
                
            end
            
            %             if ~isempty(fileNameResults)
            %                 save(fileNameResults,'results','-v7.3');
            %             end
        end
        
        function results = stimulationEffectMTL(obj,runData,fileNameResults, whatToRun)
            %
            % In MTL channels the events that are of interest are: ripples, slow waves, spindles, and triplets of
            % ripples-slow waves-spindles (ripple-SW couplet is defined as a a distance between ripple to SW of
            % minDistRippleSW to maxDistRippleSW, and distance of SW to spindles between
            % minDistSWSpindle to maxDistSWSpindle ms after a slow wave).
            % The method calculates the rates of these events around stimulations and around control points. The controls
            % are random point in the quiet epochs between stimulations. The number of controls is set by the property
            % nControls. The rate function is created by using movsum on the events scatter, using a window of
            % avgFigBinning ms. The duration of the window around the stimulations is defined by timeToCheckBeforeSW to
            % timeToCheckAfterSW for slow waves, timeToCheckBeforeSp to timeToCheckAfterSp for spindles and
            % timeToCheckBeforeRip to timeToCheckAfterRip for ripples.
            %
            % In addition to the rate function, the method calculates: a. the probability of an event post stimulation:
            % #stimulations after which the event occurred / #total number of stimulations, calculated similarly for
            % control points, b. the event rate post stimulation: number of times the event occurs in the post stimulation
            % epoch / total duration of post stimulation epoch (in minutes). The baseline (“control”) for this variable is
            % the number of times the event occurred in the quiet epoch divided by that epoch’s duration. The post
            % stimulation window for calculating both the probability and the rate is set by timeToCheckAfterSW for slow
            % waves, timeToCheckBeforeSp for spindles and timeToCheckAfterRip for ripples.
            % The input runData is a struct in the length of number of patients (for which the analysis is required).
            % In addition it receives the input parameter fileNameResults which includes the file name into which the
            % results will be saved (optional).
            % Each element (=patient) in runData should include the fields:
            % patientName
            % MTLChans – list of channel indices for which to perform the analysis.
            % ExpDataFileName – name (including path) of the EXP_DATA for the patient.
            % macroMontageFileName - the file name (including path) of the macromontage.
            % DataFolder – The folder in which the raw data files are saved (the method assumes the prefix for the files is
            % CSC, can be changed by the property dataFilePrefix).
            % SpindlesFileNames - name (including path) of the spindle mat files in which the spindle times for the macro
            % channels are saved (the method assumes the name of the file is SpindlesFileNames<#channel index>).
            % SWStaresinaFileName - name (including path) of the slow wave mat files in which the slow waves times for the
            % macro channels are saved (the method assumes the name of the file is SWStaresinaFileName <#channel index>).
            % RipplesFileNames - name (including path) of the ripples mat files in which the ripples times for the macro
            % channels are saved (the method assumes the name of the file is RipplesFileNames<#channel index>).
            %
            % sleepScoringFileName – file name (including path) of the sleep scoring mat file.
            %
            % The output struct results includes all the results of the analysis, which can then be plotted using
            % plotStimEffectCortical. The output struct is a struct with the length of the number of patients (=the length
            % of runData), where each element includes:
            % patientName
            % resultsPerChan – a struct in the length of the number of channels required for the analysis per the patient.
            % Each element in resultsPerChan includes the fields:
            % channelNum
            % nStims – The number of stimulations.
            % baselineSWRate – rate of slow waves in the epoch before the stimulations (#sw/second).
            % baselineSpRate – rate of spindles in the epoch before the stimulations (#spindles/second).
            % baselineRipRate - rate of ripples in the epoch before the stimulations (#sw-spindles/second).
            % baselineRipSWSpRate – rate of ripples-slow waves-spindles triplets in the epoch before the stimulations
            % (#sw-spindles/second).
            % timesAroundStimSW – a cell array with two elements: 1. A matrix where each row is the sw rate function around
            % a stimulation, 2. A matrix where each row is the sw rate function around a control.
            % timesAroundStimSp – a cell array with two elements: 1. A matrix where each row is the spindle rate function
            % around a stimulation, 2. A matrix where each row is the spindle rate function around a control.
            % timesAroundStimRipSWSp - a cell array with two elements: 1. A matrix where each row is the Ripple-SW-spindle
            % triplet rate function around a stimulation, 2. A matrix where each row is the Ripple-SW-spindle triplet rate
            % function around a control.
            % dataAroundStims - a cell array with two elements: 1. A matrix where each row is the raw data around a
            % stimulation, 2. A matrix where each row is the raw data around a control.
            % probsStimRipple – an array with two elements: 1. The probability of ripple post stimulation, 2. The
            % probability of ripple post control.
            % probsStimSW – an array with two elements: 1. The probability of slow wave post stimulation, 2. The
            % probability of slow wave post control.
            % probsStimSp - an array with two elements: 1. The probability of spindle post stimulation, 2. The probability
            % of spindle post control.
            % probsStimRipRSWSp - an array with two elements: 1. The probability of Ripple-SW-spindle post stimulation,
            % 2. The probability of Ripple-SW-spindle post control.
            % eventsPerMinRip – an array with two elements: 1. The event rate of ripple post stimulation, 2. The event rate
            % of ripple post control.
            % eventsPerMinSW – an array with two elements: 1. The event rate of slow waves post stimulation, 2. The event
            % rate of slow waves post control.
            % eventsPerMinSp - an array with two elements: 1. The event rate of spindles post stimulation, 2. The event
            % rate of spindles post control.
            % eventsPerMinRipSWSp - an array with two elements: 1. The event rate of ripple-SW-spindle post stimulation,
            % 2. The event rate of ripple-SW-spindle post control.
            % dataSTD – the STD of the raw data, required for the figures production.
            
            
            if nargin < 3
                fileNameResults = '';
            end
            
            minDistSWSpindle = obj.minDistSWSpindle*obj.samplingRate/1000;
            maxDistSWSpindle = obj.maxDistSWSpindle*obj.samplingRate/1000;
            minDistRippleSW = obj.minDistRippleSW*obj.samplingRate/1000;
            maxDistRippleSW = obj.maxDistRippleSW*obj.samplingRate/1000;
            avgFigBinning = obj.avgFigBinning*obj.samplingRate/1000;
            timeToCheckAfterSW = obj.timeToCheckAfterSW*obj.samplingRate/1000;
            timeToCheckAfterSp = obj.timeToCheckAfterSp*obj.samplingRate/1000;
            timeToCheckAfterRip = obj.timeToCheckAfterRip*obj.samplingRate/1000;
            timeToCheckBeforeSW = obj.timeToCheckBeforeSW*obj.samplingRate/1000;
            timeToCheckBeforeSp = obj.timeToCheckBeforeSp*obj.samplingRate/1000;
            timeToCheckBeforeRip = obj.timeToCheckBeforeRip*obj.samplingRate/1000;
            timeToCheckMax = max([timeToCheckAfterSW,timeToCheckAfterSp,timeToCheckAfterRip]);
            midTimeRangeAfterStim = obj.midTimeRangeAfterStim*obj.samplingRate;
            
            controlDist = obj.controlDistRip*obj.samplingRate/1000;
            
            %go over required patients
            nPatients = length(runData);
            for iPatient = 1:nPatients
                
                disp(['Patient ',runData(iPatient).patientName,' ',num2str(iPatient),'/',num2str(nPatients)]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX;
                
                %load macro montage for area name
                macroMontage = load(runData(iPatient).macroMontageFileName);
                macroMontage = macroMontage.MacroMontage;
                
                %load sleep scoring
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                % dataDuration = stimTimes(end)+timeToCheckMax;
                % extend duration to include post-stim session
                dataDuration = floor(expData.Session_start_end_msec(end,2));
                if sum(stimTimes > dataDuration)
                    warning('%d stims occur after end of session', sum(stimTimes > dataDuration))
                end
                stimTimes(stimTimes > dataDuration) = [];
                
                %analyze only for stim times during NREM
                if ~isempty(sleepScoring)
                    stimTimes = stimTimes(sleepScoring(round(stimTimes))==1);
                end
                
                nStims = length(stimTimes);
                
                resultsPerChan = [];
                nChans = length(runData(iPatient).MTLChans);
                
                
                
                %find indices of quiet epochs - times between stimulations
                %beyond midTimeRangeAfterStim
                quietEpochsBetween = zeros(1,dataDuration);
                stimDiffs = diff(stimTimes);
                stimIndsWithMidPauseAfter = [find(stimDiffs >= midTimeRangeAfterStim) length(stimTimes)];
                stimIndsWithMidPauseAfterTimes = stimTimes(stimIndsWithMidPauseAfter);
                
                stimIndsAfterPauseTimes = stimTimes([1,stimIndsWithMidPauseAfter(1:end-1)+1]); % end of pause block

                
                for iStim = 1:length(stimIndsWithMidPauseAfter)-1
                    quietEpochsBetween(stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim+1:stimTimes(stimIndsWithMidPauseAfter(iStim)+1)-timeToCheckMax-obj.timeBeforeStim) = 1;
                end
                quietEpochsBetweenTimes = find(quietEpochsBetween);
                
                % March 21, add midRange estimate - times that are
                % immediatly after end of stimulation block
                quietEpochsMidRange = zeros(1,dataDuration);
                for iStim = 1:length(stimIndsWithMidPauseAfter)-1
                    quietEpochsMidRange(stimIndsWithMidPauseAfterTimes(iStim)+timeToCheckMax+1:stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim) = 1;
                end
                
                % March 2022 - add an *equal* time for CONTROL of midRange
                % estimate - the last chunck of the quiet epochs.
                %                 quietEpochsMidRange_CONTROL = zeros(1,dataDuration);
                %                 for iStim = 1:length(stimIndsWithMidPauseAfter)-1
                %                     quietEpochsMidRange_CONTROL(stimIndsWithMidPauseAfterTimes(iStim+1)-midTimeRangeAfterStim+1:stimIndsWithMidPauseAfterTimes(iStim+1)-timeToCheckMax) = 1;
                %                 end
                quietEpochsMidRange_CONTROL = zeros(1,dataDuration);
                for iStim = 1:length(stimIndsAfterPauseTimes)-1
                    quietEpochsMidRange_CONTROL(stimIndsAfterPauseTimes(iStim)-midTimeRangeAfterStim+1:stimIndsAfterPauseTimes(iStim)-timeToCheckMax) = 1;
                end
                quietEpochsMidRange_CONTROLTimes = find(quietEpochsMidRange_CONTROL);
                
                
                firstStim = stimTimes(1)-obj.timeBeforeStim;
                durationBeforeStim = sum(sleepScoring(1:firstStim)==obj.sleepEpochs)/obj.samplingRate; %seconds
                
                if ~whatToRun.INCLUDE_BOTH_HEM
                    probe_str = macroMontage(runData(iPatient).probeChan).Area;
                    hemisphere_str = probe_str(1);
                end
                
                %go over required channels
                for iChan = 1:nChans
                    
                    currChan = runData(iPatient).MTLChans(iChan);
                    currArea = macroMontage(currChan).Area;
                    
                    % If the run-mode includes only the PROBE's hemisphere
                    % to focus on local effects
                    if ~whatToRun.INCLUDE_BOTH_HEM
                        if ~strcmp(currArea(1),hemisphere_str)
                            continue;
                        end
                    end
                    
                    resultsPerChan(iChan).channelNum = currChan;
                    resultsPerChan(iChan).nStims = nStims;
                    resultsPerChan(iChan).area = currArea;
                    
                    
                    %load data (to show average data around stim)
                    try
                        data = [runData(iPatient).DataFolder '\',obj.dataFilePrefix, num2str(currChan) '.mat'];
                        data = load(data);
                        data = data.data;
                    catch
                        disp([runData(iPatient).DataFolder '\',obj.dataFilePrefix, num2str(currChan) '.mat doesn''t exist']);
                        continue;
                    end
                    
                    %load spindles
                    if isfield(runData(iPatient), 'SpindlesFileNames') && ~isempty(runData(iPatient).SpindlesFileNames)
                        try
                            spindlesTimes = [runData(iPatient).SpindlesFileNames num2str(currChan) '.mat'];
                            spindlesTimes = load(spindlesTimes);
                            spindlesTimes = spindlesTimes.spindlesTimes;
                        catch
                            disp([runData(iPatient).SpindlesFileNames num2str(currChan) '.mat doesn''t exist']);
                            spindlesTimes = [];
                        end
                    end
                    
                    %load slow waves
                    if isfield(runData(iPatient), 'SWStaresinaFileName') && ~isempty(runData(iPatient).SWStaresinaFileName)
                        try
                            slowWavesTimes = [runData(iPatient).SWStaresinaFileName num2str(currChan) '.mat'];
                            slowWavesTimes = load(slowWavesTimes);
                            slowWavesTimes = slowWavesTimes.slowWavesTimes;
                        catch
                            disp([runData(iPatient).SWStaresinaFileName num2str(currChan) '.mat doesn''t exist']);
                            slowWavesTimes = [];
                        end
                    end
                    
                    %load ripples
                    if isfield(runData(iPatient), 'RipplesFileNames') && ~isempty(runData(iPatient).RipplesFileNames)
                        try
                            ripplesTimes = [runData(iPatient).RipplesFileNames num2str(currChan) '.mat'];
                            ripplesTimes = load(ripplesTimes);
                            ripplesTimes = ripplesTimes.ripplesTimes;
                        catch
                            disp([runData(iPatient).RipplesFileNames num2str(currChan) '.mat doesn''t exist']);
                            ripplesTimes = [];
                        end
                    end
                    
                    spindlesTimes = spindlesTimes(spindlesTimes<(dataDuration+maxDistSWSpindle+maxDistRippleSW));
                    % slowWavesTimes = slowWavesTimes(slowWavesTimes<dataDuration+maxDistRippleSW);
                    slowWavesTimes = slowWavesTimes(slowWavesTimes<dataDuration); % bug fix - jan 2021 - can't take timepoints beyond dataDuration
                    ripplesTimes = ripplesTimes(ripplesTimes<dataDuration);
                    
                    spIndices = zeros(1,dataDuration+maxDistSWSpindle+maxDistRippleSW);
                    spIndices(spindlesTimes) = 1;
                    swSpTimes = [];
                    
                    %first find coupled slow wave - spindles
                    for iSW = 1:length(slowWavesTimes)
                        if sum(spIndices(slowWavesTimes(iSW)+minDistSWSpindle:slowWavesTimes(iSW)+maxDistSWSpindle))>0
                            swSpTimes(end+1) = slowWavesTimes(iSW);
                        end
                    end
                    
                    SWSpIndices = zeros(1,dataDuration+maxDistRippleSW);
                    SWSpIndices(swSpTimes) = 1;
                    RipSWSpTimes = [];
                    %now find trios - ripples-sw-sp
                    for iRipple = 1:length(ripplesTimes)
                        if sum(SWSpIndices(ripplesTimes(iRipple)+minDistRippleSW:ripplesTimes(iRipple)+maxDistRippleSW))>0
                            RipSWSpTimes(end+1) = ripplesTimes(iRipple);
                        end
                    end
                    
                    %find baseline rates
                    baselineRipRate = length(ripplesTimes(ripplesTimes<firstStim))/durationBeforeStim;
                    baselineSWRate = length(slowWavesTimes(slowWavesTimes<firstStim))/durationBeforeStim;
                    baselineSpRate = length(spindlesTimes(spindlesTimes<firstStim))/durationBeforeStim;
                    baselineRipSWSpRate = length(RipSWSpTimes(RipSWSpTimes<firstStim))/durationBeforeStim;
                    resultsPerChan(iChan).baselineRipRate = baselineRipRate;
                    resultsPerChan(iChan).baselineSWRate = baselineSWRate;
                    resultsPerChan(iChan).baselineSpRate = baselineSpRate;
                    resultsPerChan(iChan).baselineRipSWSpRate = baselineRipSWSpRate;
                    
                    %build stimulation triggered rates for all four events
                    %types (ripples, sw, spindles, rip-sw-sp)
                    
                    spRate = movsum(spIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    swIndices = zeros(1,dataDuration);
                    swIndices(slowWavesTimes) = 1;
                    swRate = movsum(swIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    ripIndices = zeros(1,dataDuration);
                    ripIndices(ripplesTimes) = 1;
                    ripRate = movsum(ripIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    RipSWSpIndices = zeros(1,dataDuration);
                    RipSWSpIndices(RipSWSpTimes) = 1;
                    RipSWSpRate = movsum(RipSWSpIndices,avgFigBinning,'Endpoints','fill')/(obj.avgFigBinning/1000);
                    
                    timesAroundStimRip = nan(nStims,timeToCheckBeforeRip+timeToCheckAfterRip+1);
                    probsStimRip = [0 0];
                    eventsNRip = [0 0];
                    
                    timesAroundStimSW = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    probsStimSW = [0 0];
                    eventsNSW = [0 0];
                    
                    timesAroundStimSp = nan(nStims,timeToCheckBeforeSp+timeToCheckAfterSp+1);
                    probsStimSp = [0 0];
                    eventsNSp = [0 0];
                    
                    timesAroundStimRipSWSp = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    probsStimRipSWSp = [0 0];
                    eventsNRipSWSp = [0 0];
                    
                    dataAroundStims = nan(nStims,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    
                    for iStim = 1:nStims
                        
                        %get raw data around stims
                        dataAroundStims(iStim,:) = data(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        
                        timesAroundStimRip(iStim,:) = ripRate(stimTimes(iStim)-timeToCheckBeforeRip:stimTimes(iStim)+timeToCheckAfterRip);
                        currNRipStim = sum(ripIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterRip));
                        if currNRipStim>0
                            probsStimRip(1) = probsStimRip(1)+1;
                            eventsNRip(1) = eventsNRip(1)+currNRipStim;
                        end
                        
                        timesAroundStimSW(iStim,:) = swRate(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        currNSWStim = sum(swIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSW));
                        if currNSWStim>0
                            probsStimSW(1) = probsStimSW(1)+1;
                            eventsNSW(1) = eventsNSW(1)+currNSWStim;
                        end
                        
                        timesAroundStimSp(iStim,:) = spRate(stimTimes(iStim)-timeToCheckBeforeSp:stimTimes(iStim)+timeToCheckAfterSp);
                        currNSpStim = sum(spIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSp));
                        if currNSpStim>0
                            probsStimSp(1) = probsStimSp(1)+1;
                            eventsNSp(1) = eventsNSp(1)+currNSpStim;
                        end
                        
                        timesAroundStimRipSWSp(iStim,:) = RipSWSpRate(stimTimes(iStim)-timeToCheckBeforeSW:stimTimes(iStim)+timeToCheckAfterSW);
                        currNRipSWSpStim = sum(RipSWSpIndices(stimTimes(iStim):stimTimes(iStim)+timeToCheckAfterSW));
                        if currNRipSWSpStim>0
                            probsStimRipSWSp(1) = probsStimRipSWSp(1)+1;
                            eventsNRipSWSp(1) = eventsNRipSWSp(1)+currNRipSWSpStim;
                        end
                        
                    end
                    
                    probsStimRip(1) = probsStimRip(1)/nStims;
                    probsStimSW(1) = probsStimSW(1)/nStims;
                    probsStimSp(1) = probsStimSp(1)/nStims;
                    probsStimRipSWSp(1) = probsStimRipSWSp(1)/nStims;
                    eventsPerMinRip(1) = eventsNRip(1)/(nStims*obj.timeToCheckAfterRip/60000); %events per min
                    eventsPerMinSW(1) = eventsNSW(1)/(nStims*obj.timeToCheckAfterSW/60000); %events per min
                    eventsPerMinSp(1) = eventsNSp(1)/(nStims*obj.timeToCheckAfterSp/60000);
                    eventsPerMinRipSWSp(1) = eventsNRipSWSp(1)/(nStims*obj.timeToCheckAfterSW/60000);
                    
                    %create controls - randon points in the quiet epochs
                    timesAroundStimRipControl = nan(obj.nControls,timeToCheckBeforeRip+timeToCheckAfterRip+1);
                    timesAroundStimSWControl = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    timesAroundStimSpControl = nan(obj.nControls,timeToCheckBeforeSp+timeToCheckAfterSp+1);
                    timesAroundStimRipSWSpControl = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    dataAroundConts = nan(obj.nControls,timeToCheckBeforeSW+timeToCheckAfterSW+1);
                    
                    
                    if whatToRun.controlTimestamps
                        try
                            stim_control_sets = load(runData(iPatient).controlTimestampsFilename,'stim_control_sets');
                            controlTimestamps = stim_control_sets.stim_control_sets.shamTimestamps;
                            controlTimestamps(controlTimestamps > (length(ripRate) - timeToCheckAfterSW)) = [];
                            obj.nControls = length(controlTimestamps);
                        catch
                            warning('control timestamps not found')
                        end
                    end
                    
                    
                    for iControl = 1:obj.nControls
                        
                        if whatToRun.controlTimestamps
                            % controls are matched to stim distance to
                            % slow-waves in the quiet epochs
                            controlInd = controlTimestamps(iControl);
                        else
                            %controls are random indices in the quiet epochs
                            controlInd = quietEpochsBetweenTimes(randi(length(quietEpochsBetweenTimes)));
                        end
                        
                        %get raw data around controls
                        dataAroundConts(iControl,:) = data(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        
                        timesAroundStimRipControl(iControl,:) = ripRate(controlInd-timeToCheckBeforeRip:controlInd+timeToCheckAfterRip);
                        currNRipCont = sum(ripIndices(controlInd:controlInd+timeToCheckAfterRip));
                        if currNRipCont>0
                            probsStimRip(2) = probsStimRip(2)+1;
                            eventsNRip(2) = eventsNRip(2)+currNRipCont;
                            
                        end
                        
                        timesAroundStimSWControl(iControl,:) = swRate(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        currNSWCont = sum(swIndices(controlInd:controlInd+timeToCheckAfterSW));
                        if currNSWCont>0
                            probsStimSW(2) = probsStimSW(2)+1;
                            eventsNSW(2) = eventsNSW(2)+currNSWCont;
                        end
                        
                        timesAroundStimSpControl(iControl,:) = spRate(controlInd-timeToCheckBeforeSp:controlInd+timeToCheckAfterSp);
                        currNSpCont = sum(spIndices(controlInd:controlInd+timeToCheckAfterSp));
                        if currNSpCont>0
                            probsStimSp(2) = probsStimSp(2)+1;
                            eventsNSp(2) = eventsNSp(2)+currNSpCont;
                        end
                        
                        timesAroundStimRipSWSpControl(iControl,:) = RipSWSpRate(controlInd-timeToCheckBeforeSW:controlInd+timeToCheckAfterSW);
                        currNRipSWSpCont = sum(RipSWSpIndices(controlInd:controlInd+timeToCheckAfterSW));
                        %                         nSWSp(iStim) = currNSWSpStim;
                        if currNRipSWSpCont>0
                            probsStimRipSWSp(2) = probsStimRipSWSp(2)+1;
                            eventsNRipSWSp(2) = eventsNRipSWSp(2)+currNRipSWSpCont;
                        end
                    end
                    
                    probsStimRip(2) = probsStimRip(2)/obj.nControls;
                    probsStimSW(2) = probsStimSW(2)/obj.nControls;
                    probsStimSp(2) = probsStimSp(2)/obj.nControls;
                    probsStimRipSWSp(2) = probsStimRipSWSp(2)/obj.nControls;
                    
                    %(1) "control" baseline for the event rate is the event
                    %rate after sham-stimulation points
                    eventsPerMinRip(2) = eventsNRip(2)/(obj.nControls*obj.timeToCheckAfterRip/60000); %events per min
                    eventsPerMinSW(2) = eventsNSW(2)/(obj.nControls*obj.timeToCheckAfterSW/60000); %events per min
                    eventsPerMinSp(2) = eventsNSp(2)/(obj.nControls*obj.timeToCheckAfterSp/60000);
                    eventsPerMinRipSWSp(2) = eventsNRipSWSp(2)/(obj.nControls*obj.timeToCheckAfterSW/60000);
                    
                    %(2) "control" baseline for the event rate is the event rate over entire quiet epoch
                    % FIX - March 2022 - make the control shorter to match
                    % the intermidiate condition below
%                   eventsPerMinRip(3) = sum(ripIndices&quietEpochsBetween)/(sum(quietEpochsBetween)/60000); %events per min
%                   eventsPerMinSW(3) = sum(swIndices&quietEpochsBetween)/(sum(quietEpochsBetween)/60000); %events per min
%                   eventsPerMinSp(3) = sum(spIndices(1:dataDuration)&quietEpochsBetween)/(sum(quietEpochsBetween)/60000);
%                   eventsPerMinRipSWSp(3) = sum(RipSWSpIndices&quietEpochsBetween)/(sum(quietEpochsBetween)/60000);
                    eventsPerMinRip(3) = sum(ripIndices&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000); %events per min
                    eventsPerMinSW(3) = sum(swIndices&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000); %events per min
                    eventsPerMinSp(3) = sum(spIndices(1:dataDuration)&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000);
                    eventsPerMinRipSWSp(3) = sum(RipSWSpIndices&quietEpochsMidRange_CONTROL)/(sum(quietEpochsMidRange_CONTROL)/60000);
                    
                    % event rate (per min) over short period immediatly following stimulation block)
                    eventsPerMinRip(4) =  sum(ripIndices&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000); %events per min
                    eventsPerMinSW(4) = sum(swIndices&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000); %events per min
                    eventsPerMinSp(4) = sum(spIndices(1:dataDuration)&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000);
                    eventsPerMinRipSWSp(4) = sum(RipSWSpIndices&quietEpochsMidRange)/(sum(quietEpochsMidRange)/60000);
                    
                    %get std of data for clipping later
                    resultsPerChan(iChan).dataSTD = nanstd(data(sleepScoring==1));
                    
                    resultsPerChan(iChan).timesAroundStimSW = {timesAroundStimSW, timesAroundStimSWControl};
                    resultsPerChan(iChan).timesAroundStimSp = {timesAroundStimSp, timesAroundStimSpControl};
                    resultsPerChan(iChan).timesAroundStimRip = {timesAroundStimRip, timesAroundStimRipControl};
                    resultsPerChan(iChan).timesAroundStimRipSWSp = {timesAroundStimRipSWSp, timesAroundStimRipSWSpControl};
                    resultsPerChan(iChan).dataAroundStims = {dataAroundStims, dataAroundConts};
                    resultsPerChan(iChan).probsStimSW = probsStimSW;
                    resultsPerChan(iChan).probsStimSp = probsStimSp;
                    resultsPerChan(iChan).probsStimRip = probsStimRip;
                    resultsPerChan(iChan).probsStimRipSWSp = probsStimRipSWSp;
                    resultsPerChan(iChan).eventsPerMinRip = eventsPerMinRip;
                    resultsPerChan(iChan).eventsPerMinSW = eventsPerMinSW;
                    resultsPerChan(iChan).eventsPerMinSp = eventsPerMinSp;
                    resultsPerChan(iChan).eventsPerMinRipSWSp = eventsPerMinRipSWSp;
                    
                end
                
                results(iPatient).resultsPerChan = resultsPerChan;
            end
            
            if ~isempty(fileNameResults)
                save(fileNameResults,'results','-v7.3');
            end
        end
        
        
        
        
        % PLOTTING FOR FIGURES
        
        % Creating a merged dataset to use for Fig 3 plots + plot more
        % conditions against each other for Extended Data Fig 6
         function plotPopulationNeuralPhaseLocking(obj, results, runData, outputFigureFolder, whatTorun)
          
          % For printing -
          % stimSession = obj.phaseLockStimSessionID;
          % baselineInd = obj.phaseLockbaselineInd;
          
          
          nPatients = length(results);
          pt_list = [];
          pVal_mu = []; Nspikes = []; FiringRates = [];
          mu_e = []; r_e = []; r_sq = [];
          area_list = {};
          
          isFrontal = []; isMTL = []; isPHG = []; isTG = []; isP = []; isOc = []; isProbeHem = [];
          AREA_CNT.MTL  = 0; AREA_CNT.Frontal = 0;  AREA_CNT.TG  = 0; AREA_CNT.Par = 0; AREA_CNT.Oc = 0;  AREA_CNT.other  = 0;
          AREA_CNT_SU.MTL  = 0; AREA_CNT_SU.Frontal = 0;  AREA_CNT_SU.TG  = 0; AREA_CNT_SU.Par = 0; AREA_CNT_SU.Oc = 0;  AREA_CNT_SU.other  = 0;
          
          [subGroup, cmap] = defineSubGroups_stimAnatomy();
          
          SU_CNT = zeros(1,nPatients);
          for iPatient = 1:nPatients
                SU_CNT(iPatient) = 0;
                ptNum = str2num(results(iPatient).patientName(2:end));
                depthV = [];
                if isempty(results(iPatient).areasList); continue; end
                
                if ismember(ptNum,subGroup{3}); continue; end % we only have one pt for mixed-phase
                
                try
                    % load macro montage for area name
                    macroMontage = load(runData(iPatient).macroMontageFileName);
                    macroMontage = macroMontage.MacroMontage;
                    macroCh = runData(iPatient).probeChan;
                    macroChArea = macroMontage(macroCh).Area;
                    isLeftProbe = strcmp(macroChArea(1),'L');
                catch
                    disp('probe area-name doesn''t match options')
                end
                IndStruct = results(iPatient).stimCELL.IndStruct;
                
                
                [Narea,Nchannel,Nsu] = size(results(iPatient).spikeCELL);
                M = 4;
                N = ceil(Narea/M+1);
                
                for iArea = 1:Narea
                    
                    for iChannel = 1:Nchannel
                        
                        % subplot(M,N,iArea+1)
                        if isempty(results(iPatient).spikeCELL{iArea,iChannel})
                            continue
                        end
                        
                        for ii_su = 1:Nsu
                            
                            ACTIVE_UNIT = 0;    
                            if isempty(results(iPatient).spikeCELL{iArea,iChannel})
                                continue
                            end
                            if isempty(results(iPatient).spikeCELL{iArea,iChannel,ii_su})
                                continue
                            end
                            
                            % Does this SU have some condition with
                            % decent spiking activity?
                            
                            for ii = 1:length(results(iPatient).spikeCELL{iArea,iChannel,ii_su}.phaseLockStructEpochs)
                                if ~isfield(results(iPatient).spikeCELL{iArea,iChannel,ii_su},'spikeInfo')
                                    break;
                                end
                                phaseLockStructEpochs_local = results(iPatient).spikeCELL{iArea,iChannel,ii_su}.phaseLockStructEpochs{ii};
                                if length(phaseLockStructEpochs_local) >  1
                                    ACTIVE_UNIT = 1;
                                    SU_CNT(iPatient) = SU_CNT(iPatient) + 1;
                                    break
                                end
                            end
                            
                            if ~ACTIVE_UNIT % This means this unit did not pass basic reqs in obj.neuralPhaseLocking()
                                continue
                            end
                            pt_list = [pt_list, ptNum];
                            area = classifyArea(results(iPatient).areasList{iArea});
                            L = length(area_list);
                            area_list{L+1} = results(iPatient).areasList{iArea};
                            if isempty(area); error('need to update area list classification'); end
                            isMTL = [isMTL area.isMTL];
                            isTG = [isTG area.isTG];
                            isP = [isP area.isP];
                            isOc = [isOc area.isOc];
                            isFrontal = [isFrontal area.isFrontal];
                            isProbeHem = [isProbeHem (area.isLeft == isLeftProbe)];
                         
                            unit_index = sum(SU_CNT);
                           
                            
                            % Go over all conditions (see
                            % ase.phaseLockconditions_str)
                            for ii_c = 1:length(results(iPatient).spikeCELL{iArea,iChannel,ii_su}.phaseLockStructEpochs)
   
                                phaseLockStructEpochs_local = results(iPatient).spikeCELL{iArea,iChannel,ii_su}.phaseLockStructEpochs{ii_c};
                                
                                if length(phaseLockStructEpochs_local) >  1
                                    
                                    
                                    Ns = phaseLockStructEpochs_local(IndStruct.phasesWhenSpikesOccur);
                                    Ns = length(Ns{1});
                                    teval_sec = (results(iPatient).spikeCELL{iArea,iChannel,ii_su}.duration_ms(1)/1000);
                                    FRate = Ns/teval_sec; %hz
                                    
                                    % initialize info for this unit/condition
                                    FiringRates(ii_c, unit_index) = FRate;
                                    Nspikes(ii_c,unit_index) = Ns;
                                    mu_e(ii_c,unit_index) = NaN;
                                    r_e(ii_c,unit_index) = NaN;
                                    pVal_mu(ii_c,unit_index) = NaN;
                                    r_sq(ii_c, unit_index) = NaN;

                                    if  FRate > obj.minRatetoIncludeUnitPhaseAnalysis % sufficient  # of spikes and reliable phase locking
                                        
                                        % a = phaseLockStructEpochs_local(IndStruct.mu);
                                        a = phaseLockStructEpochs_local(IndStruct.fitRes);
                                        r_sq(ii_c, unit_index) = a{1}.goodness.rsquare;
                                        ppp = phaseLockStructEpochs_local(IndStruct.pval);
                                        a = phaseLockStructEpochs_local(IndStruct.estPhase);
                                        mu_e(ii_c, unit_index) = a{1};
                                        pVal_mu(ii_c, unit_index) = ppp{1};
                                        % a = phaseLockStructEpochs_local(IndStruct.circ_length);
                                        a = phaseLockStructEpochs_local(IndStruct.r);
                                        r_e(ii_c, unit_index) = a{1}; 
                                    end
                                end
                            end % conditions 
                        end % units
                    end % channel
                end % area list
            end % pt list
            ALL_UNITS.pt_list = pt_list;
            ALL_UNITS.area_list = area_list;
            ALL_UNITS.isMTL = isMTL;
            ALL_UNITS.isTG = isTG;
            ALL_UNITS.isP = isP;
            ALL_UNITS.isOc = isOc;
            ALL_UNITS.isFrontal = isFrontal;
            ALL_UNITS.isProbeHem = isProbeHem;
            
            ALL_UNITS.FiringRates = FiringRates;
            ALL_UNITS.Nspikes = Nspikes;
            ALL_UNITS.mu_e= mu_e;
            ALL_UNITS.r_e= r_e;
            ALL_UNITS.pVal_mu= pVal_mu;
            ALL_UNITS.r_sq = r_sq;

            save(whatTorun.saveFilename,'ALL_UNITS','runData')
            
            % PLOT 
            f0 = newA4figure(sprintf('SU_compare_rsq_pV'));
            axes('position',[0.1 0.1 0.4 0.1])
            for ii = 1; semilogx(ALL_UNITS.pVal_mu(ii,:),ALL_UNITS.r_sq(ii,:),'.'); hold all; end
            line([0.01 0.01], get(gca,'ylim'),'color','k')
            line(get(gca,'xlim'),[0.25 0.25],'color','k')
            xlabel('Rayleigh test, p-value (log)')
            ylabel('R^2 value of fit')
            axis([10^(-80) 1 0 1])
            
            ind = find(ALL_UNITS.r_sq(1,:) > 0.25);
            P = ALL_UNITS.pVal_mu(1,ind);
            pp = sum(P > 0.05)/length(P);
            
            f1 = newA4figure(sprintf('Sup6c_SU_phaseLocking_compare_cond',obj.phaseLockbaselineInd,obj.phaseLockStimSessionID));
            pos1 = [0.1 0.13 0.17 0.2];
            pos2 = [0.1 0.6 0.17 0.2];
            pos3 = [0.4 0.13 0.17 0.2];
            pos4 = [0.4 0.6 0.17 0.2];

            cmap2 = brewermap(10,'set2');
            r_e_col = r_e';
            if (size(r_e_col,1) ~= length(obj.phaseLockconditions_str))
                error('')
            end
            
            clear str DATA
            for ii_a = 1:4
                if ii_a == 1
                    ax = axes('position',pos1);
                    title_str = 'non-mtl units, phase locking sig, cond vs pre-stim';
                    baseline = 1;
                elseif ii_a == 2
                    ax = axes('position',pos2);
                    title_str = 'mtl units, phase locking sig, cond vs pre-stim';
                    baseline = 1;
                elseif ii_a == 3
                    ax = axes('position',pos3);
                    title_str = 'non-mtl units, phase locking sig, cond vs last min of pause';
                    baseline = 6;
                elseif ii_a == 4
                    ax = axes('position',pos4);
                    title_str = 'mtl units, phase locking sig, cond vs last min of pause';
                    baseline = 6;
                end
                clear str DATA A
                COND = [1,3:4,6];
                COND_COLOR = [8,1,4,6,2,8];
                plotCnt = 0;
                for ii = COND
                    
                    if ii == baseline; str{ii} = []; continue; end
                    plotCnt = plotCnt + 1;    
                    A = (r_e_col(:,ii) - r_e_col(:,baseline))./(r_e_col(:,ii) + r_e_col(:,baseline));
                    
                    
                    % exclusion criteria applied in spike-phase
                    % calculation = (1) stability of ISI histogram, (2) N
                    % spikes per period, (3) minimal firing rate per
                    % period.
                    rmv_ind = isnan(r_e_col(:,ii)) | isnan(r_e_col(:,baseline)) | ...
                              ALL_UNITS.pVal_mu(ii,:)' > 0.05 | ALL_UNITS.pVal_mu(baseline,:)' > 0.05 | ...
                              isnan(ALL_UNITS.r_sq(ii,:)') | isnan(ALL_UNITS.r_sq(baseline,:)') | ...
                              r_e_col(:,ii) == 0 | r_e_col(:,baseline) == 0  ;
                          
                     if ismember(ii_a,[1,3]) % remove MTL units
                            rmv_ind = rmv_ind | ALL_UNITS.isMTL';
                     else
                            rmv_ind = rmv_ind | ~ALL_UNITS.isMTL';
                     end
                    
%                     if ismember(ii_a,[1,3]) % top row - remove MTL units
%                         rmv_ind = (isnan(r_e_col(:,ii)) | isnan(r_e_col(:,baseline)) | ...
%                             (r_e_col(:,ii) == 0) | (r_e_col(:,baseline) == 0 ) | ALL_UNITS.isMTL' | ALL_UNITS.r_sq(ii,:)' > 0.07 | ...
%                              ALL_UNITS.pVal_mu(baseline,:)' > 0.05);
%                     else % bottom row - remove ~MTL units
%                         rmv_ind = (isnan(r_e_col(:,ii)) | isnan(r_e_col(:,baseline)) | ...
%                         	 (r_e_col(:,ii) == 0) | (r_e_col(:,baseline) == 0) | ~ALL_UNITS.isMTL' | ALL_UNITS.pVal_mu(ii,:)' > 0.05 | ...
%                              ALL_UNITS.pVal_mu(baseline,:)' > 0.05 );
%                     end
                    A(rmv_ind) = [];
                    DATA{ii} = A;
                    hold all
                    % h = distributionPlot(DATA{ii},'histOri','right','color',cmap2(5+9*ii,:),...
                    %    'showMM',2, 'xValues',ii, 'histOpt',0,'divFactor',[-1.05:0.1:1.05] )
                    violinplot(plotCnt, DATA{ii},'color',cmap2(COND_COLOR(ii),:));
                    str{plotCnt} = sprintf('(%d) N=%d / p=%1.1e',ii,length(A),signrank(A));
                    
                end
                set(gca,'xlim',[0.5 3.5],'xtick',[1:length(COND)],'xticklabels',str,'XTickLabelRotation',35,'ylim',[-.8 .8],'fontsize',12,...
                'ytick',[-0.8,0 0.8])
                plot(get(gca,'xlim'),[0 0],'--k')
                title(title_str,'fontsize',6)
                
            end
            
            XLIM = get(gca,'xlim');
            YLIM = get(gca,'ylim');
            clear str
            str{1} = 'conditions:';
            for ii = COND
                str{ii+1} = sprintf('cond %d: %s',ii,obj.phaseLockconditions_str{ii});
            end
            n = length(str);
            str{end+1} = sprintf('pval th = %1.2f',obj.pThresh);
            str{end+1} = sprintf('min firing rate = %1.1fHz',obj.minRatetoIncludeUnitPhaseAnalysis);
            str{end+1} = sprintf('total of units (before cherry picking) -%d',sum(SU_CNT));
            text(XLIM(2)+diff(XLIM)/5,YLIM(1),str);            
            
            % Plot without excluding random phase (no p-value exclusion)
            % PLOT 
            f1 = newA4figure(sprintf('SU_phaseLocking_compare_cond_NoPvalueExclusion',obj.phaseLockbaselineInd,obj.phaseLockStimSessionID));
            pos1 = [0.1 0.13 0.2 0.2];
            pos2 = [0.1 0.6 0.2 0.2];
            pos3 = [0.4 0.13 0.2 0.2];
            pos4 = [0.4 0.6 0.2 0.2];

            cmap2 = colormap(hsv);
            r_e_col = r_e';
            if (size(r_e_col,1) ~= length(obj.phaseLockconditions_str))
                error('')
            end
            
            clear str DATA
            for ii_a = 1:4
                if ii_a == 1
                    ax = axes('position',pos1);
                    title_str = 'non-mtl units, phase locking sig, cond vs pre-stim';
                    baseline = 1;
                elseif ii_a == 2
                    ax = axes('position',pos2);
                    title_str = 'mtl units, phase locking sig, cond vs pre-stim';
                    baseline = 1;
                    elseif ii_a == 3
                    ax = axes('position',pos3);
                    title_str = 'non-mtl units, phase locking sig, cond vs last min of pause';
                    baseline = 6;
                    elseif ii_a == 4
                    ax = axes('position',pos4);
                    title_str = 'mtl units, phase locking sig, cond vs last min of pause';
                    baseline = 6;
                end
                clear str DATA A
                for ii = 1:length(obj.phaseLockconditions_str)
                    
                    if ii == baseline; str{ii} = []; continue; end
                        
                    A = (r_e_col(:,ii) - r_e_col(:,baseline))./(r_e_col(:,ii) + r_e_col(:,baseline));

                    if ismember(ii_a,[1,3]) % top row - remove MTL units
                        rmv_ind = (isnan(r_e_col(:,ii)) | isnan(r_e_col(:,baseline)) | ...
                            (r_e_col(:,ii) == 0) | (r_e_col(:,baseline) == 0 )| ALL_UNITS.isMTL');
                    else % bottom row - remove ~MTL units
                        rmv_ind = (isnan(r_e_col(:,ii)) | isnan(r_e_col(:,baseline)) | ...
                        	 (r_e_col(:,ii) == 0) | (r_e_col(:,baseline) == 0) | ~ALL_UNITS.isMTL' );
                    end
                    A(rmv_ind) = [];
                    DATA{ii} = A;
                    hold all
                    % h = distributionPlot(DATA{ii},'histOri','right','color',cmap2(5+9*ii,:),...
                    %    'showMM',2, 'xValues',ii, 'histOpt',0,'divFactor',[-1.05:0.1:1.05] )
                    violinplot(ii, DATA{ii},'color',cmap2(5+9*ii,:));
                    str{ii} = sprintf('(%d) N=%d / p=%1.1e',ii,length(A),signrank(A));
                    
                end
                set(gca,'xtick',[1:6],'xticklabels',str,'XTickLabelRotation',35,'ylim',[-1 1.2],'fontsize',12,...
                'ytick',[-1,0 1.2])
                title(title_str,'fontsize',7)
                plot(get(gca,'xlim'),[0 0],'--k')

            end
            
            XLIM = get(gca,'xlim');
            YLIM = get(gca,'ylim');
            clear str
            str{1} = 'conditions:';
            for ii = 1:length(obj.phaseLockconditions_str)
                str{ii+1} = sprintf('cond %d: %s',ii,obj.phaseLockconditions_str{ii});
            end
            n = length(str);
            str{end+1} = sprintf('pval th = none');
            str{end+1} = sprintf('min firing rate = %1.1fHz',obj.minRatetoIncludeUnitPhaseAnalysis);
            str{end+1} = sprintf('total of units (before cherry picking) -%d',sum(SU_CNT));
            text(XLIM(2)+diff(XLIM)/5,YLIM(1),str);  
            
            PrintActiveFigs(outputFigureFolder)
            
      end %func

      
      % Main panels for Fig 3 - also additional conditions are plotted
      % here
      function plotFig3_panels_NeuralPhaseLocking(obj, resultsFile, runData, outputFigureFolder, whatTorun)
          
          mm = matfile(resultsFile);
          ALL_UNITS = mm.ALL_UNITS;
%           isTG = ALL_UNITS.isTG & ~ALL_UNITS.isProbeHem;
%           isP = ALL_UNITS.isP & ~ALL_UNITS.isProbeHem;
%           isFrontal = ALL_UNITS.isFrontal & ~ALL_UNITS.isProbeHem;
%           isOc = ALL_UNITS.isOc & ~ALL_UNITS.isProbeHem;
%           isMTL =  ALL_UNITS.isMTL & ALL_UNITS.isProbeHem;
          isTG = ALL_UNITS.isTG ;
          isP = ALL_UNITS.isP ;
          isFrontal = ALL_UNITS.isFrontal ;
          isOc = ALL_UNITS.isOc ;
          isMTL =  ALL_UNITS.isMTL ;
          
          % PLOT
          for ii_f = 1:3
              
              if ii_f == 1
                  obj.phaseLockbaselineInd = 1;
                  obj.phaseLockStimSessionID = 2 ;
              elseif ii_f == 2
                  obj.phaseLockbaselineInd = 6;
                  obj.phaseLockStimSessionID = 3 ;
              else
                  obj.phaseLockbaselineInd = 1;
                  obj.phaseLockStimSessionID = 3 ;
              end
              f1 = newA4figure(sprintf('Fig3_phaseLocking_revision_cond_%d vs %d',obj.phaseLockbaselineInd,obj.phaseLockStimSessionID));
              
              
              N = 4; M = 5;
              
              %% panel C
              % 1 - 'on', 2 - 'off', '3' - P > 0.05
              % 4 - stim, 5 - baseline
              
              cmap = getUnifiedColors;

              
              pVal_mu_pre = ALL_UNITS.pVal_mu(obj.phaseLockbaselineInd,:);
              mu_e_pre = ALL_UNITS.mu_e(obj.phaseLockbaselineInd,:);
              r_e_pre = ALL_UNITS.r_e(obj.phaseLockbaselineInd,:);
              
              pVal_mu_post = ALL_UNITS.pVal_mu(obj.phaseLockStimSessionID,:);
              mu_e_post = ALL_UNITS.mu_e(obj.phaseLockStimSessionID,:);
              r_e_post = ALL_UNITS.r_e(obj.phaseLockStimSessionID,:);
              
              r_change_vec = (r_e_post - r_e_pre)./(r_e_post + r_e_pre);
              
              invalid_indices_pre = logical(isnan(mu_e_pre) | mu_e_pre == 0);
              invalid_indices_post = logical(isnan(mu_e_post) | mu_e_post == 0);
              
              clear DATA
              for ii_p = 1:4
                  if ii_p == 1
                  axes('position',[0.1,0.85,0.1,0.1])
                  % pie chart - MTL - baseline
                  Nc = sum(logical(isMTL) & ~invalid_indices_pre);
                  cell_indices = logical(isMTL) &  ~invalid_indices_pre & (pVal_mu_pre < obj.pval_th);
                  DATA{1} = wrapTo2Pi(mu_e_pre(cell_indices));
                  
                  pie_data(3) = sum( pVal_mu_pre(logical(isMTL)) > obj.pval_th );
                  
              elseif ii_p == 2
                  % pie chart - MTL - stim effect
                  axes('position',[0.2,0.85,0.1,0.1])
                  Nc = sum(logical(isMTL)& ~invalid_indices_post );
                  cell_indices = logical(isMTL) & ~invalid_indices_post &  (pVal_mu_post < obj.pval_th);
                  DATA{1} = wrapTo2Pi(mu_e_post(cell_indices));
                  
                  pie_data(3) = sum(pVal_mu_post(logical(isMTL)) > obj.pval_th  );
                  
              elseif ii_p == 3
                  % pie chart - ~MTL - baseline
                  axes('position',[0.1,0.75,0.1,0.1])
                  Nc = sum(logical(~isMTL)& ~invalid_indices_pre);
                  cell_indices = logical(~isMTL) & ~invalid_indices_pre & (pVal_mu_pre < obj.pval_th);
                  DATA{1} = wrapTo2Pi(mu_e_pre(cell_indices) );
                  
                  pie_data(3) = sum( pVal_mu_pre(logical(~isMTL)) > obj.pval_th );
                  
              elseif ii_p == 4
                  % pie chart - ~MTL - stim effect
                  
                  axes('position',[0.2,0.75,0.1,0.1])
                  Nc = sum(logical(~isMTL) & ~invalid_indices_post);
                  cell_indices = logical(~isMTL) & ~invalid_indices_post &  logical(pVal_mu_post < obj.pval_th);
                  DATA{1} = wrapTo2Pi(mu_e_post(cell_indices) );
                  
                  pie_data(3) = sum( pVal_mu_post(logical(~isMTL)) > obj.pval_th  );
                  
                  
              end
              % 'ON' phase
              pie_data(1) = sum((obj.minUpPhase_360 <= DATA{1}) & (DATA{1} <= obj.maxUpPhase_360) & ~isnan(DATA{1}));
              % 'OFF' phase
              pie_data(2) = sum((obj.minUpPhase_360 > DATA{1}) | (DATA{1} > obj.maxUpPhase_360) & ~isnan(DATA{1}));
              
              pie_data(3) = Nc - sum(pie_data(1:2));
              pie_data(pie_data == 0) = eps;
              str{ii_p} = sprintf('(%d) pie data - ON - %2.1f, OFF - %2.1f, rand - %2.1f, N = %d',ii_p,100*pie_data/sum(pie_data),sum(pie_data));
              pp = pie(pie_data, [1,1,1],{'','',''});
              pp(1).FaceColor = cmap(1,:);
              pp(1).EdgeColor = [0,0,0];
              pp(3).FaceColor = cmap(2,:);
              pp(3).EdgeColor = [0,0,0];
              
              if length(pp)>5
                  pp(5).EdgeColor = [0,0,0];
                  pp(5).FaceColor = cmap(3,:);
              end
              
              axis off
          end
          XLIM = get(gca,'xlim');
          YLIM = get(gca,'ylim');
          text(XLIM(2)+diff(XLIM)/5,YLIM(2),str);
       
          %% Panel D
          % R-change histogram per brain area

              
            bin_step = 0.1;
            bins1 = [0:-bin_step:-1];
            bins2 = [bin_step:bin_step:1];
            bins = [fliplr(bins1) bins2];
                          
            cmap2 = colormap(hsv);

            % Breaking down to anatomical areas
            % DATAx{1} = r_change_vec(logical(isMTL | isTG));
            % DATAx{2} = r_change_vec(logical(isFrontal | isP | isOc));
            for ii_a = 1:3
                if ii_a == 1
                    % ids = logical(isMTL | isTG);
                    ids = logical(isMTL);
                    
                    pos = [0.7813    0.4    0.1237    0.1];
                    bar_color = cmap2(45,:);
                elseif ii_a == 2
                    ids = logical(~isMTL);
                    
                    pos = [0.7813    0.24    0.1237    0.1];
                    bar_color = 'k';
                else
                    ids = logical(1:length(isMTL));
                    pos = [0.7813    0.08    0.1237    0.1];
                    pos2 = [0.6  0.08    0.1237    0.1];
                    if ii_f == 2
                        bar_color = 'm';
                    else
                        bar_color = 'y'
                    end
                end
                
                ids = logical(ids) & (pVal_mu_pre < obj.pval_th) & (pVal_mu_post < obj.pval_th);
                A = r_e_pre(ids);
                B = r_e_post(ids);
                rmv_ind = isnan(A) | isnan(B) | A == 0 | B == 0;
                A(rmv_ind) = [];
                B(rmv_ind) = [];
                DATAx{ii_a} = (B-A)./(B+A);
                
                axes('position',pos)
                if length(A) > 5
                shuffledDist = shuffleData(A, B,10000);
                [p_shuff(ii_a),h,stats] = ranksum(shuffledDist, DATAx{ii_a});
                end
                clear counts
                [counts,edges] = histcounts(shuffledDist,bins);
                h2 =  histogram('BinEdges',edges,'BinCounts',counts,'facecolor',[0.8 0.8 0.8],'Normalization','probability');
                hold on
                clear counts
                [counts,edges] = histcounts(DATAx{ii_a}',bins);               
                h1 = histogram('BinEdges',edges,'BinCounts',counts,'facecolor',bar_color,'Normalization','probability');
                h1.FaceAlpha = 0.5;
                
                box off
                if ii_a == 1
                    set(gca,'ylim',[0 .5],'ytick',[ 0 .5])
                else
                    set(gca,'ylim',[0 .5],'ytick',[ 0 .5])
                end
                set(gca,'XLim',[-1 1],'xtick' ,[-1:1:1]);
                box off
                xrange = get(gca,'XLim'); yrange = get(gca,'YLim');
                lineh = line([0 0], [0 yrange(2)],'color','k','linewidth',2);
                lineh = line(median(DATAx{ii_a})*ones(1,2), [0 yrange(2)],'color','r','linewidth',2);
                
                
                if ii_a == 3 & ii_f == 2
                    % plot another shuffled istribution 
                    NeuralPhaseLocking_BS_OutputFile = 'E:\Data_p\ClosedLoopDataset\SWTriggeredSpikeHistResults\SWTriggeredSpikeHist_allSes_AllPts_stableUnits_BS_1000.mat';
                    mm = matfile(NeuralPhaseLocking_BS_OutputFile);
                    bs_results = mm.bs_results;
                    gt_diff = bs_results.gt_diff;
                    diff_bs = bs_results.diff_bs;
                    pval_bs = bs_results.pval_bs;
                    axes('position',pos2)
                    ind =  find(sum(pval_bs(1:2,:,1) < 0.05) == 2);
                    B = gt_diff(ind);
                    A = diff_bs(ind,:);
                    [counts,edges] = histcounts(A(:),bins);
                    h2 =  histogram('BinEdges',edges,'BinCounts',counts,'facecolor',[0.8 0.8 0.8],'Normalization','probability');
                    hold on
                    clear counts
                    [counts,edges] = histcounts(DATAx{ii_a}',bins);
                    h1 = histogram('BinEdges',edges,'BinCounts',counts,'facecolor',bar_color,'Normalization','probability');
                    h1.FaceAlpha = 0.5;
                    [h,p_bs] = kstest2(DATAx{ii_a},A(:));
                    disp(sprintf('kstest2 - bootstrap vs orig (all cells) - p = %2.2e',p_bs))
                    
                    box off
                    if ii_a == 1
                        set(gca,'ylim',[0 .5],'ytick',[ 0 .5])
                    else
                        set(gca,'ylim',[0 .5],'ytick',[ 0 .5])
                    end
                    set(gca,'XLim',[-1 1],'xtick' ,[-1:1:1]);
                    box off
                    xrange = get(gca,'XLim'); yrange = get(gca,'YLim');
                    lineh = line([0 0], [0 yrange(2)],'color','k','linewidth',2);
                    
                
                end

            end
            

            pos = [0.7813    0.65    0.1237    0.1];
            axes('position',pos)
            bins_plot = edges + diff(edges(1:2))/2;
            bins_plot(end) = [];
            clear counts
            [counts(:,1),edges] = histcounts(DATAx{1}',bins);
            [counts(:,2),edges] = histcounts(DATAx{2}',bins);
            barh = bar(bins_plot,counts,'stacked');
            barh(1).FaceColor = cmap2(46,:);
            barh(2).FaceColor = [0.3 0.3 0.3];
            
            box off
            hold on
            if ii_f == 3
                set(gca,'ylim',[0 25],'ytick',[ 0 25])
            end
            xrange = get(gca,'XLim'); yrange = get(gca,'YLim');
            lineh = line([0 0], [0 yrange(2)],'color','k','linewidth',2);
            
            % Signtest used in main test
            pval(1) = signrank(DATAx{1});
            pval(2) = signrank(DATAx{2});
            [pval(3), h] = signrank([DATAx{1}, DATAx{2}]);
            [pval(4), h] = ranksum(DATAx{1}, DATAx{2});
            N3 = sum(~isnan([DATAx{1}, DATAx{2}]));
            disp(sprintf('stats for histogram, p = %2.2e N = %d (MTL)',pval(1),sum(~isnan(DATAx{1}))))
            disp(sprintf('stats for histogram, p = %2.2e N = %d (OF,Par,OC, TG)',pval(2),sum(~isnan(DATAx{2}))))
            disp(sprintf('stats for all units, p = %2.2e rnaksum, N = %d',pval(3),N3))
            disp(sprintf('shuffle test - mtl P = %2.2e, ~MTL P = %2.2e',p_shuff(1),p_shuff(2)))
            
            
            
            XLIM = get(gca,'xlim');
            YLIM = get(gca,'ylim');
            t1 = obj.phaseLockconditions_str(obj.phaseLockbaselineInd);
            t2 = obj.phaseLockconditions_str(obj.phaseLockStimSessionID);
            clear str
            str{1} = sprintf('baseline = %s, stim = %s',t1{1},t2{1});
            str{2} = sprintf('MTL - p = %1.1e (%d)',pval(1),length(DATAx{1}));
            str{3} = sprintf('~MTL - p = %1.1e (%d)',pval(2),length(DATAx{2}));
            str{4} = sprintf('all units - p = %1.1e',pval(3));
            str{5} = sprintf('shuffle test - mtl P = %2.2e, ~MTL P = %2.2e',p_shuff(1),p_shuff(2));
            str{6} = sprintf('shuffle test - all P = %2.2e',p_shuff(3));
            text(XLIM(1)-diff(XLIM)/5,YLIM(1)-diff(YLIM)/1.7,str,'fontsize',7.5);

            %% Sup figure panels

            % Breaking down to anatomical areas
            % R-change histogram per brain area
            axes('position',[0.43,0.25,0.2866,0.1577])
            
            cmap2 = colormap(hsv);
          
            for ii = 1:5
                if ii == 1
                    ids = logical(isMTL);
                elseif ii == 2
                    ids = logical(isTG);
                elseif ii == 3
                    ids = logical(isFrontal);
                elseif ii == 4
                    ids = logical(isP);
                elseif ii == 5
                    ids = logical(isOc);
                end
                A = r_e_pre(logical(ids) & (pVal_mu_post <= obj.pval_th));
                B = r_e_post(logical(ids) & (pVal_mu_post <= obj.pval_th));
                rmv_ind = isnan(A) | isnan(B) | A == 0 | B == 0;
                A(rmv_ind) = [];
                B(rmv_ind) = [];
                DATA{ii} = (B-A)./(B+A);
                if isempty(DATA{ii}); continue; end
                pval_area(ii) = signrank(DATA{ii});

                if ~isempty(DATA{ii})
                    distributionPlot(DATA{ii}','histOri','right','color',cmap2(ii*9,:),'showMM',0, 'xValues',ii, 'histOpt',0,'divFactor',[-1.05:0.1:1.05] )
                end
            end
            area_str_cell = {'MTL','Am,HG','OF,AF,In','Par','Occ'};
            for ii = 1:5
                if isempty(DATA{ii}); continue; end
                a = area_str_cell(ii);
                str_text{ii} = sprintf('signrank %s p = %2.2f',a{1},pval_area(ii));
            end
            axis tight
            set(gca,'xtick',[1:5],'xticklabels',area_str_cell)
            xtickangle( 45 )
            hold on
            line(get(gca,'xlim'),[0,0],'color','k')
            line([1.45,1.45],get(gca,'ylim'),'color','k')
            for ii = 1:length(DATA)
                line([ii+0.55,ii+0.55],get(gca,'ylim'),'color','k')
            end
            set(gca,'ylim',[-1,1.05],'ytick',[-1:0.5:1])
            ylabel('Phase locking change')
            XLIM = get(gca,'xlim');
            YLIM = get(gca,'ylim');
            text(XLIM(1),YLIM(1)-0.75*diff(YLIM),str_text)
            
            
            %%         
            axes('position',[0.1    0.25    0.2866    0.1577])
            cla
            for ii = 1:5
                if ii == 1
                    cell_indices = logical(isMTL);
                elseif ii == 2
                    cell_indices = logical(isTG);
                elseif ii == 3
                    cell_indices = logical(isFrontal);
                elseif ii == 4
                    cell_indices = logical(isP);
                elseif ii == 5
                    cell_indices = logical(isOc);
                end
                DATA{1} = wrapTo2Pi(mu_e_pre(cell_indices &  (pVal_mu_pre <= obj.pval_th) ));
                rmv_ind = isnan(DATA{1}) | DATA{1} == 0;
                DATA{1}(rmv_ind) = [];
                
                DATA{2} = wrapTo2Pi(mu_e_post(cell_indices &  (pVal_mu_post <= obj.pval_th) ));
                rmv_ind = isnan(DATA{2}) | DATA{2} == 0;
                DATA{2}(rmv_ind) = [];
                
                
                bin_plots = [0:0.3:2*pi];
                if ~isempty(DATA{2})
                     distributionPlot(DATA{2}','histOri','right','color',cmap(4,:),'widthDiv',[2 2],'showMM',0,'histOpt',0,'divFactor',bin_plots,'globalNorm',0 , 'xValues',ii)
                end
                if ~isempty(DATA{1})
                    distributionPlot(DATA{1}','histOri','left','color',cmap(5,:),'widthDiv',[2 1],'showMM',0,'histOpt',0,'divFactor',bin_plots ,'globalNorm',0, 'xValues',ii)
                end
                
            end
                
            set(gca,'ylim',[0,2*pi],'ytick',[0,pi,2*pi],'yticklabel',{'0','180','360'})
            set(gca,'xtick',[1:5],'xticklabels',{'MTL','OF,AF,In','Am,HG','Par','Occ'})
            xtickangle( 45 )
            ylabel('Mean phase (deg)')
            
            for ii = 1:length(DATA)
                line([ii ii],get(gca,'ylim'),'color','k')
            end
            
            hold on;
            t = [-pi:0.0001:3*pi];
            example = 0.3*cos(t)+0.3;
            LW_cos = 4;
            ind = t >= obj.minUpPhase_360 & t <= obj.maxUpPhase_360;
            plot(example(ind), t(ind),'linewidth',LW_cos,'color',cmap(1,:));
            hold on
            ind = t < obj.minUpPhase_360;
            plot(example(ind), t(ind),'linewidth',LW_cos,'color',cmap(2,:));
            ind = t > obj.maxUpPhase_360;
            plot(example(ind), t(ind),'linewidth',LW_cos,'color',cmap(2,:));
            
            line(get(gca,'xlim'),[pi pi],'color','k')
            
          end % plots
          
          PrintActiveFigs(outputFigureFolder)
          
          
      end % func
        
    end