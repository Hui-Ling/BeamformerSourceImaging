function varargout = process_beamformer_con_speedup( varargin )
% PROCESS_BEAMFORMER_TEST: 

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: 
eval(macro_method)
%macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Beamformer-based connectivity imaging';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 3005;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','matrix'};
    sProcess.OutputTypes = {'results','results'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    
   %sProcess.options.sep.Type    = 'label';
%     % Definition of the options
%     sProcess.options.oriconstraint.Comment = {'Unconstrained', 'Cortical Constrained'};
%     sProcess.options.oriconstraint.Type    = 'radio';
%     sProcess.options.oriconstraint.Value   = 1;

    sProcess.options.result_comm.Comment    = 'Comment: ';
    sProcess.options.result_comm.Type       = 'text';
    sProcess.options.result_comm.Value      = '';
    % === MEASURE ===
    sProcess.options.measure.Comment    = {'SILSC (Corr)', 'BIPAC (ESC)', 'Methods (Measures):'};
    sProcess.options.measure.Type       = 'radio_line';
    sProcess.options.measure.Value      = 1;
    % Options: Time-freq
    sProcess.options.labelC.Comment = '<HTML><B><U>Options for BIPAC (ESC)</U></B>:';
    sProcess.options.labelC.Type     = 'label';
    % === NESTING FREQ
    sProcess.options.nesting.Comment = 'Nesting frequency band (low, data, Files A):';
    sProcess.options.nesting.Type    = 'range';
    sProcess.options.nesting.Value   = {[2, 30], 'Hz', 2};
    % === NESTED FREQ
    sProcess.options.nested.Comment = 'Nested frequency band (high, ref, Files B):';
    sProcess.options.nested.Type    = 'range';
    sProcess.options.nested.Value   = {[40, 150], 'Hz', 2};
    % Options: number of randomization test
    sProcess.options.nRand.Comment = 'Number of randomization (very time consuming)';
    sProcess.options.nRand.Type   = 'value';
    sProcess.options.nRand.Value   = {0, '', 0};    
    sProcess.options.nRand.Hidden   = 1;  
%     sProcess.options.ref_replicate_num.Comment = 'Number to replicate reference signals:';
%     sProcess.options.ref_replicate_num.Type    = 'value';
%     sProcess.options.ref_replicate_num.Value   = {1,'times',0};
%     
%     sProcess.options.ref_replicate_interval.Comment = 'Interval of replicated reference signals:';
%     sProcess.options.ref_replicate_interval.Type    = 'value';
%     sProcess.options.ref_replicate_interval.Value   = {0.01,'ms',1};

    % Separator
    sProcess.options.sep2.Type     = 'separator';
    sProcess.options.sep2.Comment = '  ';
    
    sProcess.options.labelB.Comment = '<HTML><B><U>Estimation options (Files A)</U></B>:';
    sProcess.options.labelB.Type = 'label';
%     sProcess.options.beamformertype.Comment = {'Vector-type beamformer', 'Scalar-type beamformer'};
%     sProcess.options.beamformertype.Type    = 'radio';
%     sProcess.options.beamformertype.Value   = 1;
    % === MINIMUM VARIANCE
       sProcess.options.minvar_time.Comment = 'Minimum variance time window: ';
       sProcess.options.minvar_time.Type    = 'timewindow';
       sProcess.options.minvar_time.Value   = [];    

    % === ACTIVE TIME RANGE
    sProcess.options.corr_range.Comment = 'Time range of interest: ';
    sProcess.options.corr_range.Type    = 'timewindow';
    sProcess.options.corr_range.Value   = [];
    % === ACTIVE TIME WINDOW SIZE
    sProcess.options.active_window_size.Comment = 'Sliding window size: ';
    sProcess.options.active_window_size.Type    = 'value';
    sProcess.options.active_window_size.Value   = {0.1, 'ms', 1};
    % === ACTIVE TEMPORAL RESOLUTION
    sProcess.options.corr_tresolution.Comment = 'Temporal resolution: ';
    sProcess.options.corr_tresolution.Type    = 'value';
    sProcess.options.corr_tresolution.Value   = {0.01, 'ms', 1};
    % === REGULARIZATION
    sProcess.options.reg.Comment = 'Regularization parameter: ';
    sProcess.options.reg.Type    = 'value';
    sProcess.options.reg.Value   = {0.1, '%', 4};
    % === Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    
    % Separator
    sProcess.options.sep.Type     = 'separator';
    sProcess.options.sep.Comment = '  ';   
    
    sProcess.options.labelA.Comment = '<HTML><B><U>Reference signal options (Files B)</U></B>:';
    sProcess.options.labelA.Type     = 'label';
    % === REF TIME TYPE ===
    sProcess.options.ref_range_type.Comment    = {'Fixed time range (Absolute)', 'Range of time delay (Relative)', 'Time range type:'};
    sProcess.options.ref_range_type.Type       = 'radio_line';
    sProcess.options.ref_range_type.Value      = 1;    
    
    % === TIME WINDOW ===
    sProcess.options.ref_range.Comment = 'Time range to extract ref signal: ';
    sProcess.options.ref_range.Type    = 'timewindow';
    sProcess.options.ref_range.Value   = [];
    
    sProcess.options.ref_replicate_tresolution.Comment = 'Temporal resolution to extract ref signal: ';
    sProcess.options.ref_replicate_tresolution.Type    = 'value';
    sProcess.options.ref_replicate_tresolution.Value   = {0.01, 'ms', 1}; 
    
    % === BASELINE TIME RANGE FOR Z STATISTICS CALCULATION
    sProcess.options.baseline_time.Comment = 'Baseline: ';
    sProcess.options.baseline_time.Type    = 'baseline';
    sProcess.options.baseline_time.Value   = [];
    
    % Options: Hanning window
    sProcess.options.isHann.Comment = 'Apply hanning window';
    sProcess.options.isHann.Type    = 'checkbox';
    sProcess.options.isHann.Value   = 0;
    

    %sProcess.options.sep2.Type    = 'label';
%     % Options: Mirror
%     sProcess.options.mirror.Comment = 'Estimating in frequency domain';
%     sProcess.options.mirror.Type    = 'checkbox';
%     sProcess.options.mirror.Value   = 0;
%     sProcess.options.mirror.InputTypes = {'data'};


    
%     % === Low bound
%     sProcess.options.highpassA.Comment = 'Lower cutoff frequency (Data):';
%     sProcess.options.highpassA.Type    = 'value';
%     sProcess.options.highpassA.Value   = {15,'Hz ',2};
%     sProcess.options.highpassA.InputTypes = {'data'};
%     % === High bound
%     sProcess.options.lowpassA.Comment = 'Upper cutoff frequency (Data):';
%     sProcess.options.lowpassA.Type    = 'value';
%     sProcess.options.lowpassA.Value   = {29,'Hz ',2};
%     sProcess.options.lowpassA.InputTypes = {'data'};
%      % === Low bound
%     sProcess.options.highpassB.Comment = 'Lower cutoff frequency (Ref):';
%     sProcess.options.highpassB.Type    = 'value';
%     sProcess.options.highpassB.Value   = {0,'Hz ',2};
%     sProcess.options.highpassB.InputTypes = {'data'};
%     % === High bound
%     sProcess.options.lowpassB.Comment = 'Upper cutoff frequency (Ref):';
%     sProcess.options.lowpassB.Type    = 'value';
%     sProcess.options.lowpassB.Value   = {0,'Hz ',2};
%     sProcess.options.lowpassB.InputTypes = {'data'};  
%     % === TF METHOD  ===
%     sProcess.options.tfmethod.Comment    = {'Hilbert', 'Wavelet', 'STFT', 'TF method:'};
%     sProcess.options.tfmethod.Type       = 'radio_line';
%     sProcess.options.tfmethod.Value      = 1;
%     % === Width 
%     sProcess.options.width.Comment = 'Cycles of wavelet:';
%     sProcess.options.width.Type    = 'value';
%     sProcess.options.width.Value   = {7, 'Cycles', 0}; 
%     % === WINDOW LENGTH
%     sProcess.options.winlength.Comment = 'Estimator window length: ';
%     sProcess.options.winlength.Type    = 'value';
%     sProcess.options.winlength.Value   = {0.128, 'ms ', 1};
%     sProcess.options.winlength.InputTypes = {'data'};
%     % === Low bound
%     sProcess.options.winoverlap.Comment = 'Overlap percentage: ';
%     sProcess.options.winoverlap.Type    = 'value';
%     sProcess.options.winoverlap.Value   = {0.75, '% ', 1};
%     sProcess.options.winoverlap.InputTypes = {'data'};
%     % === Freq band
%     sProcess.options.freqband.Comment = 'Frequency band of interested:';
%     sProcess.options.freqband.Type    = 'text';
%     sProcess.options.freqband.Value   = 'beta';
%     sProcess.options.freqband.InputTypes = {'timefreq'};
    
   
%     % === SCOUTS ===
%     sProcess.options.scouts.Comment = 'Use scouts';
%     sProcess.options.scouts.Type = 'scout';
%     sProcess.options.scouts.Value      = [];
    
%     % === Freq band
%     sProcess.options.DICStest.Comment = 'DICS test';
%     sProcess.options.DICStest.Type    = 'text';
%     sProcess.options.DICStest.Value   = '1';
    
    % === CONNECT INPUT
    %sProcess = process_corr1n('DefineConnectOptions', sProcess, 1);
end

%% ===== DEFINE SCOUT OPTIONS =====
function sProcess = DefineConnectOptions(sProcess, isConnNN) %#ok<DEFNU>
    % === TIME WINDOW ===
    sProcess.options.label1.Comment = '<HTML><B><U>Input options</U></B>:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % === FROM: CONNECTIVITY [1xN] ===
    if ~isConnNN
        % === FROM: REFERENCE CHANNELS ===
        sProcess.options.src_channel.Comment    = 'Source channel: ';
        sProcess.options.src_channel.Type       = 'channelname';
        sProcess.options.src_channel.Value      = 'name';
        sProcess.options.src_channel.InputTypes = {'data'};
        % === FROM: ROW NAME ===
        sProcess.options.src_rowname.Comment    = 'Source rows (names or indices): ';
        sProcess.options.src_rowname.Type       = 'text';
        sProcess.options.src_rowname.Value      = '';
        sProcess.options.src_rowname.InputTypes = {'timefreq', 'matrix'};
    end
    % === TO: SENSOR SELECTION ===
    sProcess.options.dest_sensors.Comment    = 'Sensor types or names (empty=all): ';
    sProcess.options.dest_sensors.Type       = 'text';
    sProcess.options.dest_sensors.Value      = 'MEG, EEG';
    sProcess.options.dest_sensors.InputTypes = {'data'};
    % === SCOUTS ===
    sProcess.options.scouts.Comment = 'Use scouts';
    if isConnNN
        sProcess.options.scouts.Type = 'scout_confirm';
    else
        sProcess.options.scouts.Type = 'scout';
    end
    sProcess.options.scouts.Value      = [];
    sProcess.options.scouts.InputTypes = {'results'};
    % Atlas: surface/volume
    sProcess.options.isvolume.Comment = '';
    sProcess.options.isvolume.Type    = 'checkbox';
    sProcess.options.isvolume.Value   = 0;
    sProcess.options.isvolume.Hidden  = 1;
    % === SCOUT FUNCTION ===
    sProcess.options.scoutfunc.Comment    = {'Mean', 'Max', 'PCA', 'Std', 'All', 'Scout function:'};
    sProcess.options.scoutfunc.Type       = 'radio_line';
    sProcess.options.scoutfunc.Value      = 2;
    sProcess.options.scoutfunc.InputTypes = {'results'};
    % === SCOUT TIME ===
    sProcess.options.scouttime.Comment    = {'Before', 'After', 'When to apply the scout function:'};
    sProcess.options.scouttime.Type       = 'radio_line';
    sProcess.options.scouttime.Value      = 2;
    sProcess.options.scouttime.InputTypes = {'results'};
end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>
    % Initialize returned list of files
    OutputFiles = {};

    if isempty(sInputsA) || isempty(sInputsB)
        bst_report('Error', sProcess, sInputsA, 'No inputs');
        return;
    end
        
    InputsData = sInputsA ;
    InputsRef = sInputsB ;

    % Get option values
    %isFreq = sProcess.options.mirror.Value;
    OPTIONS.isHann = sProcess.options.isHann.Value;
    OPTIONS.BeamformerType = 2;%sProcess.options.beamformertype.Value;
    %ActiveTime  = sProcess.options.active_time.Value{1};
    % Get option values
    OPTIONS.BaselineTime  = sProcess.options.baseline_time.Value{1};
    OPTIONS.MinVarTime  = sProcess.options.minvar_time.Value{1};
    %OPTIONS.RefTime   = sProcess.options.ref_range.Value{1};
    OPTIONS.CORRrange   = sProcess.options.corr_range.Value{1};
    OPTIONS.Reg         = sProcess.options.reg.Value{1};
    OPTIONS.SensorTypes = sProcess.options.sensortypes.Value; 
    OPTIONS.CORRTResolu = sProcess.options.corr_tresolution.Value{1};
    OPTIONS.WinSize = sProcess.options.active_window_size.Value{1};
%     OPTIONS.BandBoundsData = sProcess.options.nesting.Value;
%     OPTIONS.BandBoundsRef =  sProcess.options.nested.Value;
%     OPTIONS.WinOverlap = sProcess.options.winoverlap.Value{1};
%     OPTIONS.SegmentLength = sProcess.options.winlength.Value{1};
%     OPTIONS.Width = sProcess.options.width.Value{1};
    % Get and check frequencies
    OPTIONS.BandBoundsData = sProcess.options.nesting.Value{1};
    OPTIONS.BandBoundsRef  = sProcess.options.nested.Value{1};
    OPTIONS.RefRangeType = sProcess.options.ref_range_type.Value;
    OPTIONS.Width = 256;
    %OPTIONS.Target = sProcess.options.scouts.Value;
    OPTIONS.nRand = sProcess.options.nRand.Value{1};
 
    OPTIONS.TFmethod = 'hilbert';
%     switch (sProcess.options.tfmethod.Value)
%         case 1, OPTIONS.TFmethod = 'hilbert';
%         case 2, OPTIONS.TFmethod = 'wavelet';
%         case 3, OPTIONS.TFmethod = 'stft';
%     end
    switch (sProcess.options.measure.Value)
        case 1, OPTIONS.measure = 'cor'; OPTIONS.method = 'SILSC';
        case 2, OPTIONS.measure = 'esc'; OPTIONS.method = 'BIPAC';
    end
    
    if strcmp(OPTIONS.measure,'esc') || strcmp(OPTIONS.measure,'cfc') 
        if (min(OPTIONS.BandBoundsData) < 0)
            bst_report('Error', sProcess, [], 'This function cannot be used to estimate PAC for nesting frequencies below 0 Hz.');
            return;
        end
        if (max(OPTIONS.BandBoundsData) > min(OPTIONS.BandBoundsRef)) && min(OPTIONS.BandBoundsRef)~=0 && max(OPTIONS.BandBoundsRef)~=0
            bst_report('Error', sProcess, [], 'The low and high frequency band cannot overlap.');
            return;
        end  
    end
    
    OPTIONS.RepRefNum = 0;%sProcess.options.ref_replicate_num.Value{1};
    OPTIONS.RepRefInterval = 0;%sProcess.options.ref_replicate_interval.Value{1};
    OPTIONS.RefDelayInterval = sProcess.options.ref_replicate_tresolution.Value{1};
    OPTIONS.RefDelayRange = sProcess.options.ref_range.Value{1};
    
    result_comment = sProcess.options.result_comm.Value;
    if ~isempty(result_comment)
        result_comment = [result_comment ': '];
    end
%     MinVarTime = OPTIONS.CORRrange;
    % ===== LOAD THE REFRENCE DATA =====
    % Read the first file in the list, to initialize the loop
    RefMat = in_bst(InputsRef(1).FileName, [], 0);
    %nRefTrials = size(RefMat.Value,1);
    RefFs = 1/(RefMat.Time(2)-RefMat.Time(1));
%     if RefMat.Time(1) > OPTIONS.RefTime(1)
%         bst_report('Warning', sProcess, InputsData, 'The start for time range of reference signal is reset to the first time point of data');
%         OPTIONS.RefTime(1) = RefMat.Time(1);
%     end
%     if RefMat.Time(end) < OPTIONS.RefTime(2)
%         bst_report('Warning', sProcess, InputsData, 'The end for time range of reference signal is reset to the last time point of data');
%         OPTIONS.RefTime(2) = RefMat.Time(end);
%     end
%     RefTimePoints= panel_time('GetTimeIndices', RefMat.Time, [OPTIONS.RefTime(1) OPTIONS.RefTime(2)]);
%     RefLength = length(RefTimePoints);
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
    DataMat = in_bst(InputsData(1).FileName, [], 0);
    if strcmpi(InputsData(1).FileType,'data')
        nChannels = size(DataMat.F,1);   
    else
        nChannels = size(DataMat.TF,1);
    end
    Time = DataMat.Time;
    nTime  = length(Time);
    DataFs = 1/(DataMat.Time(2)-DataMat.Time(1));
    
    CorrWindowSize     = OPTIONS.WinSize;%.RefTime(end)-OPTIONS.RefTime(1);%RefLength;
    if CorrWindowSize <= 0 || OPTIONS.CORRTResolu <= 0
        CorrWindowSize = OPTIONS.CORRrange(2)-OPTIONS.CORRrange(1);
        OPTIONS.CORRTResolu = 0;
    end
    
    if OPTIONS.RepRefNum > 0 && OPTIONS.RepRefInterval > 0
        CorrWindowSize = CorrWindowSize + (OPTIONS.RepRefNum-1)*OPTIONS.RepRefInterval;
    end
    HalfCorrWindowSize = CorrWindowSize/2;
    
    % ===== PROCESS THE TIME WINDOWS =====
    if round(RefFs) ~= round(DataFs)
        bst_report('Error', sProcess, InputsData, 'The sampling rates of FileA and FileB are different.');
    else
        Fs = RefFs;
    end
    if OPTIONS.CORRrange(1) > OPTIONS.CORRrange(2)
        bst_report('Error', sProcess, InputsData, 'The setting of time range of interest is incorrect.');
    end
    
    if OPTIONS.CORRrange(1) < Time(1)
        bst_report('Warning', sProcess, InputsData, 'The start for time range of interest is reset to the first time point of data');
        OPTIONS.CORRrange(1) = Time(1);
    end
    
    if OPTIONS.CORRrange(2) > Time(end)
        bst_report('Warning', sProcess, InputsData, 'The end for time range of interest is reset to the end point of data');
        OPTIONS.CORRrange(2) = Time(end);
    end
    
    if (OPTIONS.CORRrange(1)+CorrWindowSize) > OPTIONS.CORRrange(2) || OPTIONS.CORRTResolu <= 0;
        bst_report('Warning', sProcess, InputsData, 'The active window size is reset to the same as the time range of interest.');
        CorrWindowSize = OPTIONS.CORRrange(2) - OPTIONS.CORRrange(1);
        OPTIONS.CORRTResolu = 0;
    end     
    
    CorrRangePoints = panel_time('GetTimeIndices', Time, [OPTIONS.CORRrange(1)+CorrWindowSize OPTIONS.CORRrange(2)]);
    if length(CorrRangePoints)<=1
        nCORR = 1;
    else 
        nCORR = length((OPTIONS.CORRrange(1)+CorrWindowSize):OPTIONS.CORRTResolu:OPTIONS.CORRrange(2));
        if nCORR == 0  
            bst_report('Error', sProcess, InputsData, 'No correlation windows.');
        end
    end
       
    CorrTimeList= zeros(nCORR,2);
    for i=1:nCORR
        CorrTimeList(i,:) = OPTIONS.CORRrange(1) + OPTIONS.CORRTResolu*(i-1) + [0 CorrWindowSize] ;
    end 
    MinVarRangePoint = panel_time('GetTimeIndices', Time, [OPTIONS.MinVarTime(1) OPTIONS.MinVarTime(2)]);


    % ===== LOAD CHANNEL FILE =====
    % Load channel file
    ChannelMat = in_bst_channel(InputsData(1).ChannelFile);
    % Find the MEG channels
    iChannels = channel_find(ChannelMat.Channel, OPTIONS.SensorTypes);
    
    % ===== LOAD HEAD MODEL =====
    % Get channel study
    [sChannelStudy, iChannelStudy] = bst_get('ChannelFile', InputsData(1).ChannelFile);
    % Load the default head model
    HeadModelFile = sChannelStudy.HeadModel(sChannelStudy.iHeadModel).FileName;
    sHeadModel = load(file_fullpath(HeadModelFile));
    nCorrWindowPoints = length(panel_time('GetTimeIndices', Time, CorrTimeList(1,:)));
    for i=2:nCORR
        nCorrWindowPoints = min(nCorrWindowPoints,length(panel_time('GetTimeIndices', Time, CorrTimeList(i,:))));
    end
    iCorrWindowTime = zeros(nCORR, nCorrWindowPoints);
    for i = 1:nCORR
        if CorrTimeList(i,1) < Time(1) || CorrTimeList(i,2) > Time(end)
            % Add an error message to the report
            bst_report('Error', sProcess, sInputsA, 'One correlation time window is not within the data time range.');
            % Stop the process
            return;
        end            
        single_iCorrWindow = panel_time('GetTimeIndices', Time, CorrTimeList(i,:));
        if length(single_iCorrWindow) < nCorrWindowPoints
            bst_report('Error', sProcess, sInputsA, 'Wrong time range');
            % Stop the process
            return;
        end
        iCorrWindowTime(i,:) = single_iCorrWindow(1:nCorrWindowPoints);    
    end
    
    bst_progress('start', ['Applying process: ' OPTIONS.method], 'Loading reference...', 0, length(InputsRef));

    nRef = 0; nOriRef = 0;
    nTrials = length(InputsData);
    AllRefData = cell(nCORR,1);
    RowName = {};
    if OPTIONS.RefRangeType == 1
        DelayList = ((OPTIONS.RefDelayRange(1)+CorrWindowSize):OPTIONS.RefDelayInterval:OPTIONS.RefDelayRange(2))-HalfCorrWindowSize;
        nDelay = length(DelayList);
        if nDelay==0;
            nDelay=1;
            DelayList=0;
        end
    elseif OPTIONS.RefRangeType == 2
        DelayList = OPTIONS.RefDelayRange(1):OPTIONS.RefDelayInterval:OPTIONS.RefDelayRange(2);
        nDelay = length(DelayList);
        if nDelay==0;
            nDelay=1;
            DelayList=0;
        end
    end
    RefLength = nCorrWindowPoints;
    TargetList = {};
    %AllRefData = zeros(nTrials, RefLength,length(InputsRef));
    % Reading all the input files in a big matrix
    for i = 1:length(InputsRef)
        % Read the file #i
        RefMat = in_bst(InputsRef(i).FileName, [], 0);
        
        % Check the dimensions of the recordings matrix in this file
        if mod(size(RefMat.Value,1), nTrials)~=0
            % Add an error message to the report
            bst_report('Error', sProcess, InputsData, 'One reference file has a different number of trials.');
            % Stop the process
            return;
        else
            nRefSingleFile = (size(RefMat.Value,1) / nTrials);
            refData = RefMat.Value;

        end
       
        
        hw = hann(double(nCorrWindowPoints));
        for j = 1:nRefSingleFile
            fileRefData = refData(j:nRefSingleFile:end,:);
            sRate = 1 / (RefMat.Time(2) - RefMat.Time(1));

            if strcmp(OPTIONS.measure,'esc') || strcmp(OPTIONS.measure,'cfc')
                fileRefData = amp_vec(fileRefData,OPTIONS.BandBoundsRef,sRate,OPTIONS.Width,OPTIONS.TFmethod);   
                
                if strcmp(OPTIONS.measure, 'cfc')
                    fileRefData = (fileRefData / sRate).^2;
                end
            end
            
            for m=1:nDelay
                if strcmpi(InputsRef(i).FileType,'timefreq')  
                    RowName{nRef+m} = [RefMat.RowNames{j} ' (' DelayList(m) 'ms)'];
                else
                    RowName{nRef+m} = [RefMat.Description{j} ' (' DelayList(m) 'ms)'];
                end
            end
            if isfield(RefMat,'Atlas') 
                if ~isempty(RefMat.Atlas)
                    TargetList{nOriRef+j}=RefMat.Atlas.Scouts(j); 
                end
            end
            
            if OPTIONS.RefRangeType == 1
                for m=1:nDelay
                    refTime = OPTIONS.RefDelayRange(1) + OPTIONS.RefDelayInterval*(m-1) + [0 CorrWindowSize];
                    if (RefMat.Time(1) > refTime(1)) || (RefMat.Time(end) < refTime(2))
                        % Add an error message to the report
                        bst_report('Error', sProcess, InputsData, 'The selected range for extracting reference signal is out of range.');
                        % Stop the process
                        return;
                    end
                    iRefWindow = panel_time('GetTimeIndices', Time, refTime);                   
                    iRefWindow = iRefWindow(1:nCorrWindowPoints); 
                    cRefData = fileRefData(:,iRefWindow);
                    if OPTIONS.isHann
                         for n = 1:nTrials
                             cRefData(n,:) = cRefData(n,:).*hw';
                         end  
                    end

                    %else
                    %    RefMat.Value = RefMat.Value(j:nRefSingleFile:end,RefTimePoints);
                    %    RefMatAvg = mean(RefMat.Value, 2);
                    %    fileRefData = bst_bsxfun(@minus, RefMat.Value, RefMatAvg);
                    %end   
                    AllRefData{1} = cat(3,AllRefData{1}, cRefData);
                end                
            elseif OPTIONS.RefRangeType == 2
                 
                for k=1:nCORR
                    for m=1:nDelay
                        refTime = CorrTimeList(k,:) + DelayList(m);
                        if (RefMat.Time(1) > refTime(1)) || (RefMat.Time(end) < refTime(2))
                            % Add an error message to the report
                            bst_report('Error', sProcess, InputsData, 'The selected range for extracting reference signal is out of range.');
                            % Stop the process
                            return;
                        end
                        iRefWindow = panel_time('GetTimeIndices', Time, refTime);                   
                        iRefWindow = iRefWindow(1:nCorrWindowPoints); 
                        cRefData = fileRefData(:,iRefWindow);
                        if OPTIONS.isHann
                             for n = 1:nTrials
                                 cRefData(n,:) = cRefData(n,:).*hw';
                             end  
                        end

                        %else
                        %    RefMat.Value = RefMat.Value(j:nRefSingleFile:end,RefTimePoints);
                        %    RefMatAvg = mean(RefMat.Value, 2);
                        %    fileRefData = bst_bsxfun(@minus, RefMat.Value, RefMatAvg);
                        %end   
                        AllRefData{k} = cat(3,AllRefData{k}, cRefData);
                    end
                end
            end
        end
        nRef = nRef + nRefSingleFile*nDelay;
        nOriRef = nOriRef + nRefSingleFile;
        bst_progress('inc',1);
    end
     
    RepRefIntervalPoints = round(OPTIONS.RepRefInterval / (Time(2) - Time(1)));
    
    if OPTIONS.RepRefNum > 0 && RepRefIntervalPoints > 0
        for k=1:nCORR
            % nRefTrials, RefLength, nRef
            RepRefData = zeros(nTrials,RefLength + (OPTIONS.RepRefNum-1)*RepRefIntervalPoints,nRef*OPTIONS.RepRefNum);

            % Relicate Reference Signal
            for i = 1:OPTIONS.RepRefNum           
                RepRefData(:,(i-1)*RepRefIntervalPoints + (1:RefLength),(i-1)*nRef + (1:nRef)) = AllRefData{k};
            end
            nRef = nRef*OPTIONS.RepRefNum;
            AllRefData{k} = RepRefData;
        end
    else
        OPTIONS.RepRefNum = 1;
    end
    for i=1:nCORR
        AllRefData{i} = permute(AllRefData{i},[2 3 1]); %RefLength, nRef, nRefTrials 
    end
    AllCovRef = cell(nCORR,1);
    AllRefDataCell = AllRefData;
    AllRefDataCellorig = AllRefDataCell;
    clear AllRefData;
    for m=1:nCORR
        if OPTIONS.RefRangeType == 1 && m > 1
            if nRef > 1
                AllCovRef{m}.origcov = AllCovRef{1}.origcov;
                AllCovRef{m}.originv_cov=AllCovRef{1}.originv_cov;
            end
            AllCovRef{m}.cov = AllCovRef{1}.cov;
            AllCovRef{m}.inv_cov = AllCovRef{1}.inv_cov;
            AllRefDataCellorig{m} = AllRefDataCellorig{1};
            AllRefDataCell{m} = AllRefDataCell{1};
            continue;
        end
        AllRefData = AllRefDataCell{m};
        covRef = zeros(nRef,nRef);
        for i = 1:nTrials
            singleRefTrial = squeeze(AllRefData(:,:,i))';
            singleRefTrialAvg = mean(singleRefTrial, 2);
            AllRefData(:,:,i) = bst_bsxfun(@minus, singleRefTrial, singleRefTrialAvg)';
            covRef = covRef + AllRefData(:,:,i)'*AllRefData(:,:,i);
        end
        AllRefDataCellorig{m}=AllRefData;
        
        if nRef > 1
            origAllRefData = AllRefDataCell{m};
            covRef = covRef/(nTrials*RefLength-1);
            [v,d]=eig(covRef);  
            C=v*diag(diag(d).^(-0.5));
            origcovRef = covRef;
            for i = 1:nTrials
                singleRefTrial = C'*squeeze(origAllRefData(:,:,i))';
                singleRefTrialAvg = mean(singleRefTrial, 2);
                AllRefData(:,:,i) = bst_bsxfun(@minus, singleRefTrial, singleRefTrialAvg)';
                covRef = covRef + AllRefData(:,:,i)'*AllRefData(:,:,i);
            end
            AllRefDataCell{m} = AllRefData;
%             figure;
%             subplot(2,1,1);
%             imagesc(origcovRef);colorbar;
%             subplot(2,1,2);  
%             imagesc(covRef/(nTrials*RefLength-1));colorbar;
            AllCovRef{m}.origcov=origcovRef;
            % Calculate the inverse of covRef
            AllCovRef{m}.originv_cov=inv(AllCovRef{m}.origcov);
        end
        covRef = covRef / (nTrials*RefLength-1);
        
        AllCovRef{m}.cov=covRef;
        % Calculate the inverse of covRef
        AllCovRef{m}.inv_cov=inv(covRef);
        
    end

%     if nRef == 1
%         covRef_inv = inv(covRef);
%     else
%         [U,S,V] = svd(covRef);
%         S = diag(S); % Covariance = Cm = U*S*U'
%         Si = diag(1 ./ (S + S(1) * Reg / 100)); % 1/(S^2 + lambda I)
%         covRef_inv = V*Si*U'; % U * 1/(S^2 + lambda I) * U' = Cm^(-1)   
%     end
   
    
    % Find the indices for covariance calculation
     %iMinVarTime = panel_time('GetTimeIndices', Time, MinVarTime);
    
    % Initialize the covariance matrices
     MinVarCov = zeros(nChannels, nChannels);
     CorrCovOfNumerator = zeros(nChannels, nRef ,nCORR);
     CorrCovOfDenominator = zeros(nChannels, nChannels, nCORR);
     CorrCovOfNumeratorOrig = zeros(nChannels, nRef ,nCORR);

     nTotalMinVar = zeros(nChannels, nChannels, nCORR);  
     nTotalCorrCovOfNumerator = zeros(nChannels, nRef, nCORR); 
     nTotalCorrCovOfDenominator = zeros(nChannels, nChannels, nCORR); 
     nTotalCorrCovOfNumeratorOrig = zeros(nChannels, nRef, nCORR); 

    % Extract the covariance matrix of the used channels
    NoiseCovMat  = load(file_fullpath(sChannelStudy.NoiseCov.FileName));
    NoiseCov     = NoiseCovMat.NoiseCov(iChannels,iChannels);
    
    bst_progress('start', ['Applying process: ' OPTIONS.method], 'Calculating covariance matrices...', 0, length(InputsData)*nCORR);
    
    if OPTIONS.BaselineTime(1)==OPTIONS.BaselineTime(2)
        iBaselineTime = [];
    else
        iBaselineTime = panel_time('GetTimeIndices', Time, OPTIONS.BaselineTime);
    end
    % Reading all the input files in a big matrix
    for i = 1:nTrials
        
            % Read the file #i
            DataMat = in_bst(InputsData(i).FileName, [], 0);
            % Check the dimensions of the recordings matrix in this file
            if (size(DataMat.F,1) ~= nChannels) || (size(DataMat.F,2) ~= nTime)
                % Add an error message to the report
                bst_report('Error', sProcess, InputsData, 'One file has a different number of channels or a different number of time samples.');
                % Stop the process
                return;
            end
            if strcmp(OPTIONS.measure,'esc') && OPTIONS.BandBoundsData(1)~=0 && OPTIONS.BandBoundsData(2)~=0
                DataMat.F = bp_vec(DataMat.F,OPTIONS.BandBoundsData,sRate,OPTIONS.Width,OPTIONS.TFmethod);
            end

            % Get good channels
            iGoodChan = find(DataMat.ChannelFlag == 1);            
            %for j = 1:nCORR
%             j=1;           
            if ~isempty(iBaselineTime)           
                % Remove average
                FavgActive = mean(DataMat.F(iGoodChan,iBaselineTime), 2);
                DataMinVar = bst_bsxfun(@minus, DataMat.F(iGoodChan,MinVarRangePoint), FavgActive);
            else
                DataMinVar = DataMat.F(iGoodChan,MinVarRangePoint);
            end
%             % Average baseline values
%             FavgMinVar = mean(DataMat.F(iGoodChan,MinVarRangePoint), 2);
%             % Remove average
%             DataMinVar = bst_bsxfun(@minus, DataMat.F(iGoodChan,MinVarRangePoint), FavgMinVar);
            % Compute covariance for this file
            fileMinVarCov = DataMat.nAvg .* (DataMinVar * DataMinVar');               
                
            MinVarCov(iGoodChan,iGoodChan,j) = bst_bsxfun(@plus,MinVarCov(iGoodChan,iGoodChan), fileMinVarCov);
            nTotalMinVar(iGoodChan,iGoodChan,j) = nTotalMinVar(iGoodChan,iGoodChan) + RefLength*DataMat.nAvg; 
            %end
%         end  
        
%         if ~isFreq
%             % Average baseline values
%             FavgMinVar = mean(DataMat.F(iGoodChan,iMinVarTime), 2);
%             % Remove average
%             DataMinVar = bst_bsxfun(@minus, DataMat.F(iGoodChan,iMinVarTime), FavgMinVar);
%             % Compute covariance for this file
%             fileMinVarCov = DataMat.nAvg .* (DataMinVar * DataMinVar');
%         else
%             [Cxy,freq] = crossSpectrum(DataMat.F(iGoodChan,iMinVarTime), DataMat.F(iGoodChan,iMinVarTime), Fs, SegmentLength, []);
%             fileMinVarCov = zeros(size(Cxy,1),size(Cxy,2));
%             nCxy = 0;
%             for f = 1:length(freq)
%                 if (freq(f) <= BandBounds(2)) && (freq(f) >= BandBounds(1))
%                     fileMinVarCov = fileMinVarCov + Cxy(:,:,f);
%                     nCxy = nCxy + 1;
%                 end
%             end
%             if nCxy > 0
%                 fileMinVarCov = DataMat.nAvg * fileMinVarCov / nCxy;
%             end
%         end
%         
%         %fileMinVarCov = DataMinVar * DataMinVar';
%         % Add file covariance to accumulator
%         MinVarCov(iGoodChan,iGoodChan) = MinVarCov(iGoodChan,iGoodChan) + fileMinVarCov;
%         nTotalMinVar(iGoodChan,iGoodChan) = nTotalMinVar(iGoodChan,iGoodChan) + length(iMinVarTime)*DataMat.nAvg;    
    
        for j = 1:nCORR
            AllRefData = AllRefDataCell{j};
            if strcmp(OPTIONS.measure,'cfc') || strcmp(OPTIONS.measure,'coh')
                %if ~isCrossFreq
                    [Cxy,freq,Cxx] = crossSpectrum(DataMat.F(iGoodChan,iCorrWindowTime(j,:)), AllRefData(:,:,i)', Fs, OPTIONS.SegmentLength, [], OPTIONS.WinOverlap);
                %else
                %    [Cxy,freq,Cxx] = crossSpectrumCrossFreq(DataMat.F(iGoodChan,iCorrWindowTime(j,:)), AllRefData(:,:,i)', Fs, BandBoundsRef, SegmentLength, [], WinOverlapLength);
                %end
                fileDataInCorrWindowCovD = zeros(size(Cxx,1),size(Cxx,2));
                fileDataInCorrWindowCovN = zeros(size(Cxy,1),size(Cxy,2));
                nCxy = 0;
                for f = 1:length(freq)
                    if (freq(f) <= OPTIONS.BandBoundsData(2)) && (freq(f) >= OPTIONS.BandBoundsData(1))
                        fileDataInCorrWindowCovD = fileDataInCorrWindowCovD + Cxx(:,:,f);
                        fileDataInCorrWindowCovN = fileDataInCorrWindowCovN + Cxy(:,:,f);
                        nCxy = nCxy + 1;
                    end
                end
                if nCxy > 0
                    fileDataInCorrWindowCovN = fileDataInCorrWindowCovN / nCxy;
                    fileDataInCorrWindowCovD = fileDataInCorrWindowCovD / nCxy;
                else
                    bst_report('Error', sProcess, InputsData, 'Cannot find signal within this frequency range.');  
                end
               
            else
                % Average baseline values
                DataAvg = mean(DataMat.F(iGoodChan,iCorrWindowTime(j,:)), 2);         
                % Remove average
                DataInCorrWindow = bst_bsxfun(@minus, DataMat.F(iGoodChan,iCorrWindowTime(j,:)), DataAvg);  
                % Compute covariance for this file
                fileDataInCorrWindowCovN = (DataInCorrWindow * AllRefData(:,:,i)); 
                fileDataInCorrWindowCovD = (DataInCorrWindow * DataInCorrWindow'); 
                
            end
            
            if nRef > 1
                fileDataInCorrWindowCovNorig = (DataInCorrWindow * AllRefDataCellorig{j}(:,:,i)); 
                % Add file covariance to accumulator
                CorrCovOfNumeratorOrig(iGoodChan,:,j) = bst_bsxfun(@plus,CorrCovOfNumeratorOrig(iGoodChan,:,j), fileDataInCorrWindowCovNorig);      
                nTotalCorrCovOfNumeratorOrig(iGoodChan,:,j) = nTotalCorrCovOfNumeratorOrig(iGoodChan,:,j) + RefLength;

            end
            % Add file covariance to accumulator
            CorrCovOfNumerator(iGoodChan,:,j) = bst_bsxfun(@plus, CorrCovOfNumerator(iGoodChan,:,j), fileDataInCorrWindowCovN);      
            nTotalCorrCovOfNumerator(iGoodChan,:,j) = nTotalCorrCovOfNumerator(iGoodChan,:,j) + RefLength;
       
            % Add file covariance to accumulator
            CorrCovOfDenominator(iGoodChan,iGoodChan,j) = bst_bsxfun(@plus, CorrCovOfDenominator(iGoodChan,iGoodChan,j), fileDataInCorrWindowCovD); 
            nTotalCorrCovOfDenominator(iGoodChan,iGoodChan,j) = nTotalCorrCovOfDenominator(iGoodChan,iGoodChan,j) + RefLength;
            bst_progress('inc',1);
        end
    
    end
    % Remove zeros from N matrix
    nTotalMinVar(nTotalMinVar <= 1) = 2;
    nTotalCorrCovOfNumerator(nTotalCorrCovOfNumerator <= 1) = 2;
    nTotalCorrCovOfDenominator(nTotalCorrCovOfDenominator <=1 ) = 2;
    nTotalCorrCovOfNumeratorOrig(nTotalCorrCovOfNumeratorOrig <= 1) = 2;
    
    % Divide final matrix by number of samples
    MinVarCov = MinVarCov ./ (nTotalMinVar-1);
    CorrCovOfNumerator = CorrCovOfNumerator ./ (nTotalCorrCovOfNumerator-1);
    CorrCovOfDenominator = CorrCovOfDenominator ./ (nTotalCorrCovOfDenominator-1);
    CorrCovOfNumeratorOrig = CorrCovOfNumeratorOrig ./ (nTotalCorrCovOfNumeratorOrig-1);
    
    CovCovOfNumerator = zeros(nChannels,nChannels,nCORR);
    for j = 1:nCORR
        CovCovOfNumerator(:,:,j) = CorrCovOfNumerator(:,:,j)*AllCovRef{j}.inv_cov*CorrCovOfNumerator(:,:,j)';       
    end
    
    % ===== PROCESS =====
    % Processing iChannels
    %%%% TO EDIT %%%%%
    % Get Study
  
    % Extract the covariance matrix of the used channels
    %ActiveCov = ActiveCov(iChannels,iChannels);
    %NoiseCovMat = load(file_fullpath(sChannelStudy.NoiseCov.FileName));
    %NoiseCov  = NoiseCovMat.NoiseCov(iChannels,iChannels);
    MinVarCov = MinVarCov(iChannels,iChannels,:);
    CovCovOfNumerator = CovCovOfNumerator(iChannels,iChannels,:);
    CorrCovOfNumerator = CorrCovOfNumerator(iChannels,:,:);
    CorrCovOfDenominator = CorrCovOfDenominator(iChannels,iChannels,:);
    CorrCovOfNumeratorOrig = CorrCovOfNumeratorOrig(iChannels,:,:);
   
    % Normalize the regularization parameter
    % eigValues = eig(MinVarCov);
    % Reg_alpha = Reg * max(eigValues);
    
    % Calculate the inverse of (Cm+alpha*I)
%     if isFreq
%         if strcmpi(InputsData(1).FileType,'data')
%             MinVarCov = CorrCovOfDenominator;
%             MinVarCov_inv = CorrCovOfDenominator;
%         end
%         for j=1:nCORR 
            MinVarCov_inv(:,:,1) = pInv(MinVarCov(:,:,1), OPTIONS.Reg);
%         end
%     else
%         tmp = pInv(MinVarCov, Reg);
%         MinVarCov_inv = zeros(length(iChannels),length(iChannels),nCORR);
%         for j=1:nCORR 
%             MinVarCov_inv(:,:,j) = tmp;
%         end
%     end
    
    % Get forward field
    Kernel = sHeadModel.Gain(iChannels,:);

    %Kernel(abs(Kernel(:)) < eps) = eps; % Set zero elements to strictly non-zero
    [nChannels ,nSources] = size(Kernel); % size of Gain Matrix

    % ===== CALCULATE NEURAL ACTIVITY INDEX AND FILTERS =====       
    CorrMaps = zeros(nSources/3,nCORR);
    AmpMaps = zeros(nSources/3,nCORR);
    Ori = zeros(3,nSources/3,nCORR);
    spatialFilter = zeros(nSources/3, nChannels,nCORR); % calculate filter for outputs
    
    bst_progress('start', ['Applying process: ' OPTIONS.method], 'Calculating spatial filters...', 0, nCORR*nSources/3);
    %testNum = str2double(sProcess.options.DICStest.Value);
    if OPTIONS.BeamformerType == 1
    
        for r = 1:3:nSources % calculate noise power for index
            r2 = (r+2)/3;
            Gr = Kernel(:,r:r+2);

            for j=1:nCORR


                W = pinv(Gr' * MinVarCov_inv * Gr)* Gr' * MinVarCov_inv;

%               [u,s,v] = svd(real(W * MinVarCov(:,:,j) * W'));

                [u,s,v] = svd( real(W * CovCovOfNumerator(:,:,j)* W'));
                if s(1) > 100*s(2)
                    csd = s(1,1);                       
                    %w = u(:,1);
                else
                    csd = trace(s); 
                    %w = [1 1 1]/sqrt(3);
                end


                [u,s,v] = svd( real(W * CorrCovOfDenominator(:,:,j)* W'));
                if s(1) > 100*s(2)
                    pow = s(1,1);                       
                else
                    pow = trace(s);                        
                end


                %CorrMaps(r2,j) = sqrt((w*CovCovOfNumerator(:,:,j)*w')/(w*CorrCovOfDenominator(:,:,j)*w'));
                CorrMaps(r2,j) = sqrt(csd/pow);
                spatialFilter(r2,:,j) = w;
                
                bst_progress('inc',1);

            end
        end
    else
        
        usedKernel = mat2cell(Kernel,nChannels,3*ones(1,nSources/3));
        
        Cmat = cellfun(@(x)MinVarCov_inv*x, usedKernel, 'UniformOutput', false);
        
        Dmat = cellfun(@(x,y)x'*y, usedKernel, Cmat, 'UniformOutput', false);
        
        for j=1:nCORR
            Pmat = cellfun(@(x) x'*CovCovOfNumerator(:,:,j)*x, Cmat, 'UniformOutput', false);
            Qmat = cellfun(@(x) x'*CorrCovOfDenominator(:,:,j)*x, Cmat, 'UniformOutput', false);
            invQmat = cellfun(@(x)pInv(x,0.0000000001),Qmat, 'UniformOutput', false);
        
            invQPmat = cellfun(@mtimes, invQmat, Pmat, 'UniformOutput', false);
            [eigVecmat, eigValmat] = cellfun(@eig, invQPmat, 'UniformOutput', false);
            
            if min(size(eigValmat{1,1}))>1
                eigValmat = cellfun(@diag, eigValmat, 'UniformOutput', false);
            end
            [~, imaxmat] = cellfun(@max, eigValmat, 'UniformOutput', false);
            
            Orimat = cellfun(@(x,y)x(:,y), eigVecmat, imaxmat, 'UniformOutput', false);

            Ori(:,:,j) = cell2mat(Orimat);

            bst_progress('inc',nSources/6);
            
            tmpMat = cellfun(@(x,y)x'*y*x, Orimat, Dmat, 'UniformOutput', false);
            tmpMat2 = cellfun(@mtimes, Cmat, Orimat, 'UniformOutput', false);

            sfmat = cellfun(@mrdivide, tmpMat2, tmpMat, 'UniformOutput', false);
            vcMat = cellfun(@(x)x'*NoiseCov*x, sfmat,'UniformOutput', false);

            tmpMat = cellfun(@(x)x'*CovCovOfNumerator(:,:,j)*x, sfmat);
            tmpMat2 = cellfun(@(x)x'*CorrCovOfDenominator(:,:,j)*x, sfmat);
            
            CorrMat = bsxfun(@rdivide, tmpMat, tmpMat2);
            CorrMaps(:,j) = arrayfun(@sqrt, CorrMat);
            
            tmpMat = cell2mat(vcMat);
            AmpMat =bsxfun(@rdivide, tmpMat2, tmpMat);
            AmpMaps(:,j) = AmpMat;%arrayfun(@sqrt, AmpMat);

            vcMat = cellfun(@sqrt, vcMat, 'UniformOutput', false);
            ImagingKernel = cellfun(@mrdivide, sfmat, vcMat, 'UniformOutput', false);
            spatialFilter(:,:,j) = cell2mat(ImagingKernel)';
            bst_progress('inc',nSources/6);
        end
            
            
%                 Cr = MinVarCov_inv * Gr; 
%                 Dr = Gr' * Cr;
%                 
%                 Pr = Cr' * CovCovOfNumerator(:,:,j)* Cr;
%                 Qr = Cr' * CorrCovOfDenominator(:,:,j)* Cr;
%                 %Qr = Cr' * NoiseCov* Cr;
% 
%                 % Regularize the matrix Q to avoid singular problem 
%                 invQr = pInv(Qr,0.0000000001);
% 
%                 % Compute the dipole orientation 
%                 % (the eigenvector corresponding to maximum eigenvalue of inv(Q)*P)
%                 [eigVectors,eigValues] = eig(invQr*Pr);
%                 % check whether eigValues are saved as matrix or vector
%                 if min(size(eigValues))==1
%                     [tmp, imax] = max(eigValues);
%                 else
%                     [tmp, imax] = max(diag(eigValues));
%                 end
%                 DipoleOri = eigVectors(:,imax);
%                 
%                 Ori(:,r2,j) = DipoleOri;
%                 w = Cr*DipoleOri/(DipoleOri'*Dr*DipoleOri);
%                 spatialFilter(r2,:,j) = w';
%                 CorrMaps(r2,j) = sqrt((w'*CovCovOfNumerator(:,:,j)*w)/(w'*CorrCovOfDenominator(:,:,j)*w));
%                 AmpMaps(r2,j)=sqrt((w'*CorrCovOfDenominator(:,:,j)*w)/(w'*NoiseCov*w));
% %                 if imag(CorrMaps(r2,j))~=0
% %                     disp('ERROR');
% %                 end
%             
%             bst_progress('inc',1);
        
    end
    nSources = nSources/3 ;

    nRand = 0;%OPTIONS.nRand;
    %% Correction
    if nRand > 0
        % Generate simulation sources
        bst_progress('start', ['Applying process: ' OPTIONS.method], 'Calculating Spurious Corr...', 0, nCORR*nSources);    
        totalTimePoints = size(DataMat.F,2);
        simuSignalRef = randn(totalTimePoints,nRef/nDelay,nTrials);
        simuSignalTarget = randn(totalTimePoints,nTrials,nRand);
        CorrMapsRand = zeros(nSources,nRand,nCORR);
        for j=1:nCORR
            simuMeasurementsRef = zeros(nChannels,totalTimePoints,nTrials);
            for iRef = 1:nOriRef
                r = TargetList{iRef}.Seed;
                Gr = Kernel(:,(r-1)*3+(1:3))*Ori(:,r,j);
                w = squeeze(spatialFilter(r,:,j));
                amps = sqrt(w*CorrCovOfDenominator(:,:,j)*w');

                for i=1:nTrials
                    powMat = repmat(amps*simuSignalRef(:,iRef,i)',nChannels,1);
                    simuMeasurementsRef(:,:,i) = simuMeasurementsRef(:,:,i) + powMat.*repmat(Gr,1,totalTimePoints);
                end
            end

            for r=1:nSources

                Gr = Kernel(:,(r-1)*3+(1:3))*Ori(:,r,j);
                wr = squeeze(spatialFilter(r,:,j));
                amps = sqrt(wr*CorrCovOfDenominator(:,:,j)*wr');
                for n=1:nRand
                    Cba = zeros(nRef,1);
                    Cbb = zeros(nRef,nRef);
                    Caa = 0;
                    for i=1:nTrials
                        powMat = repmat(amps*simuSignalTarget(:,i,n)',nChannels,1);
                        simuMeasurementsAll = simuMeasurementsRef(:,:,i) + powMat.*repmat(Gr,1,totalTimePoints);
                        ind = 1;
                        refSignals = zeros(nRef,nCorrWindowPoints);
                        for iRef = 1:nOriRef
                            w = squeeze(spatialFilter(TargetList{iRef}.Seed,:,j));               
                            tmp=w*simuMeasurementsAll;
                            for m=1:nDelay
                                refTime = CorrTimeList(j,:) + DelayList(m);
                                iRefWindow = panel_time('GetTimeIndices', Time, refTime);                    
                                refSignals(ind,:) = tmp(iRefWindow(1:nCorrWindowPoints));
                                ind = ind+1;
                                if OPTIONS.isHann                                
                                    refSignals(ind,:) = refSignals(ind,:).*hw';                                
                                end
                            end
                        end

                        tmp = wr*simuMeasurementsAll(:,iCorrWindowTime(j,:));
                        Cba = Cba + refSignals*tmp';
                        Cbb = Cbb + refSignals*refSignals';
                        Caa = Caa + tmp*tmp';              
                    end
                    Cba = Cba / (nTrials*nCorrWindowPoints-1);
                    Cbb = Cbb / (nTrials*nCorrWindowPoints-1);
                    Caa = Caa / (nTrials*nCorrWindowPoints-1);
                    CorrMapsRand(r,n,j) = sqrt(Cba'*inv(Cbb)*Cba/Caa);

                end
                bst_progress('inc',1);
            end
        end
    end
    
    %% ===== Calculate weights =====
    if nDelay > 1
        bst_progress('start', ['Applying process: ' OPTIONS.method], 'Calculating weights ...', 0, nSources);

        CorrWeight = zeros(nSources,nCORR,nDelay,nRef/nDelay);
        CorrWeight2 = zeros(nSources,nCORR,nDelay,nRef/nDelay);

        for r = 1:1:nSources
            for j=1:nCORR
                w = squeeze(spatialFilter(r,:,j));
                Cab = w*CorrCovOfNumeratorOrig(:,:,j);
                
                for k = 1:nRef                  
                    f(k) = (w*squeeze(CorrCovOfNumeratorOrig(:,k,j)))/(sqrt(w*CorrCovOfDenominator(:,:,j)*w')*sqrt(AllCovRef{j}.origcov(k,k)));
                end
                
                f2 = Cab*AllCovRef{j}.originv_cov;
                f2 = f2/norm(f2)*CorrMaps(r,j);
                ind = 0;
                for k = 1:nDelay
                    for i = 1:nRef/nDelay
                        ind = ind + 1;
                        CorrWeight(r,j,k,i) = f(ind);
                        CorrWeight2(r,j,k,i) = f2(ind);
                        %f(ind) = Cab(ind)/AllCovRef{j}.origcov(ind,ind);
                    end
                end
                
            end
            bst_progress('inc',1);
        end
        %CorrWeight = CorrWeight/max(abs(CorrWeight(:)));
    end

    
    %% ===== ASSIGN MAPS AND KERNELS =====
    bst_progress('start', ['Applying process: ' OPTIONS.method], 'Saving results...', 0, 3+nCORR*3);
    timestring = sprintf('%d_%d',round(OPTIONS.CORRrange(1)*1000),round(OPTIONS.CORRrange(2)*1000));
    
        % == Save the spatial filter as ImagingKernel ==
%     for i=1:nCORR
%         sp_window_time_string = sprintf('%.1f-%.1f', CorrTimeList(i,1)*1000, CorrTimeList(i,2)*1000);
%         sp_window_time_string = [];
    saveSF = 1;
    if saveSF == 1
        % Create a new data file structure
        ResultsMat2 = db_template('resultsmat');
        if nCORR > 1
            [Cval,Ind]=max(CorrMaps');
            clear Cval;
            for i=1:nSources
                spatialFilter(i,:,1) = spatialFilter(i,:,Ind(i));
                Ori(:,i,1) = Ori(:,i,Ind(i));
            end
            spatialFilter = squeeze(spatialFilter(:,:,1));
            Ori = squeeze(Ori(:,:,1));
        end
        ResultsMat2.ImagingKernel = spatialFilter;
        ResultsMat2.ImageGridAmp  = [];
        ResultsMat2.EstimatedGridOrient = Ori';
        ResultsMat2.nComponents = 1; 

        if strcmp(OPTIONS.measure,'cfc') || strcmp(OPTIONS.measure,'coh')
            if ~isCrossFreq               
                if BeamformerType == 1
                    ResultsMat2.Comment   = [result_comment OPTIONS.method ': spatial filter |' timestring  'ms(' num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) 'Hz)' ];
                else
                    ResultsMat2.Comment   = [result_comment OPTIONS.method ': spatial filter |' timestring  'ms(' num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) 'Hz)' ];
                end
            else
                if BeamformerType == 1
                    ResultsMat2.Comment   = [result_comment OPTIONS.method ': spatial filter |' timestring  'ms('  num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) ','  num2str(BandBoundsRef(1))  '-' num2str(BandBoundsRef(2)) 'Hz)' ];
                else
                    ResultsMat2.Comment   = [result_comment OPTIONS.method ': spatial filter |' timestring  'ms('  num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) ','  num2str(BandBoundsRef(1))  '-' num2str(BandBoundsRef(2)) 'Hz)' ];
                end
            end
              
        else
             if OPTIONS.BeamformerType == 1
                ResultsMat2.Comment   = [result_comment OPTIONS.method ': spatial filter |' timestring  'ms' ];
            else
                ResultsMat2.Comment   = [result_comment OPTIONS.method ': spatial filter |' timestring  'ms' ];
            end
        end

        ResultsMat2.Function      = 'SILSCspatialfilter';
        ResultsMat2.Time          = [];           % Leave it empty if using ImagingKernel
        ResultsMat2.DataFile      = [];
        ResultsMat2.HeadModelFile = HeadModelFile;
        ResultsMat2.HeadModelType = sHeadModel.HeadModelType;
        ResultsMat2.ChannelFlag   = [];
        ResultsMat2.GoodChannel   = iChannels;
        if strcmp(sHeadModel.HeadModelType,'volume')     
            ResultsMat2.GridLoc       = sHeadModel.GridLoc;
            ResultsMat2.SurfaceFile   = sHeadModel.SurfaceFile;
        else
            ResultsMat2.SurfaceFile   = sHeadModel.SurfaceFile;
        end
        %ResultsMat2.GridLoc       = [];%sHeadModel.GridLoc;
        % === SHARED ==
        % Get the output study (pick the one from the first file)
        iStudy = iChannelStudy;
        %nOutputFiles = length(OutputFiles);
        % Create a default output filename 
        Outs = bst_process('GetNewFilename', fileparts(InputsData(1).ChannelFile), 'results_beam_KERNEL');
        
        % Save on disk
        save(Outs, '-struct', 'ResultsMat2');
%         bst_save(OutputFiles{nOutputFiles+1}, ResultsMat2,'v6');
        % Register in database
        db_add_data(iStudy, Outs, ResultsMat2);
        bst_progress('inc',1);
    end
%     end
    
    
    if nCORR == 1        
%         CORRRangePoints = panel_time('GetTimeIndices', Time, CORRrange);
%        bst_progress('start', 'Applying process: SISSC', 'Interpolating correlation dynamics...', 0, length(CORRRangePoints));
%         ImageGridAmp = zeros(nSources,nTime);
%         
%         for i = CORRRangePoints(1):CORRRangePoints(end) 
%             ImageGridAmp(:,i) = CorrMaps(:,1);
%             bst_progress('inc',1);
%         end

        TFMat = db_template('timefreqmat');
%         timestring2 = sprintf('(ws:%d%s)',round(CorrWindowSize*1000),'ms');
        TFMat.Time          = [OPTIONS.CORRrange(1), OPTIONS.CORRrange(end)]; 
        %TFMat.TimeBands         = {'corr',[num2str(CORRrange(1)) ',' num2str(CORRrange(end))], 'mean'};           % Leave it empty if using ImagingKernel
        TFMat.TF            = [CorrMaps CorrMaps];
        
        
        [mv, mr] = max(CorrMaps);
        
        sSubject = bst_get('Subject', 0);
        sMRI = load(file_fullpath(sSubject.Anatomy.FileName));
        lc = cs_convert(sMRI, 'scs', 'mni', sHeadModel.GridLoc(mr,:))*1000;
        disp(['BIPAC> Maximum Peak Location: ' num2str(round(lc)) ]);
        disp(['BIPAC> Maximum Peak Value: ' num2str(mv) ]);
        
    else % Interpolate the correlation maps to have the same temporal resolution as the data
        bst_progress('start', ['Applying process: ' OPTIONS.method], ['Interpolating ' OPTIONS.measure ' dynamics...'], 0, nCORR+1);

        ImageGridAmpOriginal = CorrMaps;
        ImageGridAmp = zeros(nSources,nTime);
        for i=1:(nCORR-1)
            InterpolateTimeWindow = OPTIONS.CORRrange(1)+HalfCorrWindowSize+[(i-1)*OPTIONS.CORRTResolu i*OPTIONS.CORRTResolu];
            iInterpolateTime = panel_time('GetTimeIndices', Time, InterpolateTimeWindow);
            nInterpolateTime = length(iInterpolateTime);
            
            ImageGridAmp(:,iInterpolateTime(1)) = ImageGridAmpOriginal(:,i); 
            ImageGridAmp(:,iInterpolateTime(end)) = ImageGridAmpOriginal(:,i+1); 
            
            for j=2:(nInterpolateTime-1)       
                InterpolatePercentage = (j-1)/(nInterpolateTime-1);
                ImageGridAmp(:,iInterpolateTime(j)) = ImageGridAmpOriginal(:,i+1)*InterpolatePercentage + ImageGridAmpOriginal(:,i)*(1-InterpolatePercentage);
            end
            bst_progress('inc',1);
        end
            
        
        TFMat = db_template('timefreqmat');
%         timestring2 = sprintf('(ws:%d%s,tr:%.1f%s)',round(CorrWindowSize*1000),'ms',OPTIONS.CORRTResolu*1000,'ms');     
        TFMat.Time          = Time;%((OPTIONS.CORRrange(1)+HalfCorrWindowSize):OPTIONS.CORRTResolu:(OPTIONS.CORRrange(end)-HalfCorrWindowSize))-(OPTIONS.RefDelayRange(1)+HalfCorrWindowSize); 
        %TFMat.Time          = CorrTimeList(:,1);
        TFMat.TimeBands     = {};
        TFMat.TF            = ImageGridAmp;
             
        [mv, mr] = max(max(ImageGridAmp'));
        [~, mt] = max(max(ImageGridAmp));

        sSubject = bst_get('Subject', 0);
        sMRI = load(file_fullpath(sSubject.Anatomy.FileName));
        lc = cs_convert(sMRI, 'scs', 'mni', sHeadModel.GridLoc(mr,:))*1000;
        disp(['BIPAC> Maximum Peak Location: ' num2str(round(lc)) ]);
        disp(['BIPAC> Maximum Peak Value: ' num2str(mv) ]);
        disp(['BIPAC> Maximum Peak Time: ' num2str(TFMat.Time(mt)) ]);
        
    end
    
    
%     if isFreq
%         if ~isCrossFreq
%             TFMat.Method        = 'cohere';
%             TFMat.Freqs         = {'custum',[num2str(OPTIONS.BandBoundsData(1)) ',' num2str(OPTIONS.BandBoundsData(2))],'mean'};
%             if BeamformerType == 1                
%                 TFMat.Comment   = [result_comment 'DICS: cohere|' timestring  'ms(' num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) 'Hz)' timestring2];
%             else   
%                 TFMat.Comment   = [result_comment 'SILSC: cohere|' timestring  'ms(' num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) 'Hz)' timestring2];
%             end
%         else
%             TFMat.Method        = 'crossfreq coupling';
%             TFMat.Freqs         = {'custum',[num2str(BandBoundsData(1)) ',' num2str(BandBoundsData(2))],'mean'};
%             if BeamformerType == 1                
%                 TFMat.Comment   = [result_comment 'DICS: crxfreqcoup|' timestring  'ms(' num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) ','  num2str(BandBoundsRef(1))  '-' num2str(BandBoundsRef(2)) 'Hz)' timestring2];
%             else   
%                 TFMat.Comment   = [result_comment 'SILSC: crxfreqcoup|' timestring  'ms(' num2str(BandBoundsData(1))  '-' num2str(BandBoundsData(2)) ','  num2str(BandBoundsRef(1))  '-' num2str(BandBoundsRef(2)) 'Hz)' timestring2];
%             end
%         end
%     else
%         
%         TFMat.Freqs         = 0;
%     end
    TFMat.OPTIONS       = OPTIONS;
    TFMat.Method        = 'corr';
    TFMat.TFMethod        = OPTIONS.TFmethod;
    TFMat.Comment   = [result_comment OPTIONS.method ': ' upper(OPTIONS.measure) '(' timestring  'ms)'];
%     TFMat.TF = CorrMaps;
    TFMat.DataFile      = [];
    TFMat.DataType      = 'results';
    TFMat.Measure       = 'other';
    TFMat.RefRowNames   = RowName;
    TFMat.RowNames      = 1:nSources;
              % Leave it empty if using ImagingKernel
    
    TFMat.nAvg          = length(InputsData);
    if strcmp(sHeadModel.HeadModelType,'volume')     
        TFMat.GridLoc       = sHeadModel.GridLoc;
        TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
    else
        TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
    end
%     TFMat.OPTIONS.EstimatedOri = Ori';
    TFMat.OPTIONS.SpatialFilterFileName = Outs;
    TFMat.HeadModelFile = HeadModelFile;
    TFMat.HeadModelType = sHeadModel.HeadModelType;
    % === NOT SHARED ===
    % Get the output study (pick the one from the first file)       
    [StudyContent,iStudy]=bst_get('StudyWithCondition',fileparts(InputsData(1).FileName));
 
    % Create a default output filename 
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(InputsData(1).FileName), 'timefreq_connect1_beam');

    % Save on disk
    bst_save(OutputFiles{1}, TFMat,'v6');
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, TFMat);
    bst_progress('inc',1);
    clear TFMat;
    

    
    % == Save the baseline corr
    
    if nRand > 1
        for i=1:nCORR
            timestring = sprintf('%d_%d',round(CorrTimeList(i,1)*1000),round(CorrTimeList(i,2)*1000));     
            TFMat = db_template('timefreqmat');
            TFMat.Time          = 1:nRand; 
            TFMat.TF            = CorrMapsRand(:,:,i);
            TFMat.OPTIONS       = OPTIONS;
            TFMat.Method        = OPTIONS.measure;
            TFMat.TFMethod        = OPTIONS.TFmethod;
            TFMat.Comment   = [result_comment OPTIONS.method ': ' upper(OPTIONS.measure) '(rand, ' timestring  'ms)'];
            TFMat.DataFile      = [];
            TFMat.DataType      = 'results';
            TFMat.Measure       = 'other';
            TFMat.RefRowNames   = RowName;
            TFMat.RowNames      = 1:nSources;
            % Leave it empty if using ImagingKernel   
            TFMat.nAvg          = length(InputsData);
            if strcmp(sHeadModel.HeadModelType,'volume')     
                TFMat.GridLoc       = sHeadModel.GridLoc;
                TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
            else
                TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
            end
            % === NOT SHARED ===
            % Get the output study (pick the one from the first file)

            [StudyContent,iStudy]=bst_get('StudyWithCondition',fileparts(sInputsData(1).FileName));

            % Create a default output filename 
            OutputFiles{i+1} = bst_process('GetNewFilename', fileparts(InputsData(1).FileName), 'timefreq_connect1_beam');

            % Save on disk
            bst_save(OutputFiles{i+1}, TFMat,'v6');
            % Register in database
            db_add_data(iStudy, OutputFiles{i+1}, TFMat);
            bst_progress('inc',1);
            clear TFMat
        end
    end
    
    
        % == Save the contrast
    timestring = sprintf('%d_%d',round(OPTIONS.CORRrange(1)*1000),round(OPTIONS.CORRrange(2)*1000));
    if nCORR == 1        
%         CORRRangePoints = panel_time('GetTimeIndices', Time, CORRrange);
%        bst_progress('start', 'Applying process: SISSC', 'Interpolating correlation dynamics...', 0, length(CORRRangePoints));
%         ImageGridAmp = zeros(nSources,nTime);
%         
%         for i = CORRRangePoints(1):CORRRangePoints(end) 
%             ImageGridAmp(:,i) = CorrMaps(:,1);
%             bst_progress('inc',1);
%         end
        % ===== SAVE THE RESULTS =====
        bst_progress('start', ['Applying process: ' OPTIONS.method], 'Saving results...', 0, 2);
        
        TFMat = db_template('timefreqmat');
%         timestring2 = sprintf('(ws:%d%s)',round(CorrWindowSize*1000),'ms');
        TFMat.Time          = [OPTIONS.CORRrange(1), OPTIONS.CORRrange(end)]; 
        %TFMat.TimeBands         = {'corr',[num2str(CORRrange(1)) ',' num2str(CORRrange(end))], 'mean'};           % Leave it empty if using ImagingKernel
        TFMat.TF            = AmpMaps;
    else % Interpolate the correlation maps to have the same temporal resolution as the data
        bst_progress('inc',1);
        % Create a new data file structure
        bst_progress('start', ['Applying process: ' OPTIONS.method], 'Saving results...', 0, 2);
        ImageGridAmpOriginal = AmpMaps;
        ImageGridAmp = zeros(nSources,nTime);
        for i=1:(nCORR-1)
            InterpolateTimeWindow = OPTIONS.CORRrange(1)+HalfCorrWindowSize+[(i-1)*OPTIONS.CORRTResolu i*OPTIONS.CORRTResolu];
            iInterpolateTime = panel_time('GetTimeIndices', Time, InterpolateTimeWindow);
            nInterpolateTime = length(iInterpolateTime);
            
            ImageGridAmp(:,iInterpolateTime(1)) = ImageGridAmpOriginal(:,i); 
            ImageGridAmp(:,iInterpolateTime(end)) = ImageGridAmpOriginal(:,i+1); 
            
            for j=2:(nInterpolateTime-1)       
                InterpolatePercentage = (j-1)/(nInterpolateTime-1);
                ImageGridAmp(:,iInterpolateTime(j)) = ImageGridAmpOriginal(:,i+1)*InterpolatePercentage + ImageGridAmpOriginal(:,i)*(1-InterpolatePercentage);
            end
            bst_progress('inc',1);
        end
        TFMat = db_template('timefreqmat');
%         timestring2 = sprintf('(ws:%d%s,tr:%.1f%s)',round(CorrWindowSize*1000),'ms',OPTIONS.CORRTResolu*1000,'ms');     
        TFMat.Time          = Time;%((OPTIONS.CORRrange(1)+HalfCorrWindowSize):OPTIONS.CORRTResolu:(OPTIONS.CORRrange(end)-HalfCorrWindowSize))-(OPTIONS.RefDelayRange(1)+HalfCorrWindowSize); 
        TFMat.TimeBands     = {};
        TFMat.TF            = ImageGridAmp;
        
    end
    %TFMat.TF            = AmpMaps;
    TFMat.OPTIONS       = OPTIONS;
    TFMat.Method        = OPTIONS.measure;
    TFMat.TFMethod        = OPTIONS.TFmethod;
    TFMat.Comment   = [result_comment OPTIONS.method ': Fmap(' timestring  'ms)'];
    TFMat.DataFile      = [];
    TFMat.DataType      = 'results';
    TFMat.Measure       = 'other';
    TFMat.RefRowNames   = RowName;
    TFMat.RowNames      = 1:nSources;
              % Leave it empty if using ImagingKernel
    
    TFMat.nAvg          = length(InputsData);
    if strcmp(sHeadModel.HeadModelType,'volume')     
        TFMat.GridLoc       = sHeadModel.GridLoc;
        TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
    else
        TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
    end
    TFMat.HeadModelFile = HeadModelFile;
    TFMat.HeadModelType = sHeadModel.HeadModelType;
    % === NOT SHARED ===
    % Get the output study (pick the one from the first file)
    
    nOutputFiles = length(OutputFiles);
    % Create a default output filename 
    OutputFiles{nOutputFiles+1} = bst_process('GetNewFilename', fileparts(InputsData(1).FileName), 'timefreq_connect1_beam');

    % Save on disk
    bst_save(OutputFiles{nOutputFiles+1}, TFMat,'v6');
    % Register in database
    db_add_data(iStudy, OutputFiles{nOutputFiles+1}, TFMat);
    bst_progress('inc',1);
    clear TFMat;
            
    %=== save weight===
    if nDelay > 1
        cwcomment = 'weight';
        for cw = 1:2
            if cw==2
                CorrWeight = CorrWeight2;
                cwcomment = 'normalizedweight';
            end
            CorrWeight = permute(CorrWeight,[1 3 2 4]);           
            nOutputFiles = length(OutputFiles);
            
            for k = 1:nCORR%nDelay
                timestring = sprintf('%d_%d',round(CorrTimeList(k,1)*1000),round(CorrTimeList(k,2)*1000));   
                SingleCorrWeight = CorrWeight(:,:,k,1);
    %             TFMat.Time = DelayList;
                %timestring = sprintf('%d_%d',round(OPTIONS.CORRrange(1)*1000),round(OPTIONS.CORRrange(2)*1000));
%                 if nDelay == 1        
                    % ===== SAVE THE RESULTS =====
                    bst_progress('start', ['Applying process: ' OPTIONS.method ], 'Saving results...', 0, 2);

                    TFMat = db_template('timefreqmat');
    %                 timestring2 = sprintf('(ws:%d%s)',round(CorrWindowSize*1000),'ms');
    %                 TFMat.Time          = [OPTIONS.CORRrange(1), OPTIONS.CORRrange(end)]; 
                    %TFMat.TimeBands         = {'corr',[num2str(CORRrange(1)) ',' num2str(CORRrange(end))], 'mean'};           % Leave it empty if using ImagingKernel

%                 else % Interpolate the correlation maps to have the same temporal resolution as the data
%                     bst_progress('start', ['Applying process: ' OPTIONS.method ], 'Interpolating correlation dynamics...', 0, nCORR+1);
%                     bst_progress('inc',1);
%                     % Create a new data file structure
%                     bst_progress('start', ['Applying process: ' OPTIONS.method ], 'Saving weights...', 0, 2);
% 
%                     TFMat = db_template('timefreqmat');
%     %                 timestring2 = sprintf('(ws:%d%s,tr:%.1f%s)',round(CorrWindowSize*1000),'ms',CORRTResolu*1000,'ms');     
%     %                 TFMat.Time          = (OPTIONS.CORRrange(1)+HalfCorrWindowSize):CORRTResolu:(OPTIONS.CORRrange(end)-HalfCorrWindowSize); 
%                 end
                TFMat.Time = DelayList;
                if strcmp(OPTIONS.measure,'cfc') || strcmp(OPTIONS.measure,'coh')               
    %               TFMat.Freqs         = 1:nDelay                              
                    TFMat.Comment   = [result_comment OPTIONS.method ':' cwcomment ' (' timestring  'ms,' num2str(BandBounds(1))  '-' num2str(BandBounds(2)) 'Hz)' ];                                 
                else
                    TFMat.Comment   = [result_comment OPTIONS.method ':' cwcomment ' (' timestring  'ms)' ];
                end
                TFMat.TF = SingleCorrWeight;
                TFMat.DataFile      = [];
                TFMat.Method        = 'corr';
                TFMat.DataType      = 'results';
                TFMat.Measure       = 'other';
                TFMat.RefRowNames   = RowName;
                TFMat.RowNames      = 1:nSources;
                          % Leave it empty if using ImagingKernel

                TFMat.nAvg          = length(InputsData);

                if strcmp(sHeadModel.HeadModelType,'volume')     
                    TFMat.GridLoc       = sHeadModel.GridLoc;
                else
                    TFMat.SurfaceFile   = sHeadModel.SurfaceFile;
                end
                % === NOT SHARED ===
                
                % Get the output study (pick the one from the first file)
                [StudyContent,iStudy]=bst_get('StudyWithCondition',fileparts(InputsData(1).FileName));

                % Create a default output filename 
                OutputFiles{nOutputFiles+k} = bst_process('GetNewFilename', fileparts(InputsData(1).FileName), 'timefreq_connect1_BeamWeight');

                % Save on disk
                bst_save(OutputFiles{nOutputFiles+k}, TFMat,'v6');
                % Register in database
                db_add_data(iStudy, OutputFiles{nOutputFiles+k}, TFMat);
                bst_progress('inc',1);
                clear TFMat
            end
        end
    end
    bst_progress('stop');
end
%% ===== TRUNCATED PSEUDO-INVERSE =====
function X = pInv(A,Reg)
    % Inverse of 3x3 GCG' in unconstrained beamformers.
    % Since most head models have rank 2 at each vertex, we cut all the fat and
    % just take a rank 2 inverse of all the 3x3 matrices
%     [U,S,V] = svd(A);
%     Si = diag(1 ./ (S + S(1) * Reg / 100)); % 1/(S^2 + lambda I)
%     X = V*diag(Si)*U';

     eigValues = eig(A);
     Reg_alpha = Reg / 100 * max(eigValues);
     X = inv(A+Reg_alpha*eye(size(A,1)));
end
function TF = hilberData(F,BandBounds) 
    % Band-pass filter in one frequency band
    Fband = process_bandpass('Compute', F, Fs, BandBounds(1), BandBounds(2), [], 1);
    % Apply Hilbert transform
    TF = abs(hilbert(Fband')');
end
function [spectra,freq,spectraX] = crossSpectrum(X, Y, Fs, segmentLengthTime, MaxFreqRes, WinOverlap)
    % Default options
    if (nargin < 5) || isempty(MaxFreqRes)
        MaxFreqRes = 1;
    end
    if (nargin < 4) || isempty(segmentLengthTime)
        segmentLengthTime = 128;
    end
    segmentLength = floor(segmentLengthTime * Fs);
    overlap = floor(WinOverlap * segmentLength);
    
    % Signal properties
    nX = size(X, 1); 
    nY = size(Y, 1);
    nTimes = size(X, 2);
    
    % Segment indices - discard final timepoints
    overlapLength = floor(overlap * segmentLength);
    nSegments = floor((nTimes-overlapLength) / (segmentLength-overlapLength));
    partialLength = segmentLength - overlapLength;
    segmentStart = partialLength * (0:(nSegments-1)) + 1;
    segmentEnd = segmentStart + (segmentLength-1);
    if (segmentEnd(end) > nTimes)
        segmentStart(end) = [];
        segmentEnd(end) = [];
    end
    segmentIndices = [segmentStart; segmentEnd];
    % Maximum default resolution: 1 Hz
    if ~isempty(MaxFreqRes) && (nTimes > round(Fs / MaxFreqRes))
        nFFT = 2^nextpow2( round(Fs / MaxFreqRes) );
    % Use the default for FFT 
    else
        nFFT = 2^nextpow2( nTimes );
    end
    % Output frequencies
    freq = Fs/2*linspace(0, 1, nFFT/2 + 1)';
    freq(end) = [];

    % Frequency smoother (represented as time-domain multiplication)
    smoother = window('parzenwin', segmentLength) .* tukeywin(segmentLength, 0);
    smoother = smoother / sqrt(sum(smoother.^2));
    
    %% ===== VERSION 1: for loops, full matrix =====

    isCalcAuto = isequal(X,Y);
    % Initialize variables
    spectra = zeros(nX, nY, length(freq));

    % Cross-spectrum
    if ~isCalcAuto
        spectraX = zeros(nX, nX, length(freq));
%         spectraY = zeros(nY, nY, length(freq));
        for iSeg = 1:nSegments
            % Get time indices for this segment
            iTime = segmentIndices(1,iSeg):segmentIndices(2,iSeg);
            % Frequency domain spectrum after smoothing and tapering
            fourierX = abs(fft(bst_bsxfun(@times, X(:,iTime), smoother'), nFFT, 2));
            fourierY = abs(fft(bst_bsxfun(@times, Y(:,iTime), smoother'), nFFT, 2));
            % Calculate for each frequency: fourierX * fourierY'
            for f = 1:length(freq)
                spectra(:,:,f) = spectra(:,:,f) + fourierX(:,f) * fourierY(:,f)';
                spectraX(:,:,f) = spectraX(:,:,f) + fourierX(:,f) * fourierX(:,f)';
%                 spectraY(:,:,f) = spectraY(:,:,f) + fourierY(:,f) * fourierY(:,f)';
            end             
        end
        spectraX = spectraX / (nSegments * Fs);
    else
        for iSeg = 1:nSegments
            % Get time indices for this segment
            iTime = segmentIndices(1,iSeg):segmentIndices(2,iSeg);
            % Frequency domain spectrum after smoothing and tapering
            fourierX = fft(bst_bsxfun(@times, X(:,iTime), smoother'), nFFT, 2);
            % Calculate for each frequency: fourierX * fourierY'
            for f = 1:length(freq)
                spectra(:,:,f) = spectra(:,:,f) + fourierX(:,f) * fourierX(:,f)';
            end
        end
    end

    % Normalize for segments and sampling rate
    spectra = spectra / (nSegments * Fs);
    

    % [NxN]: Auto spectrum for X is contained within cross-spectral estimation
%         bst_progress('set', round(waitStart + 0.9 * waitMax));

end
function [spectra,freq,spectraX] = crossSpectrumCrossFreq(X, Y, Fs, BandBoundsY, segmentLengthTime, MaxFreqRes, WinOverlapLength)
    % Default options
    if (nargin < 6) || isempty(MaxFreqRes)
        MaxFreqRes = 1;
    end
    if (nargin < 5) || isempty(segmentLengthTime) 
        segmentLengthTime = 128;
    end
    segmentLength = floor(segmentLengthTime * Fs);
    %overlapLength = floor(WinOverlapLengthTime * Fs);
    overlap = WinOverlapLength;
    % Signal properties
    nX = size(X, 1); 
    nY = size(Y, 1);
    nTimes = size(X, 2);
    
    % Segment indices - discard final timepoints
    overlapLength = floor(overlap * segmentLength);
    nSegments = floor((nTimes-overlapLength) / (segmentLength-overlapLength));
    partialLength = segmentLength - overlapLength;
    segmentStart = partialLength * (0:(nSegments-1)) + 1;
    segmentEnd = segmentStart + (segmentLength-1);
    if (segmentEnd(end) > nTimes)
        segmentStart(end) = [];
        segmentEnd(end) = [];
    end
    segmentIndices = [segmentStart; segmentEnd];

    % Maximum default resolution: 1 Hz
    if ~isempty(MaxFreqRes) && (nTimes > round(Fs / MaxFreqRes))
        nFFT = 2^nextpow2( round(Fs / MaxFreqRes) );
    % Use the default for FFT 
    else
        nFFT = 2^nextpow2( nTimes );
    end
    % Output frequencies
    freq = Fs/2*linspace(0, 1, nFFT/2 + 1)';
    freq(end) = [];

    % Frequency smoother (represented as time-domain multiplication)
    smoother = window('parzenwin', segmentLength) .* tukeywin(segmentLength, 0);
    smoother = smoother / sqrt(sum(smoother.^2));
    
    %% ===== VERSION 1: for loops, full matrix =====

    isCalcAuto = isequal(X,Y);
    % Initialize variables
    spectra = zeros(nX, nY, length(freq));

    % Cross-spectrum
    if ~isCalcAuto
        spectraX = zeros(nX, nX, length(freq));
        %spectraY = zeros(nY, nY, length(freq));
        for iSeg = 1:nSegments
            % Get time indices for this segment
            iTime = segmentIndices(1,iSeg):segmentIndices(2,iSeg);
            % Frequency domain spectrum after smoothing and tapering
            fourierX = fft(bst_bsxfun(@times, X(:,iTime), smoother'), nFFT, 2);
            fourierY = fft(bst_bsxfun(@times, Y(:,iTime), smoother'), nFFT, 2);
            % Calculate for each frequency: fourierX * fourierY'
            for f = 1:length(freq)
                spectraOne = zeros(nX, nY);
                nF = 0;
                for f1 = 1:length(freq)
                    if freq(f1) >= BandBoundsY(1) && freq(f1) <= BandBoundsY(2)
                        spectraOne = spectraOne + fourierX(:,f) * fourierY(:,f1)';
                        nF = nF + 1;
                    end
                end
                if nF > 0
                    spectra(:,:,f) = spectra(:,:,f) + spectraOne / nF;  
                end
                spectraX(:,:,f) = spectraX(:,:,f) + fourierX(:,f) * fourierX(:,f)';
                %spectraY(:,:,f) = spectraY(:,:,f) + fourierY(:,f) * fourierY(:,f)';
            end             
        end
        spectraX = spectraX / (nSegments * Fs);
        %spectraY = spectraY / (nSegments * Fs);
    else
        for iSeg = 1:nSegments
            % Get time indices for this segment
            iTime = segmentIndices(1,iSeg):segmentIndices(2,iSeg);
            % Frequency domain spectrum after smoothing and tapering
            fourierX = fft(bst_bsxfun(@times, X(:,iTime), smoother'), nFFT, 2);
            % Calculate for each frequency: fourierX * fourierY'
            for f = 1:length(freq)
                spectra(:,:,f) = spectra(:,:,f) + fourierX(:,f) * fourierX(:,f)';
            end
        end
    end

    % Normalize for segments and sampling rate
    spectra = spectra / (nSegments * Fs);
    

    % [NxN]: Auto spectrum for X is contained within cross-spectral estimation
%         bst_progress('set', round(waitStart + 0.9 * waitMax));

end
% function TF = amp_vec(F,BandBounds,Fs,width) 
%     lower_bin = BandBounds(1);
%     upper_bin = BandBounds(2);
%     N = size(F,1);
%     TF = [];
%     for i = 1:N
%         F1 = ampvec((lower_bin + floor((upper_bin- lower_bin)/2)),F(i,:)', Fs, width);
%         if isempty(TF)
%             TF = zeros(N,length(F1));
%         end
%         TF(i,:) = F1';
%     end
% %     % Band-pass filter in one frequency band
% %     Fband = process_bandpass('Compute', F, BandBounds(1), BandBounds(2), [], 1);
% %     % Apply Hilbert transform
% %     TF = abs(hilbert(Fband')');
% end
% function TF = bp_vec(F,BandBounds,Fs,width)
%     lower_bin = BandBounds(1);
%     upper_bin = BandBounds(2);
%     N = size(F,1);
%     TF = [];
%     for i = 1:N
%         F1 = bpvec((lower_bin + floor((upper_bin- lower_bin)/2)),F(i,:)', Fs, width);
%         if isempty(TF)
%             TF = zeros(N,length(F1));
%         end
%         TF(i,:) = F1';
%     end
% %     % Band-pass filter in one frequency band
% %     Fband = process_bandpass('Compute', F, BandBounds(1), BandBounds(2), [], 1);
%     
% end
% function TF = ph_vec(F,BandBounds,Fs,width) 
%     lower_bin = BandBounds(1);
%     upper_bin = BandBounds(2);
%     N = size(F,1);
%     TF = [];
%     for i = 1:N
%         F1 = phasevec((lower_bin + floor((upper_bin- lower_bin)/2)),F(i,:)', Fs, width);
%         if isempty(TF)
%             TF = zeros(N,length(F1));
%         end
%         TF(i,:) = F1';
%     end
% %     % Band-pass filter in one frequency band
% %     Fband = process_bandpass('Compute', F, BandBounds(1), BandBounds(2), [], 1);
% %     % Apply Hilbert transform
% %     TF = angle(hilbert(Fband')');
% end
function TF = ana_vec(F,BandBounds,Fs,width,SegmentOverlap) 
    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    N = size(F,1);
    TF = [];
    w = floor(width*Fs*(1/lower_bin));
    ow = floor(w*SegmentOverlap);
    nF = length(lower_bin:upper_bin);
    for i = 1:N              
        [S,fvec,tvec]=spectrogram([F(i,end:-1:1) F(i,:) F(i,end:-1:1)],w,ow,lower_bin:upper_bin,Fs);
        ind = (tvec > length(F)/Fs) & (tvec <= 2*length(F)/Fs);
        S = S(:,ind);

        if isempty(TF)
            TF = zeros(sum(ind),length(fvec),N);
        end
        TF(:,:,i) = S';
    end
end
function TF = amp_vec(F,BandBounds,Fs,width,TFmethod) 
    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    N = size(F,1);
    T = size(F,2);
    TF = [];
    cf = (lower_bin + floor((upper_bin- lower_bin)/2));
    w = floor(width*Fs*(1/lower_bin));
    for i = 1:N
        if strcmp(TFmethod,'wavelet')
            F1 = ampvec(cf, F(i,:)', Fs, width)';          
        elseif strcmp(TFmethod,'hilbert')  
            % Band-pass filter in one frequency band
            if lower_bin ~= 0 && upper_bin ~= 0
                Fband = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1);
            else
                Fband = F(i,:);
            end
%             Fband = F(i,:);
            %Fband = Fband(T+(1:T));
            % Apply Hilbert transform
            F1 = abs(hilbert(Fband')');
%             F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
        elseif strcmp(TFmethod,'stft')            
            [S,fvec,tvec]=spectrogram([F(i,end:-1:1) F(i,:) F(i,end:-1:1)],w,w-1,lower_bin:upper_bin,Fs);
            ind = (tvec > length(F)/Fs) & (tvec <= 2*length(F)/Fs);
            F1 = mean(abs(S(:,ind)));
        end
        if isempty(TF)
            TF = zeros(N,length(F1));
        end
        TF(i,:) = F1;
    end
end

function TF = bp_vec(F,BandBounds,Fs,width,TFmethod)
%     N = size(F,1);
%     T = size(F,2);
%     TF = [];
    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    Fcell = num2cell(F,2);
%     if N>1
    if lower_bin ==0 && upper_bin==0
       TF = F;
    else
       TFcell = cellfun(@(x)process_bandpass('Compute', x, Fs, lower_bin, upper_bin, [], 1), Fcell,'UniformOutput',false);
       TF = cell2mat(TFcell);
    end
%     else
%         for i = 1:N
%             if strcmp(TFmethod,'wavelet')
%                 F1 = bpvec((lower_bin + floor((upper_bin- lower_bin)/2)),F(i,:)', Fs, width);
% 
%             elseif strcmp(TFmethod,'hilbert') || strcmp(TFmethod,'stft') 
% 
%                 F1 = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1)';
%     %             if strcmp(TFmethod,'hilbert')
%     %                 F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
%     %             end
%             end
%             if isempty(TF)
%                 TF = zeros(N,length(F1));
%             end
%             TF(i,:) = F1';
% 
%         end
%     end
    
end
% function TF = bp_vec(F,BandBounds,Fs,width,TFmethod)
%     N = size(F,1);
%     T = size(F,2);
%     TF = [];
%     lower_bin = BandBounds(1);
%     upper_bin = BandBounds(2);
%         
%     for i = 1:N
%         if strcmp(TFmethod,'wavelet')
%             F1 = bpvec((lower_bin + floor((upper_bin- lower_bin)/2)),F(i,:)', Fs, width);
% 
%         elseif strcmp(TFmethod,'hilbert') || strcmp(TFmethod,'stft') 
%             F1 = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1)';
% %             if strcmp(TFmethod,'hilbert')
% %                 F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
% %             end
%         end
%         if isempty(TF)
%             TF = zeros(N,length(F1));
%         end
%         TF(i,:) = F1';
% 
%     end
%     
% end
function TF = ph_vec(F,BandBounds,Fs,width,TFmethod) 

    lower_bin = BandBounds(1);
    upper_bin = BandBounds(2);
    N = size(F,1);
    T = size(F,2);
    TF = [];
    cf = (lower_bin + floor((upper_bin- lower_bin)/2));
    w = floor(width*Fs*(1/lower_bin));
    for i = 1:N
        if strcmp(TFmethod,'wavelet')
            F1 = phasevec(cf, F(i,:)', Fs, width)';          
        elseif strcmp(TFmethod,'hilbert')  
            % Band-pass filter in one frequency band
            Fband = process_bandpass('Compute', F(i,:), Fs, lower_bin, upper_bin, [], 1);
            % Apply Hilbert transform
            F1 = angle(hilbert(Fband')');
%             F1 = F1(ceil(0.1*T)+(1:floor(0.8*T))-1);
        elseif strcmp(TFmethod,'stft')            
            [S,fvec,tvec]=spectrogram([F(i,end:-1:1) F(i,:) F(i,end:-1:1)],w,w-1,lower_bin:upper_bin,Fs);
            ind = (tvec > length(F)/Fs) & (tvec <= 2*length(F)/Fs);
            F1 = mean(angle(S(:,ind)));
        end
        if isempty(TF)
            TF = zeros(N,length(F1));
        end
        TF(i,:) = F1;
    end


end


