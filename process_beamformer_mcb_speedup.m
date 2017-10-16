function varargout = process_beamformer_mcb_speedup( varargin )
% PROCESS_BEAMFORMER_MCB (2017.01.16): 
% USAGE:  sInput = process_beamformer_mcb_speedup('GetDescription')                      
%         sOutput = process_beamformer_mcb_speedup('Run', sProcess, sInput, method=[])
% INPUT: 
%     - Options
%           |- result_comm :    The comment of output files
%           |- oriconstraint:   If 1, unconstrained beamformer, if 0 uses normal to cortex as dipole orientation
%           |- minvar_range :   [tStart, tStop]; range of time values used to compute minimum variance criteria 
%           |- fmaps_range :    [tStart, tStop]; range of time values used to compute f maps 
%           |- fmap_size:       length of sliding window for computing fmap
%           |- fmaps_tresolution : Step between each fmap
%           |- reg :            regularization parameter
%           |- sensortypes:     MEG, EEG, MEG GRAD, MEG MAG
%           |- savefilter    :  If 1, saves the result spatial filter, if 0 only saves f maps
% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2014 University of Southern California & McGill University
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
% Authors: Hui-Ling Chan, Francois Tadel, Sylvain Baillet, 2014-2015
eval(macro_method);
%macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Beamformer: Maximum Contrast (speedup)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 330;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % Definition of the options
    sProcess.options.result_comm.Comment = 'Comment: ';
    sProcess.options.result_comm.Type = 'text';
    sProcess.options.result_comm.Value = '';
       
    % === SENSOR TYPES
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG GRAD';
    % === CORTICAL CONSTRAINED ORIENTATION
    sProcess.options.label_o.Comment = '<HTML><BR><b><u>Options: Dipole orietation</u></b>';
    sProcess.options.label_o.Type    = 'label';
    sProcess.options.oriconstraint.Comment = {'Unconstrained  <FONT color="#777777"><I>(maximum contrast of active/noise) </I></FONT>', 'Constrained <FONT color="#777777"><I>(normal to cortex)</I></FONT>'};
    sProcess.options.oriconstraint.Type    = 'radio';
    sProcess.options.oriconstraint.Value   = 1;  
    sProcess.options.ori_range.Comment = 'Time range (only for unconstr): ';
    sProcess.options.ori_range.Type    = 'poststim';
    sProcess.options.ori_range.Value   = [];
%     % === BASELINE TIME RANGE 
%     sProcess.options.baseline_time.Comment = 'Baseline (Set [0,0] if no baseline correction): ';
%     sProcess.options.baseline_time.Type    = 'baseline';
%     sProcess.options.baseline_time.Value   = [];
    
    sProcess.options.label_sf.Comment = '<HTML><BR><b><u>Options: Spatial filter</u></b>';
    sProcess.options.label_sf.Type    = 'label';
    % === SPATIAL FILTER TIME RANGE
    sProcess.options.minvar_range.Comment = 'Time range (for minimum variance): ';
    sProcess.options.minvar_range.Type    = 'poststim';
    sProcess.options.minvar_range.Value   = [];
    % === REGULARIZATION
    sProcess.options.reg.Comment = 'Regularization parameter: ';
    sProcess.options.reg.Type    = 'value';
    sProcess.options.reg.Value   = {0.3, '%', 4};
    % === SAVE SPATIAL FILTER
    sProcess.options.savefilter.Comment = 'Save spatial filters';
    sProcess.options.savefilter.Type    = 'checkbox';
    sProcess.options.savefilter.Value   = 1;  
    
    sProcess.options.label_f.Comment = '<HTML><BR><b><u>Options: Active state (for f-statistic map)</u></b>';
    sProcess.options.label_f.Type    = 'label';
    % === F-MAP TIME RANGE
    sProcess.options.fmaps_range.Comment = 'Duration of interests (DOI): ';
    sProcess.options.fmaps_range.Type    = 'poststim';
    sProcess.options.fmaps_range.Value   = [];
    % === F-MAP TIME WINDOW SIZE
    sProcess.options.fmap_size.Comment = 'Size of sliding window (0=whole DOI): ';
    sProcess.options.fmap_size.Type    = 'value';
    sProcess.options.fmap_size.Value   = {0, 'ms', 1};
    % === F-MAP TEMPORAL RESOLUTION
    sProcess.options.fmaps_tresolution.Comment = 'Step for sliding window: ';
    sProcess.options.fmaps_tresolution.Type    = 'value';
    sProcess.options.fmaps_tresolution.Value   = {0, 'ms', 1};

    sProcess.options.label_s.Comment = '<HTML><BR><b><u>Options: Scouts</u></b>';
    sProcess.options.label_s.Type    = 'label';
    % === OPTIONS: COMPUTE IN SELECT SCOUTS
    sProcess.options.usescouts.Comment = 'Compute sources in scouts (uncheck=whole brain)';
    sProcess.options.usescouts.Type    = 'checkbox';
    sProcess.options.usescouts.Value   = 0;   
    
    % === OPTIONS: SELECTED SCOUTS
    sProcess.options.scouts.Comment = 'Use scouts (no selection=all):';
    sProcess.options.scouts.Type = 'scout';
    sProcess.options.scouts.Value      = [];
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned list of files
    OutputFiles = {};
    % Selected scouts
    sScouts = sProcess.options.scouts.Value;
    % Get option values
    BaselineTime  = [0,0];%sProcess.options.baseline_time.Value{1};
    FmapRange   = sProcess.options.fmaps_range.Value{1};
    Reg         = sProcess.options.reg.Value{1};
    SensorTypes = sProcess.options.sensortypes.Value;
    
    FmapSize    = sProcess.options.fmap_size.Value{1};
    FmapTResolu = sProcess.options.fmaps_tresolution.Value{1};
    isSaveFilter = sProcess.options.savefilter.Value;
    %SFcriteriaType = sProcess.options.sfcriteria.Value;
    MinVarRange = sProcess.options.minvar_range.Value{1};
    OriRange = sProcess.options.ori_range.Value{1};
    isScout = sProcess.options.usescouts.Value;
    result_comment = sProcess.options.result_comm.Value;
    if ~isempty(result_comment)
        result_comment = [result_comment];
    end
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
    DataMat = in_bst(sInputs(1).FileName, [], 0);
    nChannels = size(DataMat.F,1);
    nTime     = size(DataMat.F,2);
    Time = DataMat.Time;
    
    % ===== PROCESS THE TIME WINDOWS =====
    if FmapRange(1) > FmapRange(2)
        bst_report('Error', sProcess, sInputs, 'The setting of time range of active state is incorrect.');
    end
    
    if FmapRange(1) < Time(1)
        bst_report('Warning', sProcess, sInputs, 'The start for time range of active state is reset to the first time point of data');
        FmapRange(1) = Time(1);
    end
    
    if FmapRange(2) > Time(end)
        bst_report('Warning', sProcess, sInputs, 'The end for time range of active state is reset to the end point of data');
        FmapRange(2) = Time(end);
    end
    
    if (FmapRange(1)+FmapSize) > FmapRange(2)
        bst_report('Warning', sProcess, sInputs, 'The F statistic window size is too large and reset to the same as the time range of active state.');
        FmapSize = FmapRange(2) - FmapRange(1);
    end
    MinVarTime = MinVarRange;
    ActiveTime = OriRange;
    %MinVarPoints = panel_time('GetTimeIndices', Time, MinVarRange);

    if FmapTResolu == 0
        nFmaps = 1;
        FmapTResolu = 1;
        
        FmapSize = FmapRange(2) - FmapRange(1);
        if FmapRange(1)+FmapSize ~= FmapRange(2)
            bst_report('Warning', sProcess, sInputs, 'Temporal resolution should not be 0 ms. The F statistic window size is reset to the same as the time range of interest.'); 
        end
    end
    
    FmapPoints = panel_time('GetTimeIndices', Time, [FmapRange(1)+FmapSize FmapRange(2)]);
    if length(FmapPoints)<=1
        nFmaps = 1;
    else 
        nFmaps = length((FmapRange(1)+FmapSize):FmapTResolu:FmapRange(2));
    end
       
    
    HalfFmapSize= FmapSize/2;
    FmapTimeList= zeros(nFmaps,2);
    for i=1:nFmaps
        FmapTimeList(i,:) = FmapRange(1) + FmapTResolu*(i-1) + [0 FmapSize] ;
    end     
    

    
    % ===== LOAD CHANNEL FILE =====
    % Load channel file
    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    % Find the MEG channels
%     iMEG = good_channel(ChannelMat.Channel, [], 'MEG');
%     iEEG = good_channel(ChannelMat.Channel, [], 'EEG');
%     iSEEG = good_channel(ChannelMat.Channel, [], 'SEEG');
%     iECOG = good_channel(ChannelMat.Channel, [], 'ECOG');
    iChannels = channel_find(ChannelMat.Channel, SensorTypes);
    
    % ===== LOAD HEAD MODEL =====
    % Get channel study
    [sChannelStudy, iChannelStudy] = bst_get('ChannelFile', sInputs(1).ChannelFile);
    % Load the default head model
    HeadModelFile = sChannelStudy.HeadModel(sChannelStudy.iHeadModel).FileName;
    sHeadModel = load(file_fullpath(HeadModelFile));
    % Get number of sources
    nSources = size(sHeadModel.GridLoc,1);
  
    if (sProcess.options.oriconstraint.Value == 2) && (isequal(sHeadModel.HeadModelType,'volume') || isempty(sHeadModel.GridOrient))
       bst_report('Error', sProcess, sInputs, 'No dipole orientation for cortical constrained beamformer estimation.');
       % Stop the process
       return;  
    end
    
    % ===== LOAD SCOUTS =====
    if isempty(sScouts) || (isScout==0) %isequal(sHeadModel.HeadModelType,'volume') || 
        isScout = 0;
        nScoutVertex = nSources;
        sScoutVerticesList = 1:nSources;
    else
        isScout = 1;
    end
    
    if isScout
        sScoutsInfo = process_extract_scout('GetScoutsInfo', '@ Beamformer:MCB', [], sHeadModel.SurfaceFile, sScouts);
        sScoutVerticesList = unique([sScoutsInfo.Vertices]);
        nScoutVertex = length(sScoutVerticesList);
        
        if isequal(sHeadModel.HeadModelType,'volume')
%             sSurface = load(file_fullpath(sHeadModel.SurfaceFile));
%             vSurf = sSurface.Vertices(sScoutVerticesList,:);
            newScoutVerticesList = [];
            for ns = 1:nScoutVertex
                SCS = sHeadModel.GridLoc(sScoutVerticesList(ns),:);
                dist = (sHeadModel.GridLoc(:,1) - SCS(1)) .^ 2 + (sHeadModel.GridLoc(:,2) - SCS(2)) .^ 2 + (sHeadModel.GridLoc(:,3) - SCS(3)) .^ 2;
                iVertex = find(sqrt(dist) < 0.01); 
%                 [~, iVertex] = min(dist);
                newScoutVerticesList = [newScoutVerticesList iVertex];
%                 sScoutVerticesList(ns) = iVertex;
            end
            sScoutVerticesList = unique(newScoutVerticesList);
            nScoutVertex = length(sScoutVerticesList);
            
%             iAtlas = find(arrayfun(@(x)strcmp(x.Name(1:6),'Volume'),sSurface.Atlas));
%             for ns = 1:length(sScoutsInfo)
%                 iScouts(ns) = find(arrayfun(@(x)strcmp(x.Label, sScoutsInfo(ns).Label),sSurface.Atlas(iAtlas).Scouts));
%             end
%             GridAtlas = sSurface.Atlas(iAtlas);
%             GridAtlas.Scouts = GridAtlas.Scouts(unique(iScouts));
%             
%             GridAtlas.Vert2Grid = cellfun(@(x)min_dist_index(x,sHeadModel.GridLoc),mat2cell(sSurface.Vertices,ones(size(sSurface.Vertices,1),1), 3));
%             GridAtlas.Grid2Source = 1:nSources;
%             for ns = 1:length(GridAtlas.Scouts)
%                 GridAtlas.Scouts(ns).Region(2) ='V';
%                 GridAtlas.Scouts(ns).Region(3) ='C'; 
%                 GridAtlas.Scouts(ns).GridRows = GridAtlas.Vert2Grid(GridAtlas.Scouts(ns).Seed);
%             end
            
%             % Remove the vertices that are outside the list of vertices in Vert2Grid
%             iVertices(iVertices > size(GridAtlas.Vert2Grid,2)) = [];
%             % Surface.Vertices => Results.GridLoc
%             iVertices = find(any(GridAtlas.Vert2Grid(:,iVertices), 2))';
%         
%             % Get indices in the ImageGridAmp/ImagingKernel matrix
%             iSourceRows = find(any(GridAtlas.Grid2Source(:,iVertices), 2))';
%             % Find over which regions this vertex selection spans
%             if (nargout >= 2)
%                 iRegionScouts = find(~cellfun(@(c)isempty(intersect(c,iVertices)), {GridAtlas.Scouts.GridRows}));
%             end
            
        end
    end
    
    
    
    % ===== LOAD THE DATA =====
    % Find the indices for covariance calculation
    iMinVarTime = panel_time('GetTimeIndices', Time, MinVarRange);
    iActiveTime = panel_time('GetTimeIndices', Time, OriRange);
    if BaselineTime(1)==BaselineTime(2)
        iBaselineTime = [];
    else
        iBaselineTime = panel_time('GetTimeIndices', Time, BaselineTime);
    end
    % Initialize the covariance matrices
    ActiveCov = zeros(nChannels, nChannels);
    nTotalActive = zeros(nChannels, nChannels);
    MinVarCov = zeros(nChannels, nChannels);
    nTotalMinVar = zeros(nChannels, nChannels);
    
    nFmapPoints = length(panel_time('GetTimeIndices', Time, FmapTimeList(1,:)));
    iFmapTime = zeros(nFmaps, nFmapPoints);
    for i = 1:nFmaps
        if FmapTimeList(i,1) < Time(1) || FmapTimeList(i,2) > Time(end)
            % Add an error message to the report
            bst_report('Error', sProcess, sInputs, 'One fmap time window is not within the data time range.');
            % Stop the process
            return;
        end           
        
        single_iFmap = panel_time('GetTimeIndices', Time, FmapTimeList(i,:));
        % detect the round error caused by the non-interger sampling rate
        % -- added in 20141026
        if  nFmapPoints - length(single_iFmap) == 1
            single_iFmap(nFmapPoints) = single_iFmap(end) + 1;
        end
        %%%%%%%%
        
        iFmapTime(i,:) = single_iFmap(1:nFmapPoints);    
    end
    % Initialize the covariance matrices
    FmapActiveCov = zeros(nChannels, nChannels, nFmaps);
    nTotalFmapActive = zeros(nChannels, nChannels, nFmaps); 

    bst_progress('start', 'Applying process: MCB', 'Calulating covariance matrix...', 0, length(sInputs)*nFmaps*4+2);
    
    AllData = [];
    
    % Reading all the input files in a big matrix
    for i = 1:length(sInputs)
        % Read the file #i
        DataMat = in_bst(sInputs(i).FileName, [], 0);
        % Check the dimensions of the recordings matrix in this file
        if (size(DataMat.F,1) ~= nChannels) || (size(DataMat.F,2) ~= nTime)
            % Add an error message to the report
            bst_report('Error', sProcess, sInputs, 'One file has a different number of channels or a different number of time samples.');
            % Stop the process
            return;
        end

        % Get good channels
        iGoodChan = find(DataMat.ChannelFlag == 1);
    

        if ~isempty(iBaselineTime)           
            % Average baseline values
            FavgActive = mean(DataMat.F(iGoodChan,iBaselineTime), 2);
            % Remove average
            DataActive = bst_bsxfun(@minus, DataMat.F(iGoodChan,iActiveTime), FavgActive);
        else
            DataActive = DataMat.F(iGoodChan,iActiveTime);
        end

        % Compute covariance for this file
        fileActiveCov = DataMat.nAvg .* (DataActive * DataActive');        
        % Add file covariance to accumulator
        ActiveCov(iGoodChan,iGoodChan) = ActiveCov(iGoodChan,iGoodChan) + fileActiveCov;       
        nTotalActive(iGoodChan,iGoodChan) = nTotalActive(iGoodChan,iGoodChan) + length(iActiveTime);

        if ~isempty(iBaselineTime)           
            % Remove average
            DataMinVar = bst_bsxfun(@minus, DataMat.F(iGoodChan,iMinVarTime), FavgActive);
        else
            DataMinVar = DataMat.F(iGoodChan,iMinVarTime);
        end

        % Compute covariance for this file
        fileMinVarCov = DataMat.nAvg .* (DataMinVar * DataMinVar');        
        % Add file covariance to accumulator
        MinVarCov(iGoodChan,iGoodChan) = MinVarCov(iGoodChan,iGoodChan) + fileMinVarCov;       
        nTotalMinVar(iGoodChan,iGoodChan) = nTotalMinVar(iGoodChan,iGoodChan) + length(iMinVarTime);
        
        
        for j = 1:nFmaps
            % Average baseline values
            if isempty(iBaselineTime) 
                FavgFmapActive = 0;           
            else
                FavgFmapActive = mean(DataMat.F(iGoodChan,iBaselineTime), 2); 
            end
                       
            % Remove average
            DataFmapActive = bst_bsxfun(@minus, DataMat.F(iGoodChan,iFmapTime(j,:)), FavgFmapActive);    
            %DataFmapActive = sym(DataFmapActive);
            
            %DataMatnAvgSym = sym(DataMat.nAvg);
            bst_progress('inc',1);
            % Compute covariance for this file
            fileFmapActiveCov =  DataMat.nAvg .* (DataFmapActive * DataFmapActive');   
            bst_progress('inc',1);
            % Add file covariance to accumulator
            FmapActiveCov(iGoodChan,iGoodChan,j) = FmapActiveCov(iGoodChan,iGoodChan,j) + fileFmapActiveCov; 
            %FmapActiveCov(iGoodChan,iGoodChan,j) = FmapActiveCov(iGoodChan,iGoodChan,j) + fileFmapActiveCov; 
            bst_progress('inc',1);
            nTotalFmapActive(iGoodChan,iGoodChan,j) = nTotalFmapActive(iGoodChan,iGoodChan,j) + nFmapPoints;
            bst_progress('inc',1);
        end
        
        if isempty(AllData)
            AllData = zeros(nChannels, length(DataMat.Time), length(sInputs));
        end
        AllData(:,:,i)=DataMat.F;
        
    
    end
    clear DataMat;
    % Remove bad channels - Added by Hui-Ling 20170113 - start
    iUnusedChannel = find(nTotalActive(1,:) < 1);
    if ~isempty(iUnusedChannel)
        for i=1:length(iUnusedChannel)
            iChannels(iChannels == iUnusedChannel(i)) = [];
        end
    end
    % Remove bad channels - Added by Hui-Ling 20170113 - end
    
    % Remove zeros from N matrix
    if sProcess.options.oriconstraint.Value == 1,
        nTotalActive(nTotalActive <= 1) = 2;
        nTotalMinVar(nTotalMinVar <= 1) = 2;
    end
    nTotalFmapActive(nTotalFmapActive <= 1) = 2;
    bst_progress('inc',1);
    % Divide final matrix by number of samples
    if sProcess.options.oriconstraint.Value == 1,
        ActiveCov = ActiveCov ./ (nTotalActive - 1);
        MinVarCov = MinVarCov ./ (nTotalMinVar - 1);
    end
    FmapActiveCov = FmapActiveCov ./ (nTotalFmapActive - 1);
    bst_progress('inc',1);
    bst_progress('stop');
    % ===== PROCESS =====

    % Number of channels used to compute sources
    nUsedChannels = length(iChannels); 
    
    % Extract the covariance matrix of the used channels
    ActiveCov = ActiveCov(iChannels,iChannels);
    NoiseCovMat = load(file_fullpath(sChannelStudy.NoiseCov(1).FileName));
    if isempty(NoiseCovMat)
        bst_report('Error', sProcess, sInputs, 'Cannot find noise covariance matrix.');
        % Stop the process
        return;
    end
    NoiseCov  = NoiseCovMat.NoiseCov(iChannels,iChannels);
    MinVarCov = MinVarCov(iChannels,iChannels);
    FmapActiveCov = FmapActiveCov(iChannels,iChannels,:);
    AllData = AllData(iChannels,:,:);
 
    %% Calculate the inverse of (C+alpha*I)
%      eigValues = eig(MinVarCov);
%      Reg_alpha = Reg / 100 * max(eigValues);
%      invMinVarCovI = inv(MinVarCov+Reg_alpha*eye(nUsedChannels));
 
     [U,S,V] = svd(MinVarCov);
     S = diag(S); % Covariance = Cm = U*S*U'
     Si = diag(1 ./ (S + S(1) * Reg / 100)); % 1/(S + lambda I)
     invMinVarCovI = V*Si*U'; % V * 1/(S + lambda I) * U' = Cm^(-1) 
    
    % Initilize the ImagingKernel (spatial filter) and ImageGridAmp (f value)
    ImagingKernel = zeros(nSources, nUsedChannels);
    ImageGridAmp = zeros(nSources, nFmaps);
    ImagingGridOri = zeros(nSources, 3);
    
    % Set the dipole positions for the computation of sources
    Loc = sHeadModel.GridLoc;
    
    
    % Get forward field
    if sProcess.options.oriconstraint.Value == 2; % anatomically constrained beamformer
        Kernel = bst_gain_orient(sHeadModel.Gain(iChannels,:), sHeadModel.GridOrient);
%         Loc = sHeadModel.GridLoc;
    else
        Kernel = sHeadModel.Gain(iChannels,:);
%         Loc = sHeadModel.GridLoc;
    end
    %Kernel(abs(Kernel(:)) < eps) = eps; % Set zero elements to strictly non-zero
    bst_progress('start', 'Applying process: MCB', 'Calculating spatial filters and f-statistic maps...', 0, 16+nFmaps);
    
    %% Compute the spatial filter and f value for each position   
    if sProcess.options.oriconstraint.Value == 1; 
        % Obtain gain matrix
        usedKernel = arrayfun(@(x)Kernel(:,[1 2 3]+3*(x-1)),sScoutVerticesList,'UniformOutput', false);
        % Compute A = inv(C+alpha*I)*Lr
        Amat = cellfun(@(x)invMinVarCovI*x, usedKernel, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % == Compute orientation using maxium constrast criterion ==
        % Compute P = A'*Ca*A
        Pmat = cellfun(@(x) x'*ActiveCov*x, Amat, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % Compute Q = A'*Cc*A
        Qmat = cellfun(@(x) x'*NoiseCov*x, Amat, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % Regularize the matrix Q to avoid singular problem 
        invQmat = cellfun(@(x)pInv(x,0.0000000001),Qmat, 'UniformOutput', false);
        bst_progress('inc',1);
        invQPmat = cellfun(@mtimes, invQmat, Pmat, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % Compute the dipole orientation 
        % (the eigenvector corresponding to maximum eigenvalue of inv(Q)*P)
        [eigVecmat, eigValmat] = cellfun(@eig, invQPmat, 'UniformOutput', false);
        bst_progress('inc',1);
        % check whether eigValues are saved as matrix or vector
        if min(size(eigValmat{1,1}))>1
            eigValmat = cellfun(@diag, eigValmat, 'UniformOutput', false);
        end
        [~, imaxmat] = cellfun(@max, eigValmat, 'UniformOutput', false);
        bst_progress('inc',1);
        Orimat = cellfun(@(x,y)x(:,y), eigVecmat, imaxmat, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % Compute B = Lr'*A
        Bmat = cellfun(@(x,y)x'*y, usedKernel, Amat, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % Compute the spatial filter
        tmpMat = cellfun(@(x,y)x'*y*x, Orimat, Bmat, 'UniformOutput', false);
        bst_progress('inc',1);
        tmpMat2 = cellfun(@mtimes, Amat, Orimat, 'UniformOutput', false);
        bst_progress('inc',1);
        sfmat = cellfun(@mrdivide, tmpMat2, tmpMat, 'UniformOutput', false);
        bst_progress('inc',1);
        
        % Compute the source power during control state
        vcMat = cellfun(@(x)x'*NoiseCov*x, sfmat,'UniformOutput', false);
        bst_progress('inc',1);
        
        
        for j = 1:nFmaps
             % Compute the source power during active state
            vaMat = cellfun(@(x)x'*FmapActiveCov(:,:,j)*x, sfmat,'UniformOutput', false);
            
             % Compute the f-statistic value by contrasting the power
             % during active and control states
            ImageGridAmp(sScoutVerticesList,j) = cellfun(@rdivide,vaMat,vcMat);
            bst_progress('inc',1);
        end
        
        % Normalize spatial filter using the amplitude of control state
        vcMat = cellfun(@sqrt, vcMat, 'UniformOutput', false);
        bst_progress('inc',1);
        tmp = cellfun(@mrdivide, sfmat, vcMat, 'UniformOutput', false);
        bst_progress('inc',1);
        if size(tmp{1},2) == 1
            tmp = cellfun(@(x)x', tmp, 'UniformOutput', false);
        end
        if size(tmp,1) == 1
            tmp = tmp';
        end
        % Save the normalized spatial filter
        ImagingKernel(sScoutVerticesList,:) = cell2mat(tmp);
        ImagingGridOri(sScoutVerticesList,:) = cell2mat(Orimat)';
        bst_progress('inc',1);
%         for iScoutVertex = 1:nScoutVertex
%                 
%             i = sScoutVerticesList(iScoutVertex);
%             iGain = [1 2 3] + 3*(i-1);
%             if Loc(i,1) < 0.01 && Loc(i,1) > -0.01 && Loc(i,2) < 0.01 && Loc(i,2) > -0.01 && Loc(i,3) < 0.01 && Loc(i,3) > -0.01
%                 bst_progress('inc',2);
%                 
%                 continue;
%             else
%                 % Compute A = inv(C+alpha*I)*Lr
%                 A = invMinVarCovI * Kernel(:,iGain);
% 
%                 % Compute B = Lr'*A
%                 B = Kernel(:,iGain)' * A;
% 
% 
%                 % == Compute orientation using maxium constrast criterion ==
%                 % Compute P = A'*Ca*A
%                 P = A' * ActiveCov * A;
%                 % Compute Q = A'*Cc*A
%                 Q = A' * NoiseCov * A;
% 
%                 % Regularize the matrix Q to avoid singular problem 
%                 [U,S,V] = svd(Q);
%                 S = diag(S); % Covariance = Cm = V*S*U'
%                 Si = diag(1 ./ (S + S(1) * 0.0000000001)); % 1/(S + lambda I)
%                 invQ = V*Si*U'; % V * 1/(S + lambda I) * U' = Cm^(-1) 
% 
%                 % Compute the dipole orientation 
%                 % (the eigenvector corresponding to maximum eigenvalue of inv(Q)*P)
%                 [eigVectors,eigValues] = eig(invQ*P);
%                 % check whether eigValues are saved as matrix or vector
%                 if(min(size(eigValues))==1)
%                     [tmp, imax] = max(eigValues);
%                 else
%                     [tmp, imax] = max(diag(eigValues));
%                 end
%                 DipoleOri = eigVectors(:,imax);
% 
%                 % Compute the spatial filter
%                 SpatialFilter = (A * DipoleOri) / (DipoleOri' * B * DipoleOri);           
% 
%                 varControl = SpatialFilter'* NoiseCov * SpatialFilter;
%                 bst_progress('inc',1);
%                 
%                 
%                 
%                 for j = 1:nFmaps
%                     % Compute the contrast of source power during active state and control state
%                     varActive = SpatialFilter'*FmapActiveCov(:,:,j)*SpatialFilter;
%     %                 varActive = 0;
%     %                 for k = 1:length(sInputs)
%     %                     varActive = varActive + sum((SpatialFilter' * AllData(:,iFmapTime(j,:),k)).^2);
%     %                 end
%     %                 varActive = varActive / (length(sInputs)*nFmapPoints);
%                     %varActive = mean((SpatialFilter' * AllData).^2); 
%                     ImageGridAmp(i,j) = varActive / varControl;
% 
%                     bst_progress('inc',1);
%                 end
%                 % Save teh result 
%                 ImagingKernel(i,:) = SpatialFilter / sqrt(varControl);
%                 
%             end
%             
%         end

    else
        for iScoutVertex = 1:nScoutVertex
            i = sScoutVerticesList(iScoutVertex);
            %iGain = [1 2 3] + 3*(i-1);
            if Loc(i,1) < 0.01 && Loc(i,1) > -0.01 && Loc(i,2) < 0.01 && Loc(i,2) > -0.01 && Loc(i,3) < 0.01 && Loc(i,3) > -0.01
                bst_progress('inc',2);
                
                continue;
            end
            % Compute the spatial filter with cortical constrained dipole orientation

            % Compute A = inv(C+alpha*I)*Lr
            A = invMinVarCovI * Kernel(:,i);

            % Compute B = Lr'*A
            B = Kernel(:,i)'*A;

            % Compute the spatial filter
            SpatialFilter = A / B;

            % compute the variance of control state            
            varControl = SpatialFilter'* NoiseCov * SpatialFilter;
            bst_progress('inc',1);

            for j = 1:nFmaps
                % Compute the contrast of source power during active state and control state
                varActive = SpatialFilter'*FmapActiveCov(:,:,j)*SpatialFilter;
                
                ImageGridAmp(i,j) = varActive / varControl;
                %ImageGridAmp(i,j) = Fvalue;
                bst_progress('inc',1);
            end
            % Save teh result 
            ImagingKernel(i,:) = SpatialFilter /  sqrt(varControl);

        end
    end

    bst_progress('stop');
    
   
    bst_progress('start', 'Applying process: MCB', 'Interpolating results...', 0, 1);
    %bst_progress('text', ['Applying process: ' sProcess.Comment ' [Interpolating results]']);
    
    if nFmaps == 1          
        
        %%%%%%%%% ADDED BY HUI-LING May 10, 2016 -- start
        [mv, mr] = max(ImageGridAmp);
        
%         sSubject = bst_get('Subject', 0);
        
        [~, iSubject] = bst_get('SurfaceFile', sHeadModel.SurfaceFile);
        sMRI = bst_memory('LoadMri', iSubject);
        
%         sMRI = load(file_fullpath(sSubject.Anatomy.FileName));
        lc = cs_convert(sMRI, 'scs', 'mni', Loc(mr,:))*1000;
        disp(['MCB> Maximum Peak Location (MNI coordinates): ' num2str(round(lc)) ]);
        disp(['MCB> Maximum Peak Value (f-statistic): ' num2str(mv) ]);
        %%%%%%%%% ADDED BY HUI-LING May 10, 2016 -- end
        
        FmapRangePoints = panel_time('GetTimeIndices', Time, FmapRange);
        ImageGridAmpOriginal = ImageGridAmp;
        ImageGridAmp = [ImageGridAmpOriginal, ImageGridAmpOriginal];
        TimeIndex = [Time(1), Time(end)];      
        
    else % Interpolate the f maps to have the same temporal resolution as the data

        ImageGridAmpOriginal = ImageGridAmp;
        ImageGridAmp = zeros(nSources,nTime);
        for i=1:(nFmaps-1)
            InterpolateTimeWindow = FmapRange(1)+HalfFmapSize+[(i-1)*FmapTResolu i*FmapTResolu];
            iInterpolateTime = panel_time('GetTimeIndices', Time, InterpolateTimeWindow);
            nInterpolateTime = length(iInterpolateTime);
            
            ImageGridAmp(:,iInterpolateTime(1)) = ImageGridAmpOriginal(:,i); 
            ImageGridAmp(:,iInterpolateTime(end)) = ImageGridAmpOriginal(:,i+1); 
            
            for j=2:(nInterpolateTime-1)       
                InterpolatePercentage = (j-1)/(nInterpolateTime-1);
                ImageGridAmp(:,iInterpolateTime(j)) = ImageGridAmpOriginal(:,i+1)*InterpolatePercentage + ImageGridAmpOriginal(:,i)*(1-InterpolatePercentage);
            end
        end        
        %%%%%%%%%%%%
        
%         InterpolateTimeWindow = FmapRange(1)+ [0 HalfFmapSize];
%         iInterpolateTime = panel_time('GetTimeIndices', Time, InterpolateTimeWindow);
%         nInterpolateTime = length(iInterpolateTime);
% 
%         for j=1:(nInterpolateTime-1)       
%             ImageGridAmp(:,iInterpolateTime(j)) = ImageGridAmpOriginal(:,1);
%         end
%         
%         InterpolateTimeWindow = FmapRange(1)+HalfFmapSize+(nFmaps-1)*FmapTResolu+[0 HalfFmapSize];
%         iInterpolateTime = panel_time('GetTimeIndices', Time, InterpolateTimeWindow);
%         nInterpolateTime = length(iInterpolateTime);
% 
%         for j=2:nInterpolateTime       
%             ImageGridAmp(:,iInterpolateTime(j)) = ImageGridAmpOriginal(:,end);
%         end
%         
        [mv, mr] = max(max(ImageGridAmp'));
        [~, mt] = max(max(ImageGridAmp));
        
        sSubject = bst_get('Subject', 0);
        sMRI = load(file_fullpath(sSubject.Anatomy.FileName));
        lc = cs_convert(sMRI, 'scs', 'mni', Loc(mr,:)*1000);
        disp(['MCB> Maximum Peak Location (MNI coordinates): ' num2str(round(lc)) ]);
        disp(['MCB> Maximum Peak Value (f-statistic): ' num2str(mv) ]);
        disp(['MCB> Maximum Peak Time: ' num2str(Time(mt)) ' seconds']);
        
        TimeIndex = Time;
    end
    bst_progress('inc',1);
    bst_progress('stop');
    
    %bst_progress('text', ['Applying process: ' sProcess.Comment ' [Saving results]']);
    % ===== SAVE THE RESULTS =====
    % Create a new data file structure
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = ImageGridAmp;
%     ResultsMat.nComponents   = 1;   % 1 or 3
    if strcmp(sHeadModel.HeadModelType,'volume')  
        ResultsMat.nComponents   = 1;
    else
        ResultsMat.nComponents   = 1;   % 1 or 3
    end
    % Comment
    Comment = [];
    if ~isempty(result_comment)
        Comment = [result_comment ': '];
    end
    
    if ~isempty(iBaselineTime)
        strContrastType = 'bl';
    else
        strContrastType = 'no bl';
    end
    
    if sProcess.options.oriconstraint.Value == 1;
        ostr = 'Unconstr';
    else
        ostr = 'Constr';
    end
    
    if FmapRange(1) > 5 || FmapRange(2) > 5 || FmapSize > 5
        timescale = 1;
        strTimeUnit = 's';
    else
        timescale = 1000;
        strTimeUnit = 'ms';
    end
    
    if (nFmaps == 1)
        Comment1 = sprintf('%sMCB/fmap (%s, %s, %d-%d%s)', Comment, ostr, strContrastType, round(FmapRange(1)*timescale), round(FmapRange(2)*timescale),strTimeUnit);
    else
        Comment1 = sprintf('%sMCB/fmap (%s, %s , %d-%d%s, ws:%d%s, tr:%d%s)', Comment, ostr, strContrastType, round((FmapRange(1)+FmapSize/2)*timescale), ...
            round((FmapRange(1) + FmapSize/2 + (nFmaps-1)*FmapTResolu)*timescale), strTimeUnit, round(FmapSize*timescale), strTimeUnit, round(FmapTResolu*timescale),strTimeUnit);
    end
    
    ResultsMat.Function      = 'MaximumContrastBeamformerResult';
    ResultsMat.Comment       = Comment1;
    ResultsMat.Time          = TimeIndex;           % Leave it empty if using ImagingKernel
    ResultsMat.DataFile      = [];
    ResultsMat.HeadModelFile = HeadModelFile;
    ResultsMat.HeadModelType = sHeadModel.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = iChannels;
    ResultsMat.SurfaceFile   = sHeadModel.SurfaceFile;
    if strcmp(sHeadModel.HeadModelType,'volume')     
        ResultsMat.GridLoc       = Loc;
%         ResultsMat.GridAtlas     = GridAtlas;
    end

    % === NOT SHARED ===
    % Get the output study (pick the one from the first file)
    iStudy = sInputs(1).iStudy;
    % Create a default output filename 
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'results_MCB_amp');

    % Save on disk
    save(OutputFiles{1}, '-struct', 'ResultsMat');
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, ResultsMat);

        
    % ===== SPATIAL FILTER: SAVE FILE =====
    if isSaveFilter
        % == Save the spatial filter as ImagingKernel ==
        % Create a new data file structure
        ResultsMat2 = db_template('resultsmat');
        ResultsMat2.ImagingKernel = ImagingKernel;
        ResultsMat2.ImageGridAmp  = [];
        if strcmp(sHeadModel.HeadModelType,'volume')  
            ResultsMat2.nComponents   = 1;
        else
            ResultsMat2.nComponents   = 1;   % 1 or 3
        end

        timestring = sprintf('%d_%d%s',round(ActiveTime(1)*timescale),round(ActiveTime(2)*timescale),strTimeUnit);
        if ~isempty(iBaselineTime)
            ResultsMat2.Comment   = [Comment 'MCB/filter (' ostr ', bl, ' timestring ')'];
        else
            ResultsMat2.Comment   = [Comment 'MCB/filter (' ostr ', no bl, ' timestring ')'];
        end
        
        ResultsMat2.Function      = 'MaximumContrastBeamformerFilter';
        ResultsMat2.Time          = [];           % Leave it empty if using ImagingKernel
        ResultsMat2.DataFile      = [];
        ResultsMat2.HeadModelFile = HeadModelFile;
        ResultsMat2.HeadModelType = sHeadModel.HeadModelType;
        ResultsMat2.ChannelFlag   = [];
        ResultsMat2.GoodChannel   = iChannels;
        ResultsMat2.SurfaceFile   = sHeadModel.SurfaceFile;
        if strcmp(sHeadModel.HeadModelType,'volume')     
            ResultsMat2.GridLoc       = Loc;
%             ResultsMat2.GridAtlas     = GridAtlas;
        end
        if sProcess.options.oriconstraint.Value == 1
            ResultsMat2.EstimatedGridOrient = ImagingGridOri;
        end
        
        % === SHARED ==
        % Get the output study (pick the one from the first file)
        iStudy = iChannelStudy;
        % Create a default output filename 
        OutputFiles{2} = bst_process('GetNewFilename', fileparts(sInputs(1).ChannelFile), 'results_MCB_KERNEL');

        % Save on disk
        save(OutputFiles{2}, '-struct', 'ResultsMat2');
        % Register in database
        db_add_data(iStudy, OutputFiles{2}, ResultsMat2);
        %%===========
    end
%     
%     if (sProcess.options.oriconstraint.Value == 1) && (isempty(sHeadModel.GridOrient)==0)
%         % ===== SAVE THE RESULTS =====
%         % Create a new data file structure
%         ResultsMat3 = db_template('resultsmat');
%         ResultsMat3.ImagingKernel = [];
%         ResultsMat3.ImageGridAmp  = [OriDifference OriDifference];
%         ResultsMat3.nComponents   = 1;   % 1 or 3
%         ResultsMat3.Comment   = 'MCB: Orientation Difference(Unconstr)';
%         ResultsMat3.Function      = 'MaximumContrastBeamformerOriDiff';
%         ResultsMat3.Time          = [1 1];           % Leave it empty if using ImagingKernel
%         ResultsMat3.DataFile      = [];
%         ResultsMat3.HeadModelFile = HeadModelFile;
%         ResultsMat3.HeadModelType = sHeadModel.HeadModelType;
%         ResultsMat3.ChannelFlag   = [];
%         ResultsMat3.GoodChannel   = iChannels;
%         ResultsMat3.SurfaceFile   = sHeadModel.SurfaceFile;
%         ResultsMat3.GridLoc       = GridLoc;
% 
%         % === NOT SHARED ===
%         % Get the output study (pick the one from the first file)
%         iStudy = sInputs(1).iStudy;
%         % Create a default output filename 
%         OutputFiles{3} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'results_MCB_oriDiff');
%         % Save on disk
%         save(OutputFiles{3}, '-struct', 'ResultsMat3');
%         % Register in database
%         db_add_data(iStudy, OutputFiles{3}, ResultsMat3);        
%     end

end

function ind = min_dist_index(source,target)
    
    dist = (source(:,1) - target(1)) .^ 2 + (source(:,2) - target(2)) .^ 2 + (source(:,3) - target(3)) .^ 2;
    
    [~,ind] = min(dist);

end
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


