%% script to test various settings in fieldtrip used to generate STA spectra and average STS on a trial-by trial basis as well as for binned trials
% trials are binned on the basis of mean firing rate in the trial such that
% the number of spikes in each of the bins are as close as possible
% The spectra for the FSIs are calculated after multiple rounds of
% subsampling
%% setup
% Setup rng seed for reproducibility
rng(4994);

clear;
cd('E:\ADRLabData');
please = [];
please.rats = {'R117','R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.existing_res_dir = 'D:\vStr_oscillatory_switch_results\temp3';
cfg_in.output_dir = 'D:\vStr_oscillatory_switch_results\temp4';
cfg_in.incl_types = [1, 2];
cfg_in.nMinSpikes1 = 400; % For on track
cfg_in.nMinSpikes2 = 400; % For near and away
cfg_in.nMinSpikes3 = 150; % For lfr, hfr, p1 and p2
cfg_in.nControlSplits = 100;
cfg_in.num_subsamples = 1000;

% SPECIFIC CELL ANALYSIS PARAMETERS
% Option 1: Single cell (uncomment and modify the line below)
% cell_name = 'R117-2007-06-01-TT04_2.t'; % Change this to your target cell name

% Option 2: Multiple cells from text file (uncomment and modify the line below)
cell_list_file = 'D:\vStr_oscillatory_switch_results\temp4\fsi_list.txt'; % Text file with one cell name per line

% PPC METHOD TO TEST
cfg_in.ppc_method = 'ppc0'; % Try different methods: 'ppc0', 'ppc1', 'ppc2', 'plv', 'ral', 'ang', 'angin', 'angout'

% ANALYSIS PARAMETERS
cfg_in.min_freq = 2; % Minimum frequency in Hertz
cfg_in.pl_thresh = 99; % Percentile threshold to establish significance of phase locking
cfg_in.diff_thresh = 95; % Percentile threshold to establish significance of PPC difference
cfg_in.peak_freq_tol_win = 5; % Window around peak frequency for diff analysis (Hz)

%%
% Load cell list
if exist('cell_list_file', 'var') && exist(cell_list_file, 'file')
    % Read cells from text file
    fid = fopen(cell_list_file, 'r');
    cell_names = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line) && ~startsWith(line, '%')
            cell_names{end+1} = line;
        end
    end
    fclose(fid);
    fprintf('Loaded %d cells from %s\n', length(cell_names), cell_list_file);
elseif exist('cell_name', 'var')
    % Single cell mode
    cell_names = {cell_name};
    fprintf('Single cell mode: %s\n', cell_name);
else
    error('No cell specified. Please set either cell_name or cell_list_file.');
end

%%
% Top level loop which calls the main function for all the sessions
% Process each cell
for iCell = 1:length(cell_names)
    current_cell = cell_names{iCell};
    fprintf('\n=== Processing cell %d/%d: %s ===\n', iCell, length(cell_names), current_cell);
    
    % Extract session from cell name
    session_from_cell = extractBefore(current_cell, '-TT');
    fprintf('Looking for session: %s\n', session_from_cell);
    
    for iS = 1:length(cfg_in.fd) % for each session...
        tokens = split(cfg_in.fd{iS},'\');
        if strcmp(tokens{end}, session_from_cell)
            cfg_in.iS = iS;
            cfg_in.cell_name = current_cell;
            pushdir(cfg_in.fd{iS});
            analyzeSpecificCell(cfg_in); % do the business
            popdir;    
            break
        end
    end % of sessions
end % of cells

%%
% Main function to analyze a specific cell with different methods
function analyzeSpecificCell(cfg_in)

    LoadExpKeys;
    
    if isfield(ExpKeys,'goodGamma_vStr')
        cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
    elseif isfield(ExpKeys, 'goodGamma')
        cfg = []; cfg.fc = ExpKeys.goodGamma;
    else
        error('Couldn''t find LFP field name.');
    end

    csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % could locdetrend to improve STA estimate
    % Note that Fieldtrip will interpoalate data to avoid gaps. Include
    % sanity tests to ensure STA/STS segments don't include these gaps.
    ft_csc = ft_read_neuralynx_interp(cfg.fc);
    
    % hacky code to separate out only ontrack-data
    temp_tvec = [0:length(ft_csc.time{1})-1];
    temp_offset = (double(ft_csc.hdr.LastTimeStamp)/1e6 - double(ft_csc.hdr.FirstTimeStamp)/1e6)/(length(temp_tvec) - 1);
    temp_tvec = temp_tvec * temp_offset;
     
    % Modify ft_csc
    ft_csc.time{1} = temp_tvec;
    ft_csc.fsample = 1/temp_offset;
    
    % Ensure that the new data doesn't have any Nans in it
    temp_tvec = temp_tvec + double(ft_csc.hdr.FirstTimeStamp)/1e6;
    temp_start = nearest_idx3(ExpKeys.TimeOnTrack, temp_tvec);
    temp_end = nearest_idx3(ExpKeys.TimeOffTrack, temp_tvec);
    cfg_onTrack.begsample = temp_start;
    cfg_onTrack.endsample = temp_end;
    data_onTrack = ft_redefinetrial(cfg_onTrack, ft_csc);
    
    if ~isempty(find(isnan(data_onTrack.time{1}),1))
        warning('On Track data for %s has gaps', cfg.fc{1});
    end
        
    lfp_tt = regexp(cfg.fc, 'CSC\d+', 'match');
    lfp_tt = str2double(lfp_tt{1}{1}(4:end)); % need this to skip cells from same tt (could make into function)
    fprintf('LFP ttno is %d\n', lfp_tt);
    
    % params
    cfg_master = []; % overall params
    cfg_master.dt = 0.001;
    cfg_master.ccMethod = 'MvdM'; % cell type classification method
    cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
    cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
    cfg_master.iS = 1; % current session number out of fd list, get this from input cfg
    cfg_master.fd = []; % full list of session fd's, get this from input cfg
    cfg_master.fd_extra = []; % get this from input cfg
    cfg_master.write_output = 0;
    cfg_master.output_prefix = 'sts_';
    cfg_master.output_dir = 'C:\temp';
    cfg_master.incl_cell_types = [1,2]; %cell types to be included
    cfg_master.num_subsamples = 1000;
    cfg_master = ProcessConfig(cfg_master,cfg_in);
    
    % spikes
    sd.S = LoadSpikesTarget(cfg_master);
    
    % Restrict spikes to only OnTrack
    sd.S = restrict(sd.S, iv(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack));
    
    % Categorize cells and add tetrode depths
    cfg_wv = []; cfg_wv.cMethod = cfg_master.ccMethod;
    s_out = CategorizeStriatumWave(cfg_wv, sd.S);

    s_out.unit = [s_out.other s_out.msn s_out.fsi];
    s_out.ident = [zeros(1, length(s_out.other)) ones(1, length(s_out.msn)) repmat(2, 1, length(s_out.fsi))];

    cfg_tt = []; cfg_tt.verbose = 1;
    cfg_tt.this_rat_ID = cfg_master.fd_extra.ratID_num(cfg_master.iS);
    cfg_tt.this_date = cfg_master.fd_extra.fd_date_num(cfg_master.iS);

    for iC = 1:length(sd.S.t)
        % Read the spike files into field trip format
        sd.S.ft_spikes(iC) = ft_read_spike(sd.S.label{iC});
        sd.S.usr.cell_type(iC) = s_out.ident((s_out.unit == iC));
        sd.S.usr.tetrodeDepths(iC) = ExpKeys.TetrodeDepths(sd.S.usr.tt_num(iC));
        cfg_tt.ttno = sd.S.usr.tt_num(iC);
        [sd.S.usr.distanceTurned(iC), prev_fd] = DistanceTurned(cfg_tt, cfg_master.fd, cfg_master.fd_extra);
        cfg_tt.verbose = 0;
    end
     
    % correlate with previous session waveforms if available
    if isempty(prev_fd) % no previous day available
        sd.S.usr.duplicate = zeros(size(sd.S.usr.tt_num));
    else
        pushdir(prev_fd);
        S2 = LoadSpikes([]);
        nSpikes = cellfun(@length, S2.t); keep = nSpikes >= cfg_master.nMinSpikes1;
        S2 = SelectTS([], S2, keep);

        s_out2 = CategorizeStriatumWave(cfg_wv, S2);
        s_out = CalcWVDistances([], s_out, s_out2); % add comparison with previous day's waveforms

        popdir; 

        % for each cell in current session, decide if duplicate
        for iC = 1:length(sd.S.t)

            this_tt_no = sd.S.usr.tt_num(iC);
            prev_day_cells = find(S2.usr.tt_num == this_tt_no);

            if isempty(prev_day_cells) % no cells recorded fron this tt in previous session
                sd.S.usr.duplicate(iC) = 0;
            else % previous day cells found
                temp_corr = s_out.corr(iC, prev_day_cells);
                temp_peakn = s_out.peakdiffn(iC, prev_day_cells);

                if temp_corr > cfg_master.maxPrevCorr & abs(temp_peakn) < cfg_master.maxPeakn % wv correlation big, peak difference small
                    sd.S.usr.duplicate(iC) = 1;
                else
                    sd.S.usr.duplicate(iC) = 0;
                end
            end
        end
    end % of previous day available checks
    
    % PLEASE SEE: Use the sd.S.usr.duplicate field if you are calculating
    % it anyway!!
    % Keep only non duplicate cells and those which were read by FieldTrip 
    % correctly
    % Also Get rid of cells of exc_types or spikes < nMinSpikes1;
    % Also get rid of cells on the LFP tt
    keep = false(1,length(sd.S.t));
    for iK = 1:length(keep)
        for iT = 1:length(cfg_master.incl_cell_types)
           keep(iK) = keep(iK) | (cfg_master.incl_cell_types(iT) == sd.S.usr.cell_type(iK)); 
        end
        keep(iK) = keep(iK) & (length(sd.S.t{iK}) > cfg_master.nMinSpikes1) & ~(sd.S.usr.tt_num(iK) == lfp_tt);
    end
    keep = keep & ~sd.S.usr.duplicate;% & sd.S.ft_spk_valid;
    sd.S = SelectTS([], sd.S, keep);
    sd.S.ft_spikes = sd.S.ft_spikes(keep);
    od.cell_type = sd.S.usr.cell_type;
    od.tt_id = sd.S.usr.tt_num;
    od.label = sd.S.label;
   
    % Find the specific cell we want to analyze
    target_cell_idx = find(strcmp(sd.S.label, cfg_in.cell_name));
    if isempty(target_cell_idx)
        error('Target cell %s not found in session', cfg_in.cell_name);
    end
    
    iC = target_cell_idx(1);
    fprintf('Analyzing cell %s (type %d)\n', sd.S.label{iC}, od.cell_type(iC));
    
    % Load existing results if available
    existing_results = loadExistingResults(cfg_in, sd.S.label{iC});
    
    % Calculate new measures for the specific cell
    new_results = calculateNewMeasures(cfg_in, cfg_master, sd, iC, ft_csc, cfg_onTrack, ExpKeys);
    
    % Plot comparison
    plotCellComparison(cfg_in, sd.S.label{iC}, existing_results, new_results);
    
    % Save results to temp file
    [~, fp, ~] = fileparts(pwd);
    results_file = fullfile(cfg_in.output_dir, sprintf('%s_%s_results.mat', fp, cfg_in.cell_name));
    save(results_file, 'new_results', 'existing_results', 'cfg_in');
    fprintf('Results saved to: %s\n', results_file);
end

%%
% Function to load existing results
function existing_results = loadExistingResults(cfg_in, cell_label)
    existing_results = [];
    
    % Try to load the existing results file
    try
        [~, fp, ~] = fileparts(pwd);
        results_file = fullfile(cfg_in.existing_res_dir, [fp '_ft_spec.mat']);
        if exist(results_file, 'file')
            temp_load = load(results_file);
            od = temp_load.od;
            
            % Find the cell in the results
            cell_idx = find(strcmp(od.label, cell_label));
            if ~isempty(cell_idx)
                if od.cell_type(cell_idx) == 1 % MSN
                    if isfield(od, 'msn_res')
                        % Find which MSN this cell corresponds to among all MSNs
                        msn_cells = find(od.cell_type == 1);
                        msn_idx = find(msn_cells == cell_idx);
                        if ~isempty(msn_idx) && length(od.msn_res.onTrack_spec) >= msn_idx
                            existing_results.msn = od.msn_res;
                            existing_results.cell_idx = msn_idx; % Store the MSN index, not the overall cell index
                        end
                    end
                elseif od.cell_type(cell_idx) == 2 % FSI
                    if isfield(od, 'fsi_res')
                        % Find which FSI this cell corresponds to among all FSIs
                        fsi_cells = find(od.cell_type == 2);
                        fsi_idx = find(fsi_cells == cell_idx);
                        if ~isempty(fsi_idx) && length(od.fsi_res.onTrack_spec) >= fsi_idx
                            existing_results.fsi = od.fsi_res;
                            existing_results.cell_idx = fsi_idx; % Store the FSI index, not the overall cell index
                        end
                    end
                end
            end
        end
    catch ME
        warning('Could not load existing results: %s', ME.message);
    end
end

%%
% Function to calculate new measures
function new_results = calculateNewMeasures(cfg_in, cfg_master, sd, iC, ft_csc, cfg_onTrack, ExpKeys)
    new_results = [];
    
    % Calculate and save STA
    cfg_ft.timwin = [-0.5 0.5];
    cfg_ft.spikechannel = sd.S.ft_spikes(iC).label{1};
    cfg_ft.channel = ft_csc.label(1);
    this_data = ft_appendspike([], ft_csc, sd.S.ft_spikes(iC));
    
    % Restrict data to only on-track data
    on_track_data = ft_redefinetrial(cfg_onTrack, this_data);
    
    % Calculate STA
    this_sta = ft_spiketriggeredaverage(cfg_ft, on_track_data);
    new_results.sta_time = this_sta.time;
    new_results.sta_vals = this_sta.avg(:,:)';
    new_results.spk_count = sum(on_track_data.trial{1}(2,:));
    
    % Calculate STS
    cfg_sts.method = 'mtmconvol';
    cfg_sts.foi = 1:1:100;
    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
    cfg_sts.taper = 'hanning';
    cfg_sts.spikechannel = sd.S.ft_spikes(iC).label{1};
    cfg_sts.channel = on_track_data.label(1);
    cfg_sts.rejectsaturation = 'no';
    this_sts = ft_spiketriggeredspectrum(cfg_sts, on_track_data);
    
    new_results.freqs = this_sts.freq;
    new_results.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
    
    % Calculate PPC with the specified method
    cfg_ppc = [];
    cfg_ppc.method = cfg_in.ppc_method;
    cfg_ppc.spikechannel = this_sts.label;
    cfg_ppc.channel = this_sts.lfplabel;
    cfg_ppc.avgoverchan = 'weighted';
    cfg_ppc.timwin = 'all';
    
    this_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, this_sts);
    
    % Store PPC results based on method
    if strcmp(cfg_in.ppc_method, 'ppc0')
        new_results.ppc = this_ppc.ppc0';
    elseif strcmp(cfg_in.ppc_method, 'ppc1')
        new_results.ppc = this_ppc.ppc1';
    elseif strcmp(cfg_in.ppc_method, 'ppc2')
        new_results.ppc = this_ppc.ppc2';
    elseif strcmp(cfg_in.ppc_method, 'plv')
        new_results.ppc = this_ppc.plv';
    elseif strcmp(cfg_in.ppc_method, 'ral')
        new_results.ppc = this_ppc.ral';
    elseif strcmp(cfg_in.ppc_method, 'ang')
        new_results.ppc = this_ppc.ang';
    elseif strcmp(cfg_in.ppc_method, 'angin')
        new_results.ppc = this_ppc.angin';
    elseif strcmp(cfg_in.ppc_method, 'angout')
        new_results.ppc = this_ppc.angout';
    else
        new_results.ppc = this_ppc.ppc0'; % default fallback
    end
    
     % Also calculate PLV and PPC0 for comparison
     cfg_plv = [];
     cfg_plv.method = 'plv';
     cfg_plv.spikechannel = this_sts.label;
     cfg_plv.channel = this_sts.lfplabel;
     cfg_plv.avgoverchan = 'weighted';
     cfg_plv.timwin = 'all';
     cfg_plv.rejectsaturation = 'no';
     
     plv_all = ft_spiketriggeredspectrum_stat(cfg_plv, this_sts);
     new_results.plv_all = plv_all.plv';
     
     % Calculate PPC0 for all trials (for overlay comparison)
     cfg_ppc0 = [];
     cfg_ppc0.method = 'ppc0';
     cfg_ppc0.spikechannel = this_sts.label;
     cfg_ppc0.channel = this_sts.lfplabel;
     cfg_ppc0.avgoverchan = 'weighted';
     cfg_ppc0.timwin = 'all';
     cfg_ppc0.rejectsaturation = 'no';
     
     ppc0_all = ft_spiketriggeredspectrum_stat(cfg_ppc0, this_sts);
     new_results.ppc0_all = ppc0_all.ppc0';
     
     % Now calculate for near reward trials
    new_results = calculateNearRewardMeasures(cfg_in, cfg_master, sd, iC, this_data, this_sts, cfg_ft, cfg_ppc, ExpKeys, new_results);
end

%%
% Function to calculate near reward measures
function new_results = calculateNearRewardMeasures(cfg_in, cfg_master, sd, iC, this_data, this_sts, cfg_ft, cfg_ppc, ExpKeys, new_results)
    
    % Get reward times and create trials
    rt1 = getRewardTimes();
    rt1 = rt1(rt1 > ExpKeys.TimeOnTrack);
    rt2 = getRewardTimes2();
    rt2 = rt2(rt2 > ExpKeys.TimeOnTrack);
    
    % Clean up reward times (same logic as original)
    rt_dif = diff(rt1);
    rt_dif = find(rt_dif <= 5);
    valid_rt1 = true(length(rt1),1);
    valid_rt2 = true(length(rt2),1);
    for i = 1:length(rt_dif)
        valid_rt1(rt_dif(i)) = false;
        valid_rt1(rt_dif(i)+1) = false;
        valid_rt2(rt2 >= rt1(rt_dif(i)) & rt2 <= rt1(rt_dif(i)+2)) = false;
    end
    
    rt_dif = diff(rt2);
    rt_dif = find(rt_dif <= 5);
    for i = 1:length(rt_dif)
        valid_rt2(rt_dif(i)) = false;
        valid_rt2(rt_dif(i)+1) = false;
        valid_rt1(rt1 >= rt2(rt_dif(i)-1) & rt1 <= rt2(rt_dif(i)+1)) = false;
    end
    
    rt1 = rt1(valid_rt1);
    rt2 = rt2(valid_rt2);
    
    if length(rt1) ~= length(rt2)
        rt1 = rt1(1:end-1);
    end
    
    keep = (rt1 <= rt2);
    rt1 = rt1(keep);
    rt2 = rt2(keep);
    
    % Create near reward trials
    w_start = rt1 - 5;
    w_end = rt2 + 5;
    w_end(end) = min(w_end(end), ExpKeys.TimeOffTrack);
    keep = ~isoutlier(w_end - w_start, 'median');
    w_start = w_start(keep);
    w_end = w_end(keep);
    rt_iv = iv(w_start, w_end);
    
    % Break down data into near trials
    temp_tvec = this_data.time{1} + double(this_data.hdr.FirstTimeStamp)/1e6;     
    temp_start = nearest_idx3(rt_iv.tstart, temp_tvec);
    temp_end = nearest_idx3(rt_iv.tend, temp_tvec);
    cfg_near_trials.trl = [temp_start, temp_end, zeros(size(temp_start))];
    near_data = ft_redefinetrial(cfg_near_trials, this_data);
    
    % Calculate firing rates and split into HFR/LFR
    tcount = length(near_data.trial);
    mfr = zeros(tcount, 1);
    all_tspikes = cell(1,tcount);
    for iT = 1:tcount
        all_tspikes{iT} = sum(near_data.trial{iT}(2,:));
        mfr(iT) = all_tspikes{iT}/near_data.time{iT}(end); 
    end
    
    % Find firing rate threshold
    nz_trials = find(mfr ~= 0);
    nz_mfr = mfr(nz_trials);
    nz_tspikes = cell2mat(all_tspikes(nz_trials));
    ufr = unique(nz_mfr);
    dif_min = sum(nz_tspikes);
    fr_thresh = 0;
    for iF = 1:length(ufr)
        cur_thresh = ufr(iF);
        l_spikes = sum(nz_tspikes(nz_mfr <= cur_thresh));
        h_spikes = sum(nz_tspikes(nz_mfr > cur_thresh));
        cur_dif = abs(h_spikes - l_spikes);
        if dif_min > cur_dif
            dif_min = cur_dif;
            fr_thresh = cur_thresh;
        end
    end
    
    hfr_trials = mfr > fr_thresh;
    lfr_trials = ~hfr_trials;
    
    % Calculate measures for HFR trials
    hfr_trl_idx = find(hfr_trials);
    if ~isempty(hfr_trl_idx)
        cfg_hfr.trl = cfg_near_trials.trl(hfr_trials,:);
        hfr_data = ft_redefinetrial(cfg_hfr, this_data);
        
        % STA
        hfr_sta = ft_spiketriggeredaverage(cfg_ft, hfr_data);
        new_results.near_hfr.sta_time = hfr_sta.time;
        new_results.near_hfr.sta_vals = hfr_sta.avg(:,:)';
        new_results.near_hfr.spk_count = sum(cell2mat(all_tspikes(hfr_trials)));
        
        % STS and PPC
        hfr_idx = [];
        for iT = 1:length(hfr_trl_idx)
            hfr_idx = [hfr_idx, find(near_data.trial{hfr_trl_idx(iT)}(2,:))];
        end
        
        % Ensure indices are within bounds
        if ~isempty(hfr_idx)
            max_idx = size(this_sts.fourierspctrm{1}, 1);
            hfr_idx = hfr_idx(hfr_idx <= max_idx & hfr_idx > 0);
        end
        
        hfr_sts = this_sts;
        hfr_sts.fourierspctrm{1} = hfr_sts.fourierspctrm{1}(hfr_idx,:,:);
        hfr_sts.time{1} = hfr_sts.time{1}(hfr_idx,:);
        hfr_sts.trial{1} = hfr_sts.trial{1}(hfr_idx,:);
        
        new_results.near_hfr.sts_vals = nanmean(sq(abs(hfr_sts.fourierspctrm{1})));
        
        hfr_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, hfr_sts);
        if strcmp(cfg_in.ppc_method, 'ppc0')
            new_results.near_hfr.ppc = hfr_ppc.ppc0';
        elseif strcmp(cfg_in.ppc_method, 'ppc1')
            new_results.near_hfr.ppc = hfr_ppc.ppc1';
        elseif strcmp(cfg_in.ppc_method, 'ppc2')
            new_results.near_hfr.ppc = hfr_ppc.ppc2';
        elseif strcmp(cfg_in.ppc_method, 'plv')
            new_results.near_hfr.ppc = hfr_ppc.plv';
        elseif strcmp(cfg_in.ppc_method, 'ral')
            new_results.near_hfr.ppc = hfr_ppc.ral';
        elseif strcmp(cfg_in.ppc_method, 'ang')
            new_results.near_hfr.ppc = hfr_ppc.ang';
        elseif strcmp(cfg_in.ppc_method, 'angin')
            new_results.near_hfr.ppc = hfr_ppc.angin';
        elseif strcmp(cfg_in.ppc_method, 'angout')
            new_results.near_hfr.ppc = hfr_ppc.angout';
        else
            new_results.near_hfr.ppc = hfr_ppc.ppc0';
        end
        
        % Also calculate PLV for HFR trials
        cfg_plv = [];
        cfg_plv.method = 'plv';
        cfg_plv.spikechannel = hfr_sts.label;
        cfg_plv.channel = hfr_sts.lfplabel;
        cfg_plv.avgoverchan = 'weighted';
        cfg_plv.timwin = 'all';
        cfg_plv.rejectsaturation = 'no';
        hfr_plv = ft_spiketriggeredspectrum_stat(cfg_plv, hfr_sts);
        new_results.near_hfr.plv = hfr_plv.plv';
        
        % Also calculate PPC0 for HFR trials
        cfg_ppc0 = [];
        cfg_ppc0.method = 'ppc0';
        cfg_ppc0.spikechannel = hfr_sts.label;
        cfg_ppc0.channel = hfr_sts.lfplabel;
        cfg_ppc0.avgoverchan = 'weighted';
        cfg_ppc0.timwin = 'all';
        cfg_ppc0.rejectsaturation = 'no';
        hfr_ppc0 = ft_spiketriggeredspectrum_stat(cfg_ppc0, hfr_sts);
        new_results.near_hfr.ppc0 = hfr_ppc0.ppc0';
    end
    
    % Calculate measures for LFR trials
    lfr_trl_idx = find(lfr_trials);
    if ~isempty(lfr_trl_idx)
        cfg_lfr.trl = cfg_near_trials.trl(lfr_trials,:);
        lfr_data = ft_redefinetrial(cfg_lfr, this_data);
        
        % STA
        lfr_sta = ft_spiketriggeredaverage(cfg_ft, lfr_data);
        new_results.near_lfr.sta_time = lfr_sta.time;
        new_results.near_lfr.sta_vals = lfr_sta.avg(:,:)';
        new_results.near_lfr.spk_count = sum(cell2mat(all_tspikes(lfr_trials)));
        
                 % STS and PPC
         lfr_idx = [];
         for iT = 1:length(lfr_trl_idx)
             lfr_idx = [lfr_idx, find(near_data.trial{lfr_trl_idx(iT)}(2,:))];
         end
         
         % Ensure indices are within bounds
         if ~isempty(lfr_idx)
             max_idx = size(this_sts.fourierspctrm{1}, 1);
             lfr_idx = lfr_idx(lfr_idx <= max_idx & lfr_idx > 0);
         end
         
         lfr_sts = this_sts;
        lfr_sts.fourierspctrm{1} = lfr_sts.fourierspctrm{1}(lfr_idx,:,:);
        lfr_sts.time{1} = lfr_sts.time{1}(lfr_idx,:);
        lfr_sts.trial{1} = lfr_sts.trial{1}(lfr_idx,:);
        
        new_results.near_lfr.sts_vals = nanmean(sq(abs(lfr_sts.fourierspctrm{1})));
        
        lfr_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, lfr_sts);
        if strcmp(cfg_in.ppc_method, 'ppc0')
            new_results.near_lfr.ppc = lfr_ppc.ppc0';
        elseif strcmp(cfg_in.ppc_method, 'ppc1')
            new_results.near_lfr.ppc = lfr_ppc.ppc1';
        elseif strcmp(cfg_in.ppc_method, 'ppc2')
            new_results.near_lfr.ppc = lfr_ppc.ppc2';
        elseif strcmp(cfg_in.ppc_method, 'plv')
            new_results.near_lfr.ppc = lfr_ppc.plv';
        elseif strcmp(cfg_in.ppc_method, 'ral')
            new_results.near_lfr.ppc = lfr_ppc.ral';
        elseif strcmp(cfg_in.ppc_method, 'ang')
            new_results.near_lfr.ppc = lfr_ppc.ang';
        elseif strcmp(cfg_in.ppc_method, 'angin')
            new_results.near_lfr.ppc = lfr_ppc.angin';
        elseif strcmp(cfg_in.ppc_method, 'angout')
            new_results.near_lfr.ppc = lfr_ppc.angout';
        else
            new_results.near_lfr.ppc = lfr_ppc.ppc0';
        end
        
        % Also calculate PLV for LFR trials
        cfg_plv = [];
        cfg_plv.method = 'plv';
        cfg_plv.spikechannel = lfr_sts.label;
        cfg_plv.channel = lfr_sts.lfplabel;
        cfg_plv.avgoverchan = 'weighted';
        cfg_plv.timwin = 'all';
        cfg_plv.rejectsaturation = 'no';
        lfr_plv = ft_spiketriggeredspectrum_stat(cfg_plv, lfr_sts);
        new_results.near_lfr.plv = lfr_plv.plv';
        
        % Also calculate PPC0 for LFR trials
        cfg_ppc0 = [];
        cfg_ppc0.method = 'ppc0';
        cfg_ppc0.spikechannel = lfr_sts.label;
        cfg_ppc0.channel = lfr_sts.lfplabel;
        cfg_ppc0.avgoverchan = 'weighted';
        cfg_ppc0.timwin = 'all';
        cfg_ppc0.rejectsaturation = 'no';
        lfr_ppc0 = ft_spiketriggeredspectrum_stat(cfg_ppc0, lfr_sts);
        new_results.near_lfr.ppc0 = lfr_ppc0.ppc0';
    end
    
    % Store firing rate info
    new_results.near.mfr = mfr;
    new_results.near.fr_thresh = fr_thresh;
    new_results.near.trial_spk_count = cell2mat(all_tspikes);
end

%%
% Function to plot cell comparison
function plotCellComparison(cfg_in, cell_label, existing_results, new_results)
    
    % Create figure
    fig = figure('WindowState', 'maximized');
    
    % Define colors
    c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
    c2 = [26/255 255/255 26/255]; % Green for HFR
    c3 = [0.7 0.7 0.7]; % Gray for all trials
    
    % Plot STA in row 1, column 1
     ax1 = subplot(2,3,1);
     hold on;
     if isfield(new_results, 'near_lfr') && isfield(new_results, 'near_hfr')
         plot(ax1, new_results.near_lfr.sta_time, new_results.near_lfr.sta_vals, 'Color', c1, 'LineWidth', 2);
         plot(ax1, new_results.near_hfr.sta_time, new_results.near_hfr.sta_vals, 'Color', c2, 'LineWidth', 2);
     end
     plot(ax1, new_results.sta_time, new_results.sta_vals, 'Color', c3, 'LineWidth', 1);
     
     ax1.Box = 'off';
     ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
     ax1.Title.String = 'STA';
     ax1.XLabel.String = 'Time (sec)';
     ax1.YLabel.String = 'Amplitude';
     ax1.Title.FontSize = 25;
     ax1.XAxis.FontSize = 18;
     ax1.YAxis.FontSize = 18;
     ax1.TickDir = 'out';
     
     % Plot STS in row 1, column 2
     ax2 = subplot(2,3,2);
     hold on;
     if isfield(new_results, 'near_lfr') && isfield(new_results, 'near_hfr')
         plot(ax2, new_results.freqs, new_results.near_lfr.sts_vals, 'Color', c1, 'LineWidth', 2);
         plot(ax2, new_results.freqs, new_results.near_hfr.sts_vals, 'Color', c2, 'LineWidth', 2);
     end
     plot(ax2, new_results.freqs, new_results.sts_vals, 'Color', c3, 'LineWidth', 1);
     
     ax2.Box = 'off';
     ax2.XLim = [cfg_in.min_freq 100];
     ax2.XTick = [cfg_in.min_freq 10 25 50 75 100];
     ax2.Title.String = 'STS';
     ax2.XLabel.String = 'Frequency (Hz)';
     ax2.YLabel.String = 'Power';
     ax2.Title.FontSize = 25;
     ax2.XAxis.FontSize = 18;
     ax2.YAxis.FontSize = 18;
     ax2.TickDir = 'out';
     
     % Plot PPC0 with overlay in row 1, column 3
     ax3 = subplot(2,3,3);
     hold on;
     if isfield(new_results, 'near_lfr') && isfield(new_results, 'near_hfr')
         plot(ax3, new_results.freqs, new_results.near_lfr.ppc0, 'Color', c1, 'LineWidth', 2);
         plot(ax3, new_results.freqs, new_results.near_hfr.ppc0, 'Color', c2, 'LineWidth', 2);
     end
     plot(ax3, new_results.freqs, new_results.ppc0_all, 'Color', c3, 'LineWidth', 1);
     
     ax3.Box = 'off';
     ax3.XLim = [cfg_in.min_freq 100];
     ax3.XTick = [cfg_in.min_freq 10 25 50 75 100];
     ax3.Title.String = 'PPC0';
     ax3.XLabel.String = 'Frequency (Hz)';
     ax3.YLabel.String = 'PPC0';
     ax3.Title.FontSize = 25;
     ax3.XAxis.FontSize = 18;
     ax3.YAxis.FontSize = 18;
     ax3.TickDir = 'out';
     
     % Plot PLV in row 2, column 3
     ax6 = subplot(2,3,6);
     hold on;
     if isfield(new_results, 'near_lfr') && isfield(new_results, 'near_hfr')
         plot(ax6, new_results.freqs, new_results.near_lfr.plv, 'Color', c1, 'LineWidth', 2);
         plot(ax6, new_results.freqs, new_results.near_hfr.plv, 'Color', c2, 'LineWidth', 2);
     end
     plot(ax6, new_results.freqs, new_results.plv_all, 'Color', c3, 'LineWidth', 1);
     
     ax6.Box = 'off';
     ax6.XLim = [cfg_in.min_freq 100];
     ax6.XTick = [cfg_in.min_freq 10 25 50 75 100];
     ax6.Title.String = 'PLV';
     ax6.XLabel.String = 'Frequency (Hz)';
     ax6.YLabel.String = 'PLV';
     ax6.Title.FontSize = 25;
     ax6.XAxis.FontSize = 18;
     ax6.YAxis.FontSize = 18;
     ax6.TickDir = 'out';
     
     % Row 2, columns 1 and 2 are empty
     ax4 = subplot(2,3,4);
     ax4.Visible = 'off';
     ax4.Box = 'off';
     
     ax5 = subplot(2,3,5);
     ax5.Visible = 'off';
     ax5.Box = 'off';
    
    % Overlay existing subsampled measures if available (only on PPC0 plot)
    if ~isempty(existing_results)
        overlayExistingMeasures(ax3, existing_results, new_results);
    end
    
    % Set legends at the end to avoid conflicts
    legend(ax1, {'LFR', 'HFR', 'All Trials'}, 'Location', 'best');
    legend(ax2, {'LFR', 'HFR', 'All Trials'}, 'Location', 'best');
    legend(ax3, {'LFR', 'HFR', 'All Trials'}, 'Location', 'best');
    legend(ax6, {'LFR', 'HFR', 'All Trials'}, 'Location', 'best');
    
    % Save figure
    [~, fp, ~] = fileparts(pwd);
    output_file = fullfile(cfg_in.output_dir, sprintf('%s_%s_%s_analysis.png', fp, cell_label, cfg_in.ppc_method));
    print(fig, '-dpng', '-r300', output_file);
    fprintf('Figure saved to: %s\n', output_file);
    
    close(fig);
end

%%
% Function to overlay existing subsampled measures
function overlayExistingMeasures(ax_ppc0, existing_results, new_results)
    
    % Define colors for existing results
    c_existing_all = [0 0 0];        % Black for all trials
    c_existing_hfr = [1 0.5 0];      % Orange for HFR
    c_existing_lfr = [0 1 1];        % Cyan for LFR
    
    % Overlay on PPC0 plot only
    if isfield(existing_results, 'fsi')
        fsi_idx = existing_results.cell_idx;
        if length(existing_results.fsi.onTrack_spec) >= fsi_idx
            % Overlay subsampled PPC0 for all trials
            if isfield(existing_results.fsi.onTrack_spec{fsi_idx}, 'subsampled_ppc')
                plot(ax_ppc0, new_results.freqs, existing_results.fsi.onTrack_spec{fsi_idx}.subsampled_ppc, ...
                    'Color', c_existing_all, 'LineStyle', '--', 'LineWidth', 1.5);
            end
            
            % Overlay HFR subsampled PPC0
            if isfield(existing_results.fsi, 'near_hfr_spec') && length(existing_results.fsi.near_hfr_spec) >= fsi_idx
                if isfield(existing_results.fsi.near_hfr_spec{fsi_idx}, 'subsampled_ppc')
                    plot(ax_ppc0, new_results.freqs, existing_results.fsi.near_hfr_spec{fsi_idx}.subsampled_ppc, ...
                        'Color', c_existing_hfr, 'LineStyle', ':', 'LineWidth', 1.5);
                end
            end
            
            % Overlay LFR subsampled PPC0
            if isfield(existing_results.fsi, 'near_lfr_spec') && length(existing_results.fsi.near_lfr_spec) >= fsi_idx
                if isfield(existing_results.fsi.near_lfr_spec{fsi_idx}, 'subsampled_ppc')
                    plot(ax_ppc0, new_results.freqs, existing_results.fsi.near_lfr_spec{fsi_idx}.subsampled_ppc, ...
                        'Color', c_existing_lfr, 'LineStyle', ':', 'LineWidth', 1.5);
                end
            end
        end
    elseif isfield(existing_results, 'msn')
        msn_idx = existing_results.cell_idx;
        if length(existing_results.msn.onTrack_spec) >= msn_idx
            % Overlay subsampled PPC0 for all trials
            if isfield(existing_results.msn.onTrack_spec{msn_idx}, 'subsampled_ppc')
                plot(ax_ppc0, new_results.freqs, existing_results.msn.onTrack_spec{msn_idx}.subsampled_ppc, ...
                    'Color', c_existing_all, 'LineStyle', '--', 'LineWidth', 1.5);
            end
            
            % Overlay HFR subsampled PPC0
            if isfield(existing_results.msn, 'near_hfr_spec') && length(existing_results.msn.near_hfr_spec) >= msn_idx
                if isfield(existing_results.msn.near_hfr_spec{msn_idx}, 'subsampled_ppc')
                    plot(ax_ppc0, new_results.freqs, existing_results.msn.near_hfr_spec{msn_idx}.subsampled_ppc, ...
                        'Color', c_existing_hfr, 'LineStyle', ':', 'LineWidth', 1.5);
                end
            end
            
            % Overlay LFR subsampled PPC0
            if isfield(existing_results.msn, 'near_lfr_spec') && length(existing_results.msn.near_lfr_spec) >= msn_idx
                if isfield(existing_results.msn.near_lfr_spec{msn_idx}, 'subsampled_ppc')
                    plot(ax_ppc0, new_results.freqs, existing_results.msn.near_lfr_spec{msn_idx}.subsampled_ppc, ...
                        'Color', c_existing_lfr, 'LineStyle', ':', 'LineWidth', 1.5);
                end
            end
        end
    end
    
    % Update legend for PPC0 plot
    if isfield(ax_ppc0, 'Legend')
        current_legend = ax_ppc0.Legend.String;
        new_legend = [current_legend, 'Subsampled All (black --)', 'Subsampled HFR (orange :)', 'Subsampled LFR (cyan :)'];
        ax_ppc0.Legend.String = new_legend;
    end
end

%%
% Other functions (keep the original LoadSpikesTarget function)
function S = LoadSpikesTarget(cfg_in)
    if ~isfield(cfg_in, 'Target') % no target specified, load them all
        S = LoadSpikes([]);
        return;
    end
    LoadExpKeys;
    % target specified, need to do some work
    % first see if this session has more than one target
    nTargets = length(ExpKeys.Target);
    if ~iscell(ExpKeys.Target) || (iscell(ExpKeys.Target) && length(ExpKeys.Target) == 1) % one target
        target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
        if isempty(target_idx)
            S = ts;
        else
            S = LoadSpikes([]);
        end
    else % multiple targets, assume TetrodeTargets exists
        please = []; please.getTTnumbers = 1;
        S = LoadSpikes([]);
        target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
        tt_num_keep = find(ExpKeys.TetrodeTargets == target_idx);
        keep = ismember(S.usr.tt_num, tt_num_keep);
        S = SelectTS([], S, keep);
    end
end