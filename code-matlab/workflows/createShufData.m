%% Script to create pooled data for sessions to be used for significance testing

% setup
clear;

in_dir = 'D:\vStr_oscillatory_switch_results\temp';
out_dir = 'D:\vStr_oscillatory_switch_results\temp2\';
spk_dt = 0.025; % interspike interval for surrogate spike train used for spike-triggered spectrum pool
ofiles = dir(in_dir);
for iS = 1:length(ofiles)
    if ~endsWith(ofiles(iS).name,'.mat')
        continue
    end
%     % manishm edit (only redo these for files where pooled sts is nan)
%     out_res = load(strcat(out_dir, ofiles(iS).name));
%     if ~any(isnan(out_res.od.pool_sts.fourierspctrm{1}(:)))
%         continue
%     end
%     clear out_res
    % Load File
    load(strcat(ofiles(iS).folder, '\', ofiles(iS).name)); % Load a particular session

    f1 = extractBefore(ofiles(iS).name,'_');
    toks = strsplit(f1, '-');
    path = strcat('E:\ADRLabData\', toks{1}, '\', strjoin(toks(1:4),'-'));
    cd(path);

     % Load CSC, ExpKeys and the cell
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
    
    cfg = [];
    sd.S = LoadSpikes(cfg);
    % Restrict spikes to only OnTrack
    sd.S = restrict(sd.S, iv(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack));
    % Read the spike files into field trip format
    sd.S.ft_spikes = ft_read_spike(sd.S.label{1});
    
    cfg_f.begsample = cfg_onTrack.begsample;
    cfg_f.endsample = cfg_onTrack.endsample;
    f_spike = sd.S.ft_spikes;
    f_spike.timestamp{1} = double(ft_csc.hdr.FirstTimeStamp) +  ...
        10^6*(ft_csc.time{1}(1):spk_dt:ft_csc.time{1}(end));
    f_data = ft_appendspike([], ft_csc,f_spike);
    f_data = ft_redefinetrial(cfg_f,f_data);
    fprintf("Number of fake spikes for session %s is %d\n", path, sum(f_data.trial{1}(2,:)));
    cfg_f = [];
    cfg_f.method = 'mtmconvol';
    cfg_f.foi = 1:1:100;
    cfg_f.t_ftimwin = 5./cfg_f.foi;
    cfg_f.taper = 'hanning';
    cfg_f.spikechannel =  f_spike.label{1};
    cfg_f.channel = f_data.label{1};
    cfg_f.rejectsaturation = 'no';
    od.pool_sts = ft_spiketriggeredspectrum(cfg_f, f_data);
    
    fn_out = strcat(out_dir, ofiles(iS).name);
    save(fn_out,'od');
    clear f_data cfg_f f_spike pool_sts
    
end