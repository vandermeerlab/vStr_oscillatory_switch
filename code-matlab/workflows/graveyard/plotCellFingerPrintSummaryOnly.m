%% Script to collect existing results and create summary plots (no individual fingerprints)
cd('D:\vStr_oscillatory_switch_results\temp3'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
odir = 'D:\vStr_oscillatory_switch_results\summary\';
unsampled_dir = 'D:\vStr_oscillatory_switch_results\temp4\';
clean_msn = 0;
clean_fsi = 0;
c1 = [75/255 0/255 146/255];  % Violet/Purple
c2 = [26/255 255/255 26/255]; % Green
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange
min_freq = 2; % Minimum frequency in Hertz
pl_thresh = 99; % Percentile threshold to establish significance of phase locking
diff_thresh = 95; % Percentile threshold to establish significance of PPC difference
peak_freq_tol_win = 5; % Window around peak frequency for diff analysis (Hz)

% Updated headers to include unsampled results
headers = {'label', 'lfr_min', 'lfr_max', 'lfr_mean', 'hfr_min', ...
    'hfr_max', 'hfr_mean', 'lfr_sts_peak', 'lfr_sts_diff', ...
    'hfr_sts_peak', 'hfr_sts_diff', 'lfr_ppc_peak', 'lfr_ppc_diff', ...
    'hfr_ppc_peak', 'hfr_ppc_diff', 'lfr_sts_peak_mag', 'hfr_sts_peak_mag', ...
    'lfr_ppc_peak_mag', 'hfr_ppc_peak_mag', 'lfr_sts_overall_peak', ...
    'hfr_sts_overall_peak', 'lfr_sts_overall_peak_mag', 'hfr_sts_overall_peak_mag', ...
    'unsampled_lfr_sts_peak', 'unsampled_lfr_sts_peak_mag', 'unsampled_hfr_sts_peak', 'unsampled_hfr_sts_peak_mag', ...
    'unsampled_lfr_ppc', 'unsampled_hfr_ppc', 'unsampled_lfr_plv', 'unsampled_hfr_plv'};

msn_summary = cell2table(cell(0,length(headers)), 'VariableNames', headers);
fsi_summary = cell2table(cell(0,length(headers)), 'VariableNames', headers);

%% Main processing loop - collect existing results only
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        
        % Process FSI cells
        if isfield(od, 'fsi_res')
            fsi_labels = od.label(od.cell_type == 2);
            fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
            [fsi_results_table, clean_fsi_count] = collect_cell_type_results_with_unsampled(od.fsi_res, fsi_labels, 'FSI', ...
                struct('c1',c1,'c2',c2,'c3',c3,'c4',c4), ...
                struct('min_freq',min_freq,'pl_thresh',pl_thresh,'diff_thresh',diff_thresh,'peak_freq_tol_win',peak_freq_tol_win,'odir',odir), ...
                headers, unsampled_dir);
            fsi_summary = [fsi_summary; fsi_results_table];
            clean_fsi = clean_fsi + clean_fsi_count;
        end
        
        % Process MSN cells  
        if isfield(od, 'msn_res')
            msn_labels = od.label(od.cell_type == 1);
            msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
            [msn_results_table, clean_msn_count] = collect_cell_type_results(od.msn_res, msn_labels, 'MSN', ...
                struct('c1',c1,'c2',c2,'c3',c3,'c4',c4), ...
                struct('min_freq',min_freq,'pl_thresh',pl_thresh,'diff_thresh',diff_thresh,'peak_freq_tol_win',peak_freq_tol_win,'odir',odir), ...
                headers);
            msn_summary = [msn_summary; msn_results_table];
            clean_msn = clean_msn + clean_msn_count;
        end
    end
end

fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);

writetable(msn_summary, strcat(odir,'msn_summary.csv'));
writetable(fsi_summary, strcat(odir,'fsi_summary.csv'));

%% Create side-by-side FSI STS peak analysis plots (subsampled vs unsampled)
create_fsi_sts_comparison_plots(fsi_summary, c1, c2, c3, c4);

%% Helper Functions

% Function to collect results from existing cell data with unsampled results for FSIs
function [summary_table, clean_count] = collect_cell_type_results_with_unsampled(cell_res, cell_labels, cell_type, colors, params, headers, unsampled_dir)
    % Extract parameters
    min_freq = params.min_freq;
    pl_thresh = params.pl_thresh;
    diff_thresh = params.diff_thresh;
    peak_freq_tol_win = params.peak_freq_tol_win;
    
    % Initialize results
    summary_results = cell(0, 31);  % Pre-allocate with correct number of columns (23 + 8 for unsampled)
    clean_count = 0;
    
    for iC = 1:length(cell_labels)
        if isfield(cell_res.near_spec{iC}, 'flag_no_control_split') && ~cell_res.near_spec{iC}.flag_no_control_split
            clean_count = clean_count + 1;
            
            % Get frequency indices and calculate firing rate stats
            x1 = find(round(cell_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
            this_freqs = cell_res.near_spec{iC}.freqs(x1:end);
            
            [firing_rate_stats] = calculate_firing_rate_stats(cell_res, iC);
            
            % Find STS peaks and magnitudes
            [sts_peaks, sts_peak_mags, sts_overall_peaks, sts_overall_peak_mags] = find_sts_peaks(cell_res, iC, x1, this_freqs, pl_thresh, cell_type);
            
            % Find STS diff peaks
            [sts_diff_peaks] = find_sts_diff_peaks(cell_res, iC, x1, this_freqs, diff_thresh, peak_freq_tol_win, sts_peaks, cell_type);
            
            % Find PPC peaks and magnitudes
            [ppc_peaks, ppc_peak_mags] = find_ppc_peaks(cell_res, iC, x1, this_freqs, pl_thresh, cell_type);
            
            % Find PPC diff peaks
            [ppc_diff_peaks] = find_ppc_diff_peaks(cell_res, iC, x1, this_freqs, diff_thresh, peak_freq_tol_win, ppc_peaks, cell_type);
            
            % Look for unsampled results for this FSI
            unsampled_data = load_unsampled_results_for_cell(cell_labels{iC}, unsampled_dir);
            
            % Compile results row with both subsampled and unsampled data
            this_row = {cell_labels{iC}, firing_rate_stats.lfr_min, firing_rate_stats.lfr_max, ...
                firing_rate_stats.lfr_mean, firing_rate_stats.hfr_min, firing_rate_stats.hfr_max, firing_rate_stats.hfr_mean, ...
                sts_peaks.lfr, sts_diff_peaks.lfr, sts_peaks.hfr, sts_diff_peaks.hfr, ...
                ppc_peaks.lfr, ppc_diff_peaks.lfr, ppc_peaks.hfr, ppc_diff_peaks.hfr, ...
                sts_peak_mags.lfr, sts_peak_mags.hfr, ppc_peak_mags.lfr, ppc_peak_mags.hfr, ...
                sts_overall_peaks.lfr, sts_overall_peaks.hfr, sts_overall_peak_mags.lfr, sts_overall_peak_mags.hfr, ...
                unsampled_data.lfr_sts_peak, unsampled_data.lfr_sts_peak_mag, ...
                unsampled_data.hfr_sts_peak, unsampled_data.hfr_sts_peak_mag, ...
                unsampled_data.lfr_ppc, unsampled_data.hfr_ppc, ...
                unsampled_data.lfr_plv, unsampled_data.hfr_plv};
            
            summary_results(end+1, :) = this_row;
        end
    end
    
    % Convert to table with consistent data types
    if ~isempty(summary_results)
        summary_table = create_consistent_table(summary_results, headers);
    else
        summary_table = cell2table(cell(0,length(headers)), 'VariableNames', headers);
    end
end

% Function to load unsampled results for a specific cell
function [unsampled_data] = load_unsampled_results_for_cell(cell_label, unsampled_dir)
    % Initialize with NaN values
    unsampled_data.lfr_sts_peak = NaN;
    unsampled_data.lfr_sts_peak_mag = NaN;
    unsampled_data.hfr_sts_peak = NaN;
    unsampled_data.hfr_sts_peak_mag = NaN;
    unsampled_data.lfr_ppc = NaN;
    unsampled_data.hfr_ppc = NaN;
    unsampled_data.lfr_plv = NaN;
    unsampled_data.hfr_plv = NaN;
    
    % Look for the corresponding _results.mat file
    searchString = strcat(unsampled_dir, cell_label, '_results.mat');
    if exist(searchString, 'file')
        try
            temp_load = load(searchString);
            
            % Check if the file contains new_results
            if isfield(temp_load, 'new_results')
                new_res = temp_load.new_results;
                
                % Extract LFR and HFR STS data
                if isfield(new_res, 'near_lfr') && isfield(new_res.near_lfr, 'sts_vals')
                    [unsampled_data.lfr_sts_peak_mag, peak_idx] = max(new_res.near_lfr.sts_vals);
                    unsampled_data.lfr_sts_peak = new_res.freqs(peak_idx);
                end
                
                if isfield(new_res, 'near_hfr') && isfield(new_res.near_hfr, 'sts_vals')
                    [unsampled_data.hfr_sts_peak_mag, peak_idx] = max(new_res.near_hfr.sts_vals);
                    unsampled_data.hfr_sts_peak = new_res.freqs(peak_idx);
                end
                
                % Extract PPC data
                if isfield(new_res, 'near_lfr') && isfield(new_res.near_lfr, 'ppc0')
                    unsampled_data.lfr_ppc = new_res.near_lfr.ppc0;
                end
                
                if isfield(new_res, 'near_hfr') && isfield(new_res.near_hfr, 'ppc0')
                    unsampled_data.hfr_ppc = new_res.near_hfr.ppc0;
                end
                
                % Extract PLV data
                if isfield(new_res, 'near_lfr') && isfield(new_res.near_lfr, 'plv')
                    unsampled_data.lfr_plv = new_res.near_lfr.plv;
                end
                
                if isfield(new_res, 'near_hfr') && isfield(new_res.near_hfr, 'plv')
                    unsampled_data.hfr_plv = new_res.near_hfr.plv;
                end
            end
        catch ME
            warning('Could not load unsampled results for %s: %s', cell_label, ME.message);
        end
    else
        fprintf('Unsampled results do not exist for %s\n', cell_label);
    end
end

% Function to collect results from existing cell data (no plotting) - for MSNs
function [summary_table, clean_count] = collect_cell_type_results(cell_res, cell_labels, cell_type, colors, params, headers)
    % Extract parameters
    min_freq = params.min_freq;
    pl_thresh = params.pl_thresh;
    diff_thresh = params.diff_thresh;
    peak_freq_tol_win = params.peak_freq_tol_win;
    
    % Initialize results
    summary_results = cell(0, 31);  % Pre-allocate with correct number of columns
    clean_count = 0;
    
    for iC = 1:length(cell_labels)
        if isfield(cell_res.near_spec{iC}, 'flag_no_control_split') && ~cell_res.near_spec{iC}.flag_no_control_split
            clean_count = clean_count + 1;
            
            % Get frequency indices and calculate firing rate stats
            x1 = find(round(cell_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
            this_freqs = cell_res.near_spec{iC}.freqs(x1:end);
            
            [firing_rate_stats] = calculate_firing_rate_stats(cell_res, iC);
            
            % Find STS peaks and magnitudes
            [sts_peaks, sts_peak_mags, sts_overall_peaks, sts_overall_peak_mags] = find_sts_peaks(cell_res, iC, x1, this_freqs, pl_thresh, cell_type);
            
            % Find STS diff peaks
            [sts_diff_peaks] = find_sts_diff_peaks(cell_res, iC, x1, this_freqs, diff_thresh, peak_freq_tol_win, sts_peaks, cell_type);
            
            % Find PPC peaks and magnitudes
            [ppc_peaks, ppc_peak_mags] = find_ppc_peaks(cell_res, iC, x1, this_freqs, pl_thresh, cell_type);
            
            % Find PPC diff peaks
            [ppc_diff_peaks] = find_ppc_diff_peaks(cell_res, iC, x1, this_freqs, diff_thresh, peak_freq_tol_win, ppc_peaks, cell_type);
            
            % Compile results row (MSNs don't have unsampled data, so fill with NaN)
            this_row = {cell_labels{iC}, firing_rate_stats.lfr_min, firing_rate_stats.lfr_max, ...
                firing_rate_stats.lfr_mean, firing_rate_stats.hfr_min, firing_rate_stats.hfr_max, firing_rate_stats.hfr_mean, ...
                sts_peaks.lfr, sts_diff_peaks.lfr, sts_peaks.hfr, sts_diff_peaks.hfr, ...
                ppc_peaks.lfr, ppc_diff_peaks.lfr, ppc_peaks.hfr, ppc_diff_peaks.hfr, ...
                sts_peak_mags.lfr, sts_peak_mags.hfr, ppc_peak_mags.lfr, ppc_peak_mags.hfr, ...
                sts_overall_peaks.lfr, sts_overall_peaks.hfr, sts_overall_peak_mags.lfr, sts_overall_peak_mags.hfr, ...
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; % 8 NaN values for unsampled data
            
            summary_results(end+1, :) = this_row;
        end
    end
    
    % Convert to table with consistent data types
    if ~isempty(summary_results)
        summary_table = create_consistent_table(summary_results, headers);
    else
        summary_table = cell2table(cell(0,length(headers)), 'VariableNames', headers);
    end
end

% Function to create side-by-side FSI STS comparison plots (subsampled vs unsampled)
function create_fsi_sts_comparison_plots(fsi_summary, c1, c2, c3, c4)
    % Calculate the contrast index and mean frequency for subsampled FSIs
    fsi_hfr_mag = fsi_summary.hfr_sts_overall_peak_mag;
    fsi_lfr_mag = fsi_summary.lfr_sts_overall_peak_mag;
    fsi_hfr_freq = fsi_summary.hfr_sts_overall_peak;
    fsi_lfr_freq = fsi_summary.lfr_sts_overall_peak;
    
    % Contrast index: (HFR - LFR) / (HFR + LFR)
    fsi_contrast = (fsi_hfr_mag - fsi_lfr_mag) ./ (fsi_hfr_mag + fsi_lfr_mag);
    % Frequency difference
    fsi_freq_diff = fsi_hfr_freq - fsi_lfr_freq;
    
    % Calculate the same for unsampled FSIs
    fsi_unsampled_hfr_mag = fsi_summary.unsampled_hfr_sts_peak_mag;
    fsi_unsampled_lfr_mag = fsi_summary.unsampled_lfr_sts_peak_mag;
    fsi_unsampled_hfr_freq = fsi_summary.unsampled_hfr_sts_peak;
    fsi_unsampled_lfr_freq = fsi_summary.unsampled_lfr_sts_peak;
    
    % Contrast index: (HFR - LFR) / (HFR + LFR)
    fsi_unsampled_contrast = (fsi_unsampled_hfr_mag - fsi_unsampled_lfr_mag) ./ (fsi_unsampled_hfr_mag + fsi_unsampled_lfr_mag);
    % Frequency difference
    fsi_unsampled_freq_diff = fsi_unsampled_hfr_freq - fsi_unsampled_lfr_freq;
    
    % Remove NaN values for plotting
    fsi_valid = ~isnan(fsi_contrast) & ~isnan(fsi_freq_diff);
    fsi_unsampled_valid = ~isnan(fsi_unsampled_contrast) & ~isnan(fsi_unsampled_freq_diff);
    
    % Create figure with side-by-side plots
    fig = figure('WindowState', 'maximized');
    
    % Left subplot - Subsampled FSIs
    subplot(1,2,1);
    scatter(fsi_freq_diff(fsi_valid), fsi_contrast(fsi_valid), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    ylabel('STS Magnitude Contrast ((HFR-LFR)/(HFR+LFR))');
    xlabel('Peak Frequency Diff (Hz)');
    title(sprintf('Subsampled FSI STS Peak Analysis (n=%d)', sum(fsi_valid)));
    ylim([-1 1]);
    xlim([-100,100]);
    
    % Right subplot - Unsampled FSIs
    subplot(1,2,2);
    scatter(fsi_unsampled_freq_diff(fsi_unsampled_valid), fsi_unsampled_contrast(fsi_unsampled_valid), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    ylabel('STS Magnitude Contrast ((HFR-LFR)/(HFR+LFR))');
    xlabel('Peak Frequency Diff (Hz)');
    title(sprintf('Unsampled FSI STS Peak Analysis (n=%d)', sum(fsi_unsampled_valid)));
    ylim([-1 1]);
    xlim([-100,100]);
    
    % Save the comparison figure
    print(fig, '-dpng', '-r300', 'D:\vStr_oscillatory_switch_results\summary\FSI_STS_comparison_subsampled_vs_unsampled');
    close;
end

% Helper function: Calculate firing rate statistics
function [stats] = calculate_firing_rate_stats(cell_res, iC)
    nz_trials = cell_res.near_spec{iC}.mfr > 0;
    lfr_trials = cell_res.near_spec{iC}.mfr <= cell_res.near_spec{iC}.fr_thresh;
    lfr_trials = lfr_trials & nz_trials;
    hfr_trials = ~lfr_trials & nz_trials;
    
    stats.lfr_min = min(cell_res.near_spec{iC}.mfr(lfr_trials));
    stats.lfr_max = max(cell_res.near_spec{iC}.mfr(lfr_trials));
    stats.lfr_mean = sum(cell_res.near_spec{iC}.trial_spk_count(lfr_trials))/ ...
        sum(cell_res.near_spec{iC}.trial_spk_count(lfr_trials)'./cell_res.near_spec{iC}.mfr(lfr_trials));
    stats.hfr_min = min(cell_res.near_spec{iC}.mfr(hfr_trials));
    stats.hfr_max = max(cell_res.near_spec{iC}.mfr(hfr_trials));
    stats.hfr_mean = sum(cell_res.near_spec{iC}.trial_spk_count(hfr_trials))/ ...
        sum(cell_res.near_spec{iC}.trial_spk_count(hfr_trials)'./cell_res.near_spec{iC}.mfr(hfr_trials));
end

% Helper function: Find STS peaks without plotting
function [peaks, peak_mags, overall_peaks, overall_peak_mags] = find_sts_peaks(cell_res, iC, x1, this_freqs, pl_thresh, cell_type)
    % Get STS data based on cell type
    if strcmp(cell_type, 'FSI')
        lfr_sts = cell_res.near_lfr_spec{iC}.subsampled_sts;
        hfr_sts = cell_res.near_hfr_spec{iC}.subsampled_sts;
    else % MSN
        lfr_sts = cell_res.near_lfr_spec{iC}.sts_vals;
        hfr_sts = cell_res.near_hfr_spec{iC}.sts_vals;
    end
    
    % Plot thresholds
    pct_sts = prctile(cell_res.near_spec{iC}.shuf_sts, pl_thresh);
    
    % Find threshold-crossing peaks
    lfr_sts_mask = lfr_sts > pct_sts;
    hfr_sts_mask = hfr_sts > pct_sts;
    
    % Initialize threshold-crossing outputs
    peaks.lfr = NaN; peaks.hfr = NaN;
    peak_mags.lfr = NaN; peak_mags.hfr = NaN;
    
    % Find LFR threshold-crossing peak
    if sum(lfr_sts_mask(x1:end)) > 0
        [~, peak_idx] = max(lfr_sts(x1:end).*(lfr_sts_mask(x1:end)));
        peaks.lfr = this_freqs(peak_idx);
        peak_mags.lfr = lfr_sts(x1-1+peak_idx);
    end
    
    % Find HFR threshold-crossing peak
    if sum(hfr_sts_mask(x1:end)) > 0
        [~, peak_idx] = max(hfr_sts(x1:end).*(hfr_sts_mask(x1:end)));
        peaks.hfr = this_freqs(peak_idx);
        peak_mags.hfr = hfr_sts(x1-1+peak_idx);
    end
    
    % Initialize overall peak outputs
    overall_peaks.lfr = NaN; overall_peaks.hfr = NaN;
    overall_peak_mags.lfr = NaN; overall_peak_mags.hfr = NaN;
    
    % Find LFR overall peak (maximum regardless of threshold)
    [~, peak_idx] = max(lfr_sts(x1:end));
    overall_peaks.lfr = this_freqs(peak_idx);
    overall_peak_mags.lfr = lfr_sts(x1-1+peak_idx);
    
    % Find HFR overall peak (maximum regardless of threshold)
    [~, peak_idx] = max(hfr_sts(x1:end));
    overall_peaks.hfr = this_freqs(peak_idx);
    overall_peak_mags.hfr = hfr_sts(x1-1+peak_idx);
end

% Helper function: Find STS diff peaks without plotting
function [diff_peaks] = find_sts_diff_peaks(cell_res, iC, x1, this_freqs, diff_thresh, peak_freq_tol_win, sts_peaks, cell_type)
    % Get diff data based on cell type
    if strcmp(cell_type, 'FSI')
        control_sts_diff = abs(cell_res.near_p1_spec{iC}.subsampled_sts - cell_res.near_p2_spec{iC}.subsampled_sts);
        sts_diff = cell_res.near_hfr_spec{iC}.subsampled_sts - cell_res.near_lfr_spec{iC}.subsampled_sts;
    else % MSN
        control_sts_diff = abs(cell_res.near_p1_spec{iC}.sts - cell_res.near_p2_spec{iC}.sts);
        sts_diff = cell_res.near_hfr_spec{iC}.sts_vals - cell_res.near_lfr_spec{iC}.sts_vals;
    end
    
    pct_sts_diff = prctile(control_sts_diff, diff_thresh);
    
    % Initialize outputs
    diff_peaks.lfr = NaN; diff_peaks.hfr = NaN;
    
    % Analyze LFR window
    if ~isnan(sts_peaks.lfr)
        lfr_window_mask = (this_freqs >= (sts_peaks.lfr - peak_freq_tol_win)) & ...
                          (this_freqs <= (sts_peaks.lfr + peak_freq_tol_win));
        
        thresh_diff = -pct_sts_diff(x1:end);
        this_diff = sts_diff(x1:end);
        thresh_diff(~lfr_window_mask) = NaN;
        this_diff(~lfr_window_mask) = NaN;
        
        cidx = this_diff < thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(thresh_diff - this_diff);
            diff_peaks.lfr = this_freqs(peak_idx);
        end
    end
    
    % Analyze HFR window
    if ~isnan(sts_peaks.hfr)
        hfr_window_mask = (this_freqs >= (sts_peaks.hfr - peak_freq_tol_win)) & ...
                          (this_freqs <= (sts_peaks.hfr + peak_freq_tol_win));
        
        thresh_diff = pct_sts_diff(x1:end);
        this_diff = sts_diff(x1:end);
        thresh_diff(~hfr_window_mask) = NaN;
        this_diff(~hfr_window_mask) = NaN;
        
        cidx = this_diff > thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(this_diff - thresh_diff);
            diff_peaks.hfr = this_freqs(peak_idx);
        end
    end
end

% Helper function: Find PPC peaks without plotting
function [peaks, peak_mags] = find_ppc_peaks(cell_res, iC, x1, this_freqs, pl_thresh, cell_type)
    % Get PPC data based on cell type
    if strcmp(cell_type, 'FSI')
        lfr_ppc = cell_res.near_lfr_spec{iC}.subsampled_ppc;
        hfr_ppc = cell_res.near_hfr_spec{iC}.subsampled_ppc;
    else % MSN
        lfr_ppc = cell_res.near_lfr_spec{iC}.ppc';
        hfr_ppc = cell_res.near_hfr_spec{iC}.ppc';
    end
    
    % Plot thresholds
    pct_ppc = prctile(cell_res.near_spec{iC}.shuf_ppc, pl_thresh);
    
    % Find peaks
    lfr_ppc_mask = lfr_ppc > pct_ppc;
    hfr_ppc_mask = hfr_ppc > pct_ppc;
    
    % Initialize outputs
    peaks.lfr = NaN; peaks.hfr = NaN;
    peak_mags.lfr = NaN; peak_mags.hfr = NaN;
    
    % Find LFR peak
    if sum(lfr_ppc_mask(x1:end)) > 0
        if strcmp(cell_type, 'FSI')
            [~, peak_idx] = max(lfr_ppc(x1:end).*(lfr_ppc_mask(x1:end)));
            peak_mags.lfr = lfr_ppc(x1-1+peak_idx);
        else % MSN
            [~, peak_idx] = max(cell_res.near_lfr_spec{iC}.ppc(x1:end)'.*(lfr_ppc_mask(x1:end)));
            peak_mags.lfr = cell_res.near_lfr_spec{iC}.ppc(x1-1+peak_idx);
        end
        peaks.lfr = this_freqs(peak_idx);
    end
    
    % Find HFR peak
    if sum(hfr_ppc_mask(x1:end)) > 0
        if strcmp(cell_type, 'FSI')
            [~, peak_idx] = max(hfr_ppc(x1:end).*(hfr_ppc_mask(x1:end)));
            peak_mags.hfr = hfr_ppc(x1-1+peak_idx);
        else % MSN
            [~, peak_idx] = max(cell_res.near_hfr_spec{iC}.ppc(x1:end)'.*(hfr_ppc_mask(x1:end)));
            peak_mags.hfr = cell_res.near_hfr_spec{iC}.ppc(x1-1+peak_idx);
        end
        peaks.hfr = this_freqs(peak_idx);
    end
end

% Helper function: Find PPC diff peaks without plotting
function [diff_peaks] = find_ppc_diff_peaks(cell_res, iC, x1, this_freqs, diff_thresh, peak_freq_tol_win, ppc_peaks, cell_type)
    % Get diff data based on cell type
    if strcmp(cell_type, 'FSI')
        control_ppc_diff = abs(cell_res.near_p1_spec{iC}.subsampled_ppc - cell_res.near_p2_spec{iC}.subsampled_ppc);
        ppc_diff = cell_res.near_hfr_spec{iC}.subsampled_ppc - cell_res.near_lfr_spec{iC}.subsampled_ppc;
    else % MSN
        control_ppc_diff = abs(cell_res.near_p1_spec{iC}.ppc - cell_res.near_p2_spec{iC}.ppc);
        ppc_diff = cell_res.near_hfr_spec{iC}.ppc' - cell_res.near_lfr_spec{iC}.ppc';
    end
    
    pct_ppc_diff = prctile(control_ppc_diff, diff_thresh);
    
    % Initialize outputs
    diff_peaks.lfr = NaN; diff_peaks.hfr = NaN;
    
    % Analyze LFR window
    if ~isnan(ppc_peaks.lfr)
        lfr_window_mask = (this_freqs >= (ppc_peaks.lfr - peak_freq_tol_win)) & ...
                          (this_freqs <= (ppc_peaks.lfr + peak_freq_tol_win));
        
        thresh_diff = -pct_ppc_diff(x1:end);
        this_diff = ppc_diff(x1:end);
        thresh_diff(~lfr_window_mask) = NaN;
        this_diff(~lfr_window_mask) = NaN;
        
        cidx = this_diff < thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(thresh_diff - this_diff);
            diff_peaks.lfr = this_freqs(peak_idx);
        end
    end
    
    % Analyze HFR window
    if ~isnan(ppc_peaks.hfr)
        hfr_window_mask = (this_freqs >= (ppc_peaks.hfr - peak_freq_tol_win)) & ...
                          (this_freqs <= (ppc_peaks.hfr + peak_freq_tol_win));
        
        thresh_diff = pct_ppc_diff(x1:end);
        this_diff = ppc_diff(x1:end);
        thresh_diff(~hfr_window_mask) = NaN;
        this_diff(~hfr_window_mask) = NaN;
        
        cidx = this_diff > thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(this_diff - thresh_diff);
            diff_peaks.hfr = this_freqs(peak_idx);
        end
    end
end

% Helper function: Create table with consistent data types
function [summary_table] = create_consistent_table(summary_results, headers)
    % Convert cell array to table first
    temp_table = cell2table(summary_results, 'VariableNames', headers);
    
    % Ensure unsampled columns are consistently double type
    unsampled_cols = {'unsampled_lfr_sts_peak', 'unsampled_lfr_sts_peak_mag', ...
                      'unsampled_hfr_sts_peak', 'unsampled_hfr_sts_peak_mag', ...
                      'unsampled_lfr_ppc', 'unsampled_hfr_ppc', ...
                      'unsampled_lfr_plv', 'unsampled_hfr_plv'};
    
    for i = 1:length(unsampled_cols)
        if ismember(unsampled_cols{i}, headers)
            col_idx = find(strcmp(headers, unsampled_cols{i}));
            if col_idx <= size(temp_table, 2)
                % Convert column to double, handling NaN values
                col_data = temp_table.(unsampled_cols{i});
                if iscell(col_data)
                    % Convert cell array to double array
                    double_data = zeros(size(col_data));
                    for j = 1:length(col_data)
                        if iscell(col_data{j})
                            double_data(j) = NaN;
                        elseif isnan(col_data{j})
                            double_data(j) = NaN;
                        elseif isvector(col_data{j}) && length(col_data{j}) > 1
                            % If it's a vector, take the mean (or max for peaks)
                            if contains(unsampled_cols{i}, 'peak')
                                % For peak values, take the maximum
                                double_data(j) = max(col_data{j}, [], 'omitnan');
                            else
                                % For other values, take the mean
                                double_data(j) = mean(col_data{j}, 'omitnan');
                            end
                        elseif ismatrix(col_data{j}) && numel(col_data{j}) > 1
                            % If it's a matrix, take the mean
                            double_data(j) = mean(col_data{j}(:), 'omitnan');
                        else
                            % Single value, convert to double
                            double_data(j) = double(col_data{j});
                        end
                    end
                    temp_table.(unsampled_cols{i}) = double_data;
                else
                    % Already numeric, ensure it's double and handle vectors/matrices
                    if isvector(col_data) && length(col_data) > 1
                        % If it's a vector, take the mean (or max for peaks) for each row
                        if contains(unsampled_cols{i}, 'peak')
                            % For peak values, take the maximum of each row
                            if size(col_data, 1) == 1
                                % Single row vector, take max
                                temp_table.(unsampled_cols{i}) = max(col_data, [], 'omitnan');
                            else
                                % Multiple rows, take max of each row
                                temp_table.(unsampled_cols{i}) = max(col_data, [], 2, 'omitnan');
                            end
                        else
                            % For other values, take the mean of each row
                            if size(col_data, 1) == 1
                                % Single row vector, take mean
                                temp_table.(unsampled_cols{i}) = mean(col_data, 'omitnan');
                            else
                                % Multiple rows, take mean of each row
                                temp_table.(unsampled_cols{i}) = mean(col_data, 2, 'omitnan');
                            end
                        end
                    elseif ismatrix(col_data) && numel(col_data) > 1
                        % If it's a matrix, take the mean of each row
                        if size(col_data, 1) == 1
                            % Single row, take mean of all elements
                            temp_table.(unsampled_cols{i}) = mean(col_data(:), 'omitnan');
                        else
                            % Multiple rows, take mean of each row
                            temp_table.(unsampled_cols{i}) = mean(col_data, 2, 'omitnan');
                        end
                    else
                        % Single value, ensure it's double
                        temp_table.(unsampled_cols{i}) = double(col_data);
                    end
                end
            end
        end
    end
    
    summary_table = temp_table;
end
