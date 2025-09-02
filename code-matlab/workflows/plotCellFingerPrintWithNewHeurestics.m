%% Script to generate cell_finger prints and summary 
cd('D:\vStr_oscillatory_switch_results\temp3'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
odir = 'D:\vStr_oscillatory_switch_results\cellFingerPrint\';
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

% Load artifact cell list
fid = fopen('D:\vStr_oscillatory_switch\sta_artifact_cells.txt', 'r');
artifact_cell_list = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
artifact_cell_list = artifact_cell_list{1};
% Remove any empty lines
artifact_cell_list = artifact_cell_list(~cellfun(@isempty, artifact_cell_list));

% Updated headers without overall peaks
headers = {'label', 'lfr_min', 'lfr_max', 'lfr_mean', 'hfr_min', ...
    'hfr_max', 'hfr_mean', 'lfr_sts_peak', 'lfr_sts_diff', ...
    'hfr_sts_peak', 'hfr_sts_diff', 'lfr_ppc_peak', 'lfr_ppc_diff', ...
    'hfr_ppc_peak', 'hfr_ppc_diff', 'lfr_sts_peak_mag', 'hfr_sts_peak_mag', ...
    'lfr_ppc_peak_mag', 'hfr_ppc_peak_mag'};

msn_summary = cell2table(cell(0,length(headers)), 'VariableNames', headers);
fsi_summary = cell2table(cell(0,length(headers)), 'VariableNames', headers);

%% Main processing loop
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
            [fsi_results_table, clean_fsi_count] = process_cell_type(od.fsi_res, fsi_labels, 'FSI', ...
                struct('c1',c1,'c2',c2,'c3',c3,'c4',c4), ...
                struct('min_freq',min_freq,'pl_thresh',pl_thresh,'diff_thresh',diff_thresh,'peak_freq_tol_win',peak_freq_tol_win,'odir',odir), ...
                headers, artifact_cell_list);
            fsi_summary = [fsi_summary; fsi_results_table];
            clean_fsi = clean_fsi + clean_fsi_count;
        end
        
        % Process MSN cells  
        if isfield(od, 'msn_res')
            msn_labels = od.label(od.cell_type == 1);
            msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
            [msn_results_table, clean_msn_count] = process_cell_type(od.msn_res, msn_labels, 'MSN', ...
                struct('c1',c1,'c2',c2,'c3',c3,'c4',c4), ...
                struct('min_freq',min_freq,'pl_thresh',pl_thresh,'diff_thresh',diff_thresh,'peak_freq_tol_win',peak_freq_tol_win,'odir',odir), ...
                headers, artifact_cell_list);
            msn_summary = [msn_summary; msn_results_table];
            clean_msn = clean_msn + clean_msn_count;
        end
    end
end

fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
save(strcat(odir,'msn_summary.mat'), 'msn_summary');
save(strcat(odir,'fsi_summary.mat'), 'fsi_summary');

%% Summary visualization
% Define colors
c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
c2 = [26/255 255/255 26/255]; % Green for HFR
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange

% Define masks for cells with significant phase locking in both HFR and LFR
only_hfr_msn = ~isnan(msn_summary.hfr_sts_peak);
only_lfr_msn = ~isnan(msn_summary.lfr_sts_peak);
clean_msn = only_hfr_msn & only_lfr_msn;

only_hfr_fsi = ~isnan(fsi_summary.hfr_sts_peak);
only_lfr_fsi = ~isnan(fsi_summary.lfr_sts_peak);
clean_fsi = only_hfr_fsi & only_lfr_fsi;

% Define masks for significant PPC differences
dif_hfr_msn = ~isnan(msn_summary.hfr_ppc_diff);
dif_lfr_msn = ~isnan(msn_summary.lfr_ppc_diff);
dif_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_diff);
dif_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_diff);

% Create figure
fig = figure('WindowState', 'maximized');
ax1 = gca;

% Plot MSNs
msn_delta_fr = msn_summary.hfr_mean(clean_msn) - msn_summary.lfr_mean(clean_msn);
msn_delta_freq = msn_summary.hfr_sts_peak(clean_msn) - msn_summary.lfr_sts_peak(clean_msn);

% MSNs with no significant PPC diff
no_diff_msn = clean_msn & ~dif_hfr_msn & ~dif_lfr_msn;
if sum(no_diff_msn) > 0
    s1 = scatter(ax1, msn_delta_fr(no_diff_msn), msn_delta_freq(no_diff_msn), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
end

% MSNs with only HFR significant PPC diff
only_hfr_diff_msn = clean_msn & dif_hfr_msn & ~dif_lfr_msn;
if sum(only_hfr_diff_msn) > 0
    s2 = scatter(ax1, msn_delta_fr(only_hfr_diff_msn), msn_delta_freq(only_hfr_diff_msn), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    % Add '|' marker
    text(msn_delta_fr(only_hfr_diff_msn), msn_delta_freq(only_hfr_diff_msn), '|', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
end

% MSNs with only LFR significant PPC diff
only_lfr_diff_msn = clean_msn & ~dif_hfr_msn & dif_lfr_msn;
if sum(only_lfr_diff_msn) > 0
    s3 = scatter(ax1, msn_delta_fr(only_lfr_diff_msn), msn_delta_freq(only_lfr_diff_msn), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    % Add '_' marker
    text(msn_delta_fr(only_lfr_diff_msn), msn_delta_freq(only_lfr_diff_msn), '_', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
end

% MSNs with both HFR and LFR significant PPC diff
both_diff_msn = clean_msn & dif_hfr_msn & dif_lfr_msn;
if sum(both_diff_msn) > 0
    s4 = scatter(ax1, msn_delta_fr(both_diff_msn), msn_delta_freq(both_diff_msn), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    % Add both '|' and '_' markers
    text(msn_delta_fr(both_diff_msn), msn_delta_freq(both_diff_msn), '|_', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
end

hold on;

% Plot FSIs
fsi_delta_fr = fsi_summary.hfr_mean(clean_fsi) - fsi_summary.lfr_mean(clean_fsi);
fsi_delta_freq = fsi_summary.hfr_sts_peak(clean_fsi) - fsi_summary.lfr_sts_peak(clean_fsi);

% FSIs with no significant PPC diff
no_diff_fsi = clean_fsi & ~dif_hfr_fsi & ~dif_lfr_fsi;
if sum(no_diff_fsi) > 0
    s5 = scatter(ax1, fsi_delta_fr(no_diff_fsi), fsi_delta_freq(no_diff_fsi), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
end

% FSIs with only HFR significant PPC diff
only_hfr_diff_fsi = clean_fsi & dif_hfr_fsi & ~dif_lfr_fsi;
if sum(only_hfr_diff_fsi) > 0
    s6 = scatter(ax1, fsi_delta_fr(only_hfr_diff_fsi), fsi_delta_freq(only_hfr_diff_fsi), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    % Add '|' marker
    text(fsi_delta_fr(only_hfr_diff_fsi), fsi_delta_freq(only_hfr_diff_fsi), '|', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
end

% FSIs with only LFR significant PPC diff
only_lfr_diff_fsi = clean_fsi & ~dif_hfr_fsi & dif_lfr_fsi;
if sum(only_lfr_diff_fsi) > 0
    s7 = scatter(ax1, fsi_delta_fr(only_lfr_diff_fsi), fsi_delta_freq(only_lfr_diff_fsi), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    % Add '_' marker
    text(fsi_delta_fr(only_lfr_diff_fsi), fsi_delta_freq(only_lfr_diff_fsi), '_', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
end

% FSIs with both HFR and LFR significant PPC diff
both_diff_fsi = clean_fsi & dif_hfr_fsi & dif_lfr_fsi;
if sum(both_diff_fsi) > 0
    s8 = scatter(ax1, fsi_delta_fr(both_diff_fsi), fsi_delta_freq(both_diff_fsi), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    % Add both '|' and '_' markers
    text(fsi_delta_fr(both_diff_fsi), fsi_delta_freq(both_diff_fsi), '|_', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
end

ax1.XLabel.String = '\Delta Firing Rate';
ax1.YLabel.String = '\Delta Peak Frequency (Hz)';
ax1.FontSize = 16;
ax1.Title.String = 'Delta Firing Rate vs Delta Peak Frequency';

% Create legend
legend_entries = {};
legend_handles = [];

if sum(clean_msn) > 0
    legend_entries{end+1} = sprintf('MSN: %d cells', sum(clean_msn));
    legend_handles(end+1) = scatter(NaN, NaN, 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
end

if sum(clean_fsi) > 0
    legend_entries{end+1} = sprintf('FSI: %d cells', sum(clean_fsi));
    legend_handles(end+1) = scatter(NaN, NaN, 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
end

if sum(only_hfr_diff_msn) > 0 || sum(only_hfr_diff_fsi) > 0
    legend_entries{end+1} = '| = HFR PPC diff only';
    legend_handles(end+1) = text(NaN, NaN, '|', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
end

if sum(only_lfr_diff_msn) > 0 || sum(only_lfr_diff_fsi) > 0
    legend_entries{end+1} = '_ = LFR PPC diff only';
    legend_handles(end+1) = text(NaN, NaN, '_', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
end

if sum(both_diff_msn) > 0 || sum(both_diff_fsi) > 0
    legend_entries{end+1} = '|_ = Both PPC diffs';
    legend_handles(end+1) = text(NaN, NaN, '|_', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
end

if ~isempty(legend_entries)
    legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 14);
end

%% Other functions
% Modular function to process each cell type
function [summary_table, clean_count] = process_cell_type(cell_res, cell_labels, cell_type, colors, params, headers, artifact_cell_list)
    % Extract colors
    c1 = colors.c1; c2 = colors.c2; c3 = colors.c3; c4 = colors.c4;
    
    % Extract parameters
    min_freq = params.min_freq;
    pl_thresh = params.pl_thresh;
    diff_thresh = params.diff_thresh;
    peak_freq_tol_win = params.peak_freq_tol_win;
    odir = params.odir;
    
    % Initialize results
    summary_results = cell(0, 19);  % Pre-allocate with correct number of columns (removed overall peaks)
    clean_count = 0;
    
    for iC = 1:length(cell_labels)
        % Check if cell should be rejected due to artifacts
        cell_label = cell_labels{iC};
        if any(contains(artifact_cell_list, cell_label))
            continue; % Skip this cell
        end
        
        if isfield(cell_res.near_spec{iC}, 'flag_no_control_split') && ~cell_res.near_spec{iC}.flag_no_control_split
            clean_count = clean_count + 1;
            
            % Get frequency indices and calculate firing rate stats
            x1 = find(round(cell_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
            this_freqs = cell_res.near_spec{iC}.freqs(x1:end);
            
            [firing_rate_stats] = calculate_firing_rate_stats(cell_res, iC);
            
            % Create figure
            fig = figure('WindowState', 'maximized');
            
            % Plot STA
            plot_sta(cell_res, iC, c1, c2, c3);
            
            % Plot STS and find peaks
            [sts_peaks, sts_peak_mags] = plot_sts_and_find_peaks(cell_res, iC, c1, c2, c3, x1, this_freqs, pl_thresh, cell_type);
            
            % Plot STS diff with windowing
            [sts_diff_peaks] = plot_sts_diff_with_windowing(cell_res, iC, c1, c2, c3, x1, this_freqs, diff_thresh, peak_freq_tol_win, sts_peaks, cell_type);
            
            % Plot PPC and find peaks
            [ppc_peaks, ppc_peak_mags] = plot_ppc_and_find_peaks(cell_res, iC, c1, c2, c3, x1, this_freqs, pl_thresh, firing_rate_stats, cell_type);
            
            % Plot PPC diff with windowing
            [ppc_diff_peaks] = plot_ppc_diff_with_windowing(cell_res, iC, c1, c2, c3, x1, this_freqs, diff_thresh, peak_freq_tol_win, ppc_peaks, cell_type);
            
            % Add cell label
            ax = subplot(2,3,4);
            text(0.2,0.5, strcat(cell_labels{iC}, '_', cell_type), 'Interpreter', 'None', 'FontSize', 14);
            ax.Box = 'off';
            ax.Visible = 'off';
            
            % Save figure
            print(fig, '-dpng', '-r300', strcat(odir, cell_labels{iC}, '_', cell_type));
            close;
            
            % Compile results row (removed overall peaks)
            this_row = {cell_labels{iC}, firing_rate_stats.lfr_min, firing_rate_stats.lfr_max, ...
                firing_rate_stats.lfr_mean, firing_rate_stats.hfr_min, firing_rate_stats.hfr_max, firing_rate_stats.hfr_mean, ...
                sts_peaks.lfr, sts_diff_peaks.lfr, sts_peaks.hfr, sts_diff_peaks.hfr, ...
                ppc_peaks.lfr, ppc_diff_peaks.lfr, ppc_peaks.hfr, ppc_diff_peaks.hfr, ...
                sts_peak_mags.lfr, sts_peak_mags.hfr, ppc_peak_mags.lfr, ppc_peak_mags.hfr};
            
            summary_results(end+1, :) = this_row;
        end
    end
    
    % Convert to table
    if ~isempty(summary_results)
        summary_table = cell2table(summary_results, 'VariableNames', headers);
    else
        summary_table = cell2table(cell(0,length(headers)), 'VariableNames', headers);
    end
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

% Helper function: Plot STA
function plot_sta(cell_res, iC, c1, c2, c3)
    ax1 = subplot(2,3,1);
    hold on;
    plot(ax1, cell_res.near_spec{iC}.sta_time, cell_res.near_lfr_spec{iC}.sta_vals, 'Color', c1);
    plot(ax1, cell_res.near_spec{iC}.sta_time, cell_res.near_hfr_spec{iC}.sta_vals, 'Color', c2);
    plot(ax1, cell_res.near_spec{iC}.sta_time, cell_res.near_spec{iC}.sta_vals, 'Color', c3);
    
    ax1.Box = 'off';
    ax1.YTick = [];
    ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
    ax1.Title.String = 'STA';
    ax1.XLabel.String = 'Time (sec)';
    ax1.Title.FontSize = 25;
    ax1.XAxis.FontSize = 18;
    ax1.YAxis.FontSize = 18;
    ax1.XAxis.FontWeight = 'normal';
    ax1.YAxis.FontWeight = 'normal';
    ax1.TickDir = 'out';
    
    legend({sprintf('LFR trials: %d spikes',cell_res.near_lfr_spec{iC}.spk_count), ...
        sprintf('HFR trials: %d spikes',cell_res.near_hfr_spec{iC}.spk_count), ...
        sprintf('All Trials: %d spikes',cell_res.near_spec{iC}.spk_count)}, 'Location','best');
end

% Helper function: Plot STS and find peaks (removed overall peaks)
function [peaks, peak_mags] = plot_sts_and_find_peaks(cell_res, iC, c1, c2, c3, x1, this_freqs, pl_thresh, cell_type)
    ax2 = subplot(2,3,2);
    hold on;
    
    % Get STS data based on cell type
    if strcmp(cell_type, 'FSI')
        lfr_sts = cell_res.near_lfr_spec{iC}.subsampled_sts;
        hfr_sts = cell_res.near_hfr_spec{iC}.subsampled_sts;
        all_sts = cell_res.near_spec{iC}.subsampled_sts;
    else % MSN
        lfr_sts = cell_res.near_lfr_spec{iC}.sts_vals;
        hfr_sts = cell_res.near_hfr_spec{iC}.sts_vals;
        all_sts = cell_res.near_spec{iC}.sts_vals;
    end
    
    plot(ax2, cell_res.near_spec{iC}.freqs, lfr_sts, 'Color', c1);
    plot(ax2, cell_res.near_spec{iC}.freqs, hfr_sts, 'Color', c2);
    plot(ax2, cell_res.near_spec{iC}.freqs, all_sts, 'Color', c3);
    
    % Plot thresholds
    pct_sts = prctile(cell_res.near_spec{iC}.shuf_sts, pl_thresh);
    plot(ax2, cell_res.near_spec{iC}.freqs, pct_sts, '--black');
    
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
        peak_mags.lfr = lfr_sts(x1-1+peak_idx);  % Get magnitude at peak
        xline(ax2, peaks.lfr, 'Color', c1, 'LineStyle', '--');
        % Add horizontal line from peak to y-axis (solid dash for threshold-crossing)
        yline(ax2, peak_mags.lfr, 'Color', c1, 'LineStyle', '--', 'Alpha', 0.7);
    end
    
    % Find HFR threshold-crossing peak
    if sum(hfr_sts_mask(x1:end)) > 0
        [~, peak_idx] = max(hfr_sts(x1:end).*(hfr_sts_mask(x1:end)));
        peaks.hfr = this_freqs(peak_idx);
        peak_mags.hfr = hfr_sts(x1-1+peak_idx);  % Get magnitude at peak
        xline(ax2, peaks.hfr, 'Color', c2, 'LineStyle', '--');
        % Add horizontal line from peak to y-axis (solid dash for threshold-crossing)
        yline(ax2, peak_mags.hfr, 'Color', c2, 'LineStyle', '--', 'Alpha', 0.7);
    end
    
    ax2.Box = 'off';
    ax2.YTick = [];
    ax2.XLim = [x1 100];
    ax2.XTick = [x1 10 25 50 75 100];
    ax2.Title.String = 'STS';
    ax2.XLabel.String = 'Frequency (Hz)';
    ax2.Title.FontSize = 25;
    ax2.XAxis.FontSize = 18;
    ax2.XAxis.FontWeight = 'normal';
    ax2.TickDir = 'out';
    legend({'LFR','HFR','All Trials', sprintf('%2d thresh', pl_thresh)}, 'Location','best');
end

% Helper function: Plot STS diff with windowing
function [diff_peaks] = plot_sts_diff_with_windowing(cell_res, iC, c1, c2, c3, x1, this_freqs, diff_thresh, peak_freq_tol_win, sts_peaks, cell_type)
    ax2 = subplot(2,3,5);
    
    % Get diff data based on cell type
    if strcmp(cell_type, 'FSI')
        control_sts_diff = abs(cell_res.near_p1_spec{iC}.subsampled_sts - cell_res.near_p2_spec{iC}.subsampled_sts);
        sts_diff = cell_res.near_hfr_spec{iC}.subsampled_sts - cell_res.near_lfr_spec{iC}.subsampled_sts;
    else % MSN
        control_sts_diff = abs(cell_res.near_p1_spec{iC}.sts - cell_res.near_p2_spec{iC}.sts);
        sts_diff = cell_res.near_hfr_spec{iC}.sts_vals - cell_res.near_lfr_spec{iC}.sts_vals;
    end
    
    pct_sts_diff = prctile(control_sts_diff, diff_thresh);
    
    plot(ax2, this_freqs, sts_diff(x1:end), 'Color', 'black');
    hold on
    plot(ax2, this_freqs, pct_sts_diff(x1:end), 'Color', c2, 'LineStyle','--');
    plot(ax2, this_freqs, -pct_sts_diff(x1:end), 'Color', c1, 'LineStyle','--');
    
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
        plot(ax2, this_freqs, thresh_diff, 'Color', c1, 'LineWidth', 2);
        
        cidx = this_diff < thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(thresh_diff - this_diff);
            diff_peaks.lfr = this_freqs(peak_idx);
            xline(ax2, diff_peaks.lfr, 'Color', c1, 'LineStyle', '--', 'LineWidth', 2);
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
        plot(ax2, this_freqs, thresh_diff, 'Color', c2, 'LineWidth', 2);
        
        cidx = this_diff > thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(this_diff - thresh_diff);
            diff_peaks.hfr = this_freqs(peak_idx);
            xline(ax2, diff_peaks.hfr, 'Color', c2, 'LineStyle', '--', 'LineWidth', 2);
        end
    end
    
    ax2.Box = 'off';
    ax2.YTick = [];
    ax2.XLim = [x1 100];
    ax2.XTick = [x1 10 25 50 75 100];
    ax2.Title.String = 'STS diff';
    ax2.XLabel.String = 'Frequency (Hz)';
    ax2.Title.FontSize = 25;
    ax2.XAxis.FontSize = 18;
    ax2.XAxis.FontWeight = 'normal';
    ax2.TickDir = 'out';
    legend({'HFR-LFR','HFR>LFR','LFR>HFR'}, 'Location','best');
end

% Helper function: Plot PPC and find peaks
function [peaks, peak_mags] = plot_ppc_and_find_peaks(cell_res, iC, c1, c2, c3, x1, this_freqs, pl_thresh, firing_rate_stats, cell_type)
    ax3 = subplot(2,3,3);
    hold on;
    
    % Get PPC data based on cell type
    if strcmp(cell_type, 'FSI')
        lfr_ppc = cell_res.near_lfr_spec{iC}.subsampled_ppc;
        hfr_ppc = cell_res.near_hfr_spec{iC}.subsampled_ppc;
        all_ppc = cell_res.near_spec{iC}.subsampled_ppc;
    else % MSN
        lfr_ppc = cell_res.near_lfr_spec{iC}.ppc';
        hfr_ppc = cell_res.near_hfr_spec{iC}.ppc';
        all_ppc = cell_res.near_spec{iC}.ppc';
    end
    
    plot(ax3, cell_res.near_spec{iC}.freqs, lfr_ppc, 'Color', c1);
    plot(ax3, cell_res.near_spec{iC}.freqs, hfr_ppc, 'Color', c2);
    plot(ax3, cell_res.near_spec{iC}.freqs, all_ppc, 'Color', c3);
    
    % Plot thresholds
    pct_ppc = prctile(cell_res.near_spec{iC}.shuf_ppc, pl_thresh);
    if any(isnan(pct_ppc))
        dummy = 1;
    end
    plot(ax3, cell_res.near_spec{iC}.freqs, pct_ppc, '--black');
    
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
        xline(ax3, peaks.lfr, 'Color', c1, 'LineStyle', '--');
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
        xline(ax3, peaks.hfr, 'Color', c2, 'LineStyle', '--');
    end
    
    ax3.Box = 'off';
    ax3.XLim = [x1 100];
    ax3.XTick = [x1 10 25 50 75 100];
    ax3.Title.String = 'PPC';
    ax3.XLabel.String = 'Frequency (Hz)';
    ax3.Title.FontSize = 25;
    ax3.XAxis.FontSize = 18;
    ax3.XAxis.FontWeight = 'normal';
    ax3.YAxis.FontSize = 18;
    ax3.YAxis.FontWeight = 'normal';
    ax3.TickDir = 'out';
    legend({sprintf('%.2f spks/s', firing_rate_stats.lfr_mean), sprintf('%.2f spks/s', firing_rate_stats.hfr_mean), ...
        'All Trials', sprintf('%2d thresh', pl_thresh)}, 'Location','best');
end

% Helper function: Plot PPC diff with windowing
function [diff_peaks] = plot_ppc_diff_with_windowing(cell_res, iC, c1, c2, c3, x1, this_freqs, diff_thresh, peak_freq_tol_win, ppc_peaks, cell_type)
    ax3 = subplot(2,3,6);
    
    % Get diff data based on cell type
    if strcmp(cell_type, 'FSI')
        control_ppc_diff = abs(cell_res.near_p1_spec{iC}.subsampled_ppc - cell_res.near_p2_spec{iC}.subsampled_ppc);
        ppc_diff = cell_res.near_hfr_spec{iC}.subsampled_ppc - cell_res.near_lfr_spec{iC}.subsampled_ppc;
    else % MSN
        control_ppc_diff = abs(cell_res.near_p1_spec{iC}.ppc - cell_res.near_p2_spec{iC}.ppc);
        ppc_diff = cell_res.near_hfr_spec{iC}.ppc' - cell_res.near_lfr_spec{iC}.ppc';
    end
    
    pct_ppc_diff = prctile(control_ppc_diff, diff_thresh);
    if any(isnan(pct_ppc_diff))
        dummy = 1;
    end
    
    plot(ax3, this_freqs, ppc_diff(x1:end), 'Color', 'black');
    hold on
    plot(ax3, this_freqs, pct_ppc_diff(x1:end), 'Color', c2, 'LineStyle','--');
    plot(ax3, this_freqs, -pct_ppc_diff(x1:end), 'Color', c1, 'LineStyle','--');
    
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
        plot(ax3, this_freqs, thresh_diff, 'Color', c1, 'LineWidth', 2);
        
        cidx = this_diff < thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(thresh_diff - this_diff);
            diff_peaks.lfr = this_freqs(peak_idx);
            xline(ax3, diff_peaks.lfr, 'Color', c1, 'LineStyle', '--', 'LineWidth', 2);
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
        plot(ax3, this_freqs, thresh_diff, 'Color', c2, 'LineWidth', 2);
        
        cidx = this_diff > thresh_diff;
        if sum(cidx) ~= 0
            this_diff(~cidx) = NaN;
            thresh_diff(~cidx) = NaN;
            [~, peak_idx] = max(this_diff - thresh_diff);
            diff_peaks.hfr = this_freqs(peak_idx);
            xline(ax3, diff_peaks.hfr, 'Color', c2, 'LineStyle', '--', 'LineWidth', 2);
        end
    end
    
    ax3.Box = 'off';
    ax3.YTick = [];
    ax3.XLim = [x1 100];
    ax3.XTick = [x1 10 25 50 75 100];
    ax3.Title.String = 'PPC diff';
    ax3.XLabel.String = 'Frequency (Hz)';
    ax3.Title.FontSize = 25;
    ax3.XAxis.FontSize = 18;
    ax3.XAxis.FontWeight = 'normal';
    ax3.TickDir = 'out';
    legend({'HFR-LFR','HFR>LFR','LFR>HFR'}, 'Location','best');
end
