%% Script to generate cell_finger prints with correlation analysis
cd('D:\vStr_oscillatory_switch_results\temp3'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
odir = 'D:\vStr_oscillatory_switch_results\cellFingerPrintsWithCorrelation\';
clean_msn = 0;
clean_fsi = 0;
c1 = [75/255 0/255 146/255];  % Violet/Purple
c2 = [26/255 255/255 26/255]; % Green
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange
min_freq = 2; % Minimum frequency in Hertz

% Load artifact cell list
fid = fopen('D:\vStr_oscillatory_switch\sta_artifact_cells.txt', 'r');
artifact_cell_list = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
artifact_cell_list = artifact_cell_list{1};
% Remove any empty lines
artifact_cell_list = artifact_cell_list(~cellfun(@isempty, artifact_cell_list));

% Updated headers for correlation analysis
headers = {'label', 'lfr_min', 'lfr_max', 'lfr_mean', 'hfr_min', ...
    'hfr_max', 'hfr_mean', 'sts_lfr_hfr_corr', 'ppc_lfr_hfr_corr', ...
    'sts_control_corr_mean', 'sts_control_corr_std', 'ppc_control_corr_mean', 'ppc_control_corr_std', ...
    'sts_control_corrs', 'ppc_control_corrs'};

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
                         [fsi_results_table, clean_fsi_count] = process_cell_type_with_correlation(od.fsi_res, fsi_labels, 'FSI', ...
                 struct('c1',c1,'c2',c2,'c3',c3,'c4',c4), ...
                 struct('min_freq',min_freq,'odir',odir), ...
                 headers, artifact_cell_list);
            fsi_summary = [fsi_summary; fsi_results_table];
            clean_fsi = clean_fsi + clean_fsi_count;
        end
        
        % Process MSN cells  
        if isfield(od, 'msn_res')
            msn_labels = od.label(od.cell_type == 1);
            msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
                         [msn_results_table, clean_msn_count] = process_cell_type_with_correlation(od.msn_res, msn_labels, 'MSN', ...
                 struct('c1',c1,'c2',c2,'c3',c3,'c4',c4), ...
                 struct('min_freq',min_freq,'odir',odir), ...
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


%% Other functions
% Modular function to process each cell type with correlation analysis
function [summary_table, clean_count] = process_cell_type_with_correlation(cell_res, cell_labels, cell_type, colors, params, headers, artifact_cell_list)
    % Extract colors
    c1 = colors.c1; c2 = colors.c2; c3 = colors.c3;
    
    % Extract parameters
    min_freq = params.min_freq;
    odir = params.odir;
    
    % Initialize results
    summary_results = cell(0, 15);  % Pre-allocate with correct number of columns
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
            
            % Calculate correlations
            [corr_stats] = calculate_correlations(cell_res, iC, x1, cell_type);
            
            % Create figure
            fig = figure('WindowState', 'maximized');
            
            % Add cell label at the top
            sgtitle(strcat(cell_labels{iC}, '_', cell_type), 'Interpreter', 'None', 'FontSize', 16, 'FontWeight', 'bold');
            
            % Plot STA
            plot_sta(cell_res, iC, c1, c2, c3);
            
            % Plot STS with correlation
            plot_sts_with_correlation(cell_res, iC, c1, c2, c3, x1, this_freqs, cell_type, corr_stats);
            
            % Plot PPC with correlation
            plot_ppc_with_correlation(cell_res, iC, c1, c2, c3, x1, this_freqs, firing_rate_stats, cell_type, corr_stats);
            
            % Save figure
            print(fig, '-dpng', '-r300', strcat(odir, cell_labels{iC}, '_', cell_type));
            close;
            
            % Compile results row
            this_row = {cell_labels{iC}, firing_rate_stats.lfr_min, firing_rate_stats.lfr_max, ...
                firing_rate_stats.lfr_mean, firing_rate_stats.hfr_min, firing_rate_stats.hfr_max, firing_rate_stats.hfr_mean, ...
                corr_stats.sts_lfr_hfr_corr, corr_stats.ppc_lfr_hfr_corr, ...
                corr_stats.sts_control_corr_mean, corr_stats.sts_control_corr_std, ...
                corr_stats.ppc_control_corr_mean, corr_stats.ppc_control_corr_std, ...
                {corr_stats.sts_control_corrs}, {corr_stats.ppc_control_corrs}};
            
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

% Helper function: Calculate correlations
function [corr_stats] = calculate_correlations(cell_res, iC, x1, cell_type)
    % Get data based on cell type
    if strcmp(cell_type, 'FSI')
        lfr_sts = cell_res.near_lfr_spec{iC}.subsampled_sts(x1:end);
        hfr_sts = cell_res.near_hfr_spec{iC}.subsampled_sts(x1:end);
        lfr_ppc = cell_res.near_lfr_spec{iC}.subsampled_ppc(x1:end);
        hfr_ppc = cell_res.near_hfr_spec{iC}.subsampled_ppc(x1:end);
        
        % Control partition correlations for STS
        sts_control_corrs = zeros(100, 1);
        ppc_control_corrs = zeros(100, 1);
        for p = 1:100
            p1_sts = cell_res.near_p1_spec{iC}.subsampled_sts(p,x1:end);
            p2_sts = cell_res.near_p2_spec{iC}.subsampled_sts(p,x1:end);
            p1_ppc = cell_res.near_p1_spec{iC}.subsampled_ppc(p,x1:end);
            p2_ppc = cell_res.near_p2_spec{iC}.subsampled_ppc(p,x1:end);
            
            sts_control_corrs(p) = corr(p1_sts', p2_sts');
            ppc_control_corrs(p) = corr(p1_ppc', p2_ppc');
        end
    else % MSN
        lfr_sts = cell_res.near_lfr_spec{iC}.sts_vals(x1:end);
        hfr_sts = cell_res.near_hfr_spec{iC}.sts_vals(x1:end);
        lfr_ppc = cell_res.near_lfr_spec{iC}.ppc(x1:end)';
        hfr_ppc = cell_res.near_hfr_spec{iC}.ppc(x1:end)';
        
        % Control partition correlations for STS
        sts_control_corrs = zeros(100, 1);
        ppc_control_corrs = zeros(100, 1);
        for p = 1:100
            p1_sts = cell_res.near_p1_spec{iC}.sts(p,x1:end);
            p2_sts = cell_res.near_p2_spec{iC}.sts(p,x1:end);
            p1_ppc = cell_res.near_p1_spec{iC}.ppc(p,x1:end);
            p2_ppc = cell_res.near_p2_spec{iC}.ppc(p,x1:end);
            
            sts_control_corrs(p) = corr(p1_sts', p2_sts');
            ppc_control_corrs(p) = corr(p1_ppc', p2_ppc');
        end
    end
    
    % Calculate correlations
    corr_stats.sts_lfr_hfr_corr = corr(lfr_sts', hfr_sts');
    corr_stats.ppc_lfr_hfr_corr = corr(lfr_ppc', hfr_ppc');
    corr_stats.sts_control_corr_mean = mean(sts_control_corrs);
    corr_stats.sts_control_corr_std = std(sts_control_corrs);
    corr_stats.ppc_control_corr_mean = mean(ppc_control_corrs);
    corr_stats.ppc_control_corr_std = std(ppc_control_corrs);
    corr_stats.sts_control_corrs = sts_control_corrs;
    corr_stats.ppc_control_corrs = ppc_control_corrs;
end

% Helper function: Plot STA
function plot_sta(cell_res, iC, c1, c2, c3)
    ax1 = subplot(1,3,1);
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

% Helper function: Plot STS with correlation
function plot_sts_with_correlation(cell_res, iC, c1, c2, c3, x1, this_freqs, cell_type, corr_stats)
    ax2 = subplot(1,3,2);
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
    
    % Add correlation text
    text(0.05, 0.95, sprintf('LFR-HFR corr: %.3f', corr_stats.sts_lfr_hfr_corr), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
    % Add control split correlation text below
    text(0.05, 0.85, sprintf('Control: %.3f ± %.3f', corr_stats.sts_control_corr_mean, corr_stats.sts_control_corr_std), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
    
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
    legend({'LFR','HFR','All Trials'}, 'Location','southeast');
end

% Helper function: Plot PPC with correlation
function plot_ppc_with_correlation(cell_res, iC, c1, c2, c3, x1, this_freqs, firing_rate_stats, cell_type, corr_stats)
    ax3 = subplot(1,3,3);
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
    
    % Add correlation text
    text(0.05, 0.95, sprintf('LFR-HFR corr: %.3f', corr_stats.ppc_lfr_hfr_corr), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
    % Add control split correlation text below
    text(0.05, 0.85, sprintf('Control: %.3f ± %.3f', corr_stats.ppc_control_corr_mean, corr_stats.ppc_control_corr_std), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
    
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
        'All Trials'}, 'Location','southeast');
end
