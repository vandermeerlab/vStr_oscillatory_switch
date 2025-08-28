%% Script to plot results from diagnoseSTAcalc.m
% This script loads the saved results and creates the comparison plots

clear;
close all;

% PARAMETERS
results_dir = 'D:\vStr_oscillatory_switch_results\temp4'; % Directory with saved results
cell_name = 'R132-2007-10-21-TT12_4.t'; % Cell to plot (change this)
session_name = 'R132-2007-10-21'; % Session name (extracted from cell name)

% Load results
results_file = fullfile(results_dir, sprintf('%s_%s_results.mat', session_name, cell_name));
if ~exist(results_file, 'file')
    error('Results file not found: %s', results_file);
end

load(results_file);
fprintf('Loaded results for cell: %s\n', cell_name);

% Create plots
plotCellComparison(cfg_in, cell_name, existing_results, new_results);

%%
% Function to plot cell comparison (same as in diagnoseSTAcalc.m)
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
     
     % Add cell label as title
     sgtitle(sprintf('%s - Method: %s (No subsampling)', cell_label, cfg_in.ppc_method), ...
         'FontSize', 16, 'Interpreter', 'None');
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
