%% Main Summary Plot Script
% Load summary files
odir = 'D:\vStr_oscillatory_switch_results\cellFingerPrint\';
load(strcat(odir,'msn_summary.mat'));
load(strcat(odir,'fsi_summary.mat'));

% Define colors
c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
c2 = [26/255 255/255 26/255]; % Green for HFR
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange

% Plotting parameters (following same heuristics as plot_correlation_summary.m)
circle_marker_size = 200;
symbol_marker_size = 30;
marker_alpha = 0.4;
box_plot_width = 0.75;

% Define masks for cells with significant PPC peaks in both HFR and LFR
only_hfr_msn = ~isnan(msn_summary.hfr_ppc_peak);
only_lfr_msn = ~isnan(msn_summary.lfr_ppc_peak);
clean_msn = only_hfr_msn & only_lfr_msn;

only_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_peak);
only_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_peak);
clean_fsi = only_hfr_fsi & only_lfr_fsi;

% Define masks for significant PPC differences
dif_hfr_msn = ~isnan(msn_summary.hfr_ppc_diff);
dif_lfr_msn = ~isnan(msn_summary.lfr_ppc_diff);
dif_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_diff);
dif_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_diff);

%% ===== DATA PREPARATION =====
% Calculate delta FR and delta freq for all eligible cells
msn_delta_fr = msn_summary.hfr_mean(clean_msn) - msn_summary.lfr_mean(clean_msn);
msn_delta_freq = msn_summary.hfr_ppc_peak(clean_msn) - msn_summary.lfr_ppc_peak(clean_msn);
fsi_delta_fr = fsi_summary.hfr_mean(clean_fsi) - fsi_summary.lfr_mean(clean_fsi);
fsi_delta_freq = fsi_summary.hfr_ppc_peak(clean_fsi) - fsi_summary.lfr_ppc_peak(clean_fsi);

%% ===== CORRELATION CALCULATIONS =====

% Calculate correlations for MSNs
if sum(clean_msn) > 1
    [msn_r, msn_p] = corrcoef(msn_delta_fr, msn_delta_freq);
    msn_r_squared = msn_r(1,2)^2;
    msn_corr_p = msn_p(1,2);
else
    msn_r_squared = NaN;
    msn_corr_p = NaN;
end

% Calculate correlations for FSIs
if sum(clean_fsi) > 1
    [fsi_r, fsi_p] = corrcoef(fsi_delta_fr, fsi_delta_freq);
    fsi_r_squared = fsi_r(1,2)^2;
    fsi_corr_p = fsi_p(1,2);
else
    fsi_r_squared = NaN;
    fsi_corr_p = NaN;
end

% Calculate correlations for each condition separately
% HFR diff only
only_hfr_diff_msn = clean_msn & dif_hfr_msn & ~dif_lfr_msn;
only_hfr_diff_fsi = clean_fsi & dif_hfr_fsi & ~dif_lfr_fsi;

if sum(only_hfr_diff_msn) > 1
    msn_hfr_delta_fr = msn_summary.hfr_mean(only_hfr_diff_msn) - msn_summary.lfr_mean(only_hfr_diff_msn);
    msn_hfr_delta_freq = msn_summary.hfr_ppc_peak(only_hfr_diff_msn) - msn_summary.lfr_ppc_peak(only_hfr_diff_msn);
    [msn_hfr_r, msn_hfr_p] = corrcoef(msn_hfr_delta_fr, msn_hfr_delta_freq);
    msn_hfr_r_squared = msn_hfr_r(1,2)^2;
    msn_hfr_corr_p = msn_hfr_p(1,2);
else
    msn_hfr_r_squared = NaN;
    msn_hfr_corr_p = NaN;
end

if sum(only_hfr_diff_fsi) > 1
    fsi_hfr_delta_fr = fsi_summary.hfr_mean(only_hfr_diff_fsi) - fsi_summary.lfr_mean(only_hfr_diff_fsi);
    fsi_hfr_delta_freq = fsi_summary.hfr_ppc_peak(only_hfr_diff_fsi) - fsi_summary.lfr_ppc_peak(only_hfr_diff_fsi);
    [fsi_hfr_r, fsi_hfr_p] = corrcoef(fsi_hfr_delta_fr, fsi_hfr_delta_freq);
    fsi_hfr_r_squared = fsi_hfr_r(1,2)^2;
    fsi_hfr_corr_p = fsi_hfr_p(1,2);
else
    fsi_hfr_r_squared = NaN;
    fsi_hfr_corr_p = NaN;
end

% LFR diff only
only_lfr_diff_msn = clean_msn & ~dif_hfr_msn & dif_lfr_msn;
only_lfr_diff_fsi = clean_fsi & ~dif_hfr_fsi & dif_lfr_fsi;

if sum(only_lfr_diff_msn) > 1
    msn_lfr_delta_fr = msn_summary.hfr_mean(only_lfr_diff_msn) - msn_summary.lfr_mean(only_lfr_diff_msn);
    msn_lfr_delta_freq = msn_summary.hfr_ppc_peak(only_lfr_diff_msn) - msn_summary.lfr_ppc_peak(only_lfr_diff_msn);
    [msn_lfr_r, msn_lfr_p] = corrcoef(msn_lfr_delta_fr, msn_lfr_delta_freq);
    msn_lfr_r_squared = msn_lfr_r(1,2)^2;
    msn_lfr_corr_p = msn_lfr_p(1,2);
else
    msn_lfr_r_squared = NaN;
    msn_lfr_corr_p = NaN;
end

if sum(only_lfr_diff_fsi) > 1
    fsi_lfr_delta_fr = fsi_summary.hfr_mean(only_lfr_diff_fsi) - fsi_summary.lfr_mean(only_lfr_diff_fsi);
    fsi_lfr_delta_freq = fsi_summary.hfr_ppc_peak(only_lfr_diff_fsi) - fsi_summary.lfr_ppc_peak(only_lfr_diff_fsi);
    [fsi_lfr_r, fsi_lfr_p] = corrcoef(fsi_lfr_delta_fr, fsi_lfr_delta_freq);
    fsi_lfr_r_squared = fsi_lfr_r(1,2)^2;
    fsi_lfr_corr_p = fsi_lfr_p(1,2);
else
    fsi_lfr_r_squared = NaN;
    fsi_lfr_corr_p = NaN;
end

% Both diffs
both_diff_msn = clean_msn & dif_hfr_msn & dif_lfr_msn;
both_diff_fsi = clean_fsi & dif_hfr_fsi & dif_lfr_fsi;

if sum(both_diff_msn) > 1
    msn_both_delta_fr = msn_summary.hfr_mean(both_diff_msn) - msn_summary.lfr_mean(both_diff_msn);
    msn_both_delta_freq = msn_summary.hfr_ppc_peak(both_diff_msn) - msn_summary.lfr_ppc_peak(both_diff_msn);
    [msn_both_r, msn_both_p] = corrcoef(msn_both_delta_fr, msn_both_delta_freq);
    msn_both_r_squared = msn_both_r(1,2)^2;
    msn_both_corr_p = msn_both_p(1,2);
else
    msn_both_r_squared = NaN;
    msn_both_corr_p = NaN;
end

if sum(both_diff_fsi) > 1
    fsi_both_delta_fr = fsi_summary.hfr_mean(both_diff_fsi) - fsi_summary.lfr_mean(both_diff_fsi);
    fsi_both_delta_freq = fsi_summary.hfr_ppc_peak(both_diff_fsi) - fsi_summary.lfr_ppc_peak(both_diff_fsi);
    [fsi_both_r, fsi_both_p] = corrcoef(fsi_both_delta_fr, fsi_both_delta_freq);
    fsi_both_r_squared = fsi_both_r(1,2)^2;
    fsi_both_corr_p = fsi_both_p(1,2);
else
    fsi_both_r_squared = NaN;
    fsi_both_corr_p = NaN;
end

%% ===== FIGURE CREATION =====

% Create figure
fig = figure('WindowState', 'maximized');

% Left Column: Main Summary Plot (All eligible cells with symbols) - spans all rows
ax1 = subplot(3,2,[1,3,5]);
hold on;

% Plot MSNs - First scatter all eligible MSNs as circles
s1 = scatter(ax1, msn_delta_fr, msn_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);

% Plot FSIs - First scatter all eligible FSIs as circles
s2 = scatter(ax1, fsi_delta_fr, fsi_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);

% Now overlay markers for significant PPC differences
% MSNs with only HFR significant PPC diff
only_hfr_diff_msn_data = only_hfr_diff_msn(clean_msn);
if sum(only_hfr_diff_msn_data) > 0
    s3 = scatter(ax1, msn_delta_fr(only_hfr_diff_msn_data), msn_delta_freq(only_hfr_diff_msn_data), symbol_marker_size,'black', 'Marker', '>');
end

% MSNs with only LFR significant PPC diff
only_lfr_diff_msn_data = only_lfr_diff_msn(clean_msn);
if sum(only_lfr_diff_msn_data) > 0
    s4 = scatter(ax1, msn_delta_fr(only_lfr_diff_msn_data), msn_delta_freq(only_lfr_diff_msn_data), symbol_marker_size,'black', 'Marker', '<');
end

% MSNs with both HFR and LFR significant PPC diff
both_diff_msn_data = both_diff_msn(clean_msn);
if sum(both_diff_msn_data) > 0
    s5 = scatter(ax1, msn_delta_fr(both_diff_msn_data), msn_delta_freq(both_diff_msn_data), symbol_marker_size,'black', 'Marker', 'x');
end

% FSIs with only HFR significant PPC diff
only_hfr_diff_fsi_data = only_hfr_diff_fsi(clean_fsi);
if sum(only_hfr_diff_fsi_data) > 0
    s6 = scatter(ax1, fsi_delta_fr(only_hfr_diff_fsi_data), fsi_delta_freq(only_hfr_diff_fsi_data), symbol_marker_size,'black', 'Marker', '>');
end

% FSIs with only LFR significant PPC diff
only_lfr_diff_fsi_data = only_lfr_diff_fsi(clean_fsi);
if sum(only_lfr_diff_fsi_data) > 0
    s7 = scatter(ax1, fsi_delta_fr(only_lfr_diff_fsi_data), fsi_delta_freq(only_lfr_diff_fsi_data), symbol_marker_size,'black', 'Marker', '<');
end

% FSIs with both HFR and LFR significant PPC diff
both_diff_fsi_data = both_diff_fsi(clean_fsi);
if sum(both_diff_fsi_data) > 0
    s8 = scatter(ax1, fsi_delta_fr(both_diff_fsi_data), fsi_delta_freq(both_diff_fsi_data), symbol_marker_size,'black', 'Marker', 'x');
end

ax1.XLabel.String = '\Delta Firing Rate';
ax1.YLabel.String = '\Delta Peak Frequency (Hz)';
ax1.FontSize = 16;
ax1.Title.String = 'All Eligible Cells';
ax1.YLim = [-100 100];

% Create legend for main plot
legend_entries = {};
legend_handles = [];

% Count total eligible cells for legend
eligible_msn = sum(clean_msn);
total_msn = length(clean_msn);
eligible_fsi = sum(clean_fsi);
total_fsi = length(clean_fsi);

% Count marker types for each cell category
only_hfr_diff_msn_count = sum(only_hfr_diff_msn);
only_lfr_diff_msn_count = sum(only_lfr_diff_msn);
both_diff_msn_count = sum(both_diff_msn);

only_hfr_diff_fsi_count = sum(only_hfr_diff_fsi);
only_lfr_diff_fsi_count = sum(only_lfr_diff_fsi);
both_diff_fsi_count = sum(both_diff_fsi);

if total_msn > 0
    if ~isnan(msn_r_squared)
        legend_entries{end+1} = sprintf('MSN: %d/%d cells (R²=%.3f, p=%.3f)', eligible_msn, total_msn, msn_r_squared, msn_corr_p);
    else
        legend_entries{end+1} = sprintf('MSN: %d/%d cells', eligible_msn, total_msn);
    end
    legend_handles(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

if total_fsi > 0
    if ~isnan(fsi_r_squared)
        legend_entries{end+1} = sprintf('FSI: %d/%d cells (R²=%.3f, p=%.3f)', eligible_fsi, total_fsi, fsi_r_squared, fsi_corr_p);
    else
        legend_entries{end+1} = sprintf('FSI: %d/%d cells', eligible_fsi, total_fsi);
    end
    legend_handles(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

% Add marker legend entries with counts
if only_hfr_diff_msn_count > 0 || only_hfr_diff_fsi_count > 0
    legend_entries{end+1} = sprintf('> = HFR PPC diff only (MSN:%d, FSI:%d)', only_hfr_diff_msn_count, only_hfr_diff_fsi_count);
    legend_handles(end+1) = scatter(NaN, NaN, symbol_marker_size, 'black', 'Marker', '>');
end

if only_lfr_diff_msn_count > 0 || only_lfr_diff_fsi_count > 0
    legend_entries{end+1} = sprintf('< = LFR PPC diff only (MSN:%d, FSI:%d)', only_lfr_diff_msn_count, only_lfr_diff_fsi_count);
    legend_handles(end+1) = scatter(NaN, NaN, symbol_marker_size, 'black', 'Marker', '<');
end

if both_diff_msn_count > 0 || both_diff_fsi_count > 0
    legend_entries{end+1} = sprintf('x = Both PPC diffs (MSN:%d, FSI:%d)', both_diff_msn_count, both_diff_fsi_count);
    legend_handles(end+1) = scatter(NaN, NaN, symbol_marker_size, 'black', 'Marker', 'x');
end

if ~isempty(legend_entries)
    legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 12);
end

% Right Column, Row 1: Only HFR significant PPC diff
ax2 = subplot(3,2,2);
hold on;

% MSNs with only HFR significant PPC diff
if sum(only_hfr_diff_msn) > 0
    scatter(ax2, msn_hfr_delta_fr, msn_hfr_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

% FSIs with only HFR significant PPC diff
if sum(only_hfr_diff_fsi) > 0
    scatter(ax2, fsi_hfr_delta_fr, fsi_hfr_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

ax2.XLabel.String = '\Delta Firing Rate';
ax2.YLabel.String = '\Delta Peak Frequency (Hz)';
ax2.FontSize = 16;
ax2.Title.String = sprintf('HFR PPC Diff Only (MSN:%d, FSI:%d)', only_hfr_diff_msn_count, only_hfr_diff_fsi_count);
ax2.YLim = [-100 100];

% Simple legend for this subplot
if sum(only_hfr_diff_msn) > 0 || sum(only_hfr_diff_fsi) > 0
    legend_handles_hfr = [];
    legend_entries_hfr = {};
    
    if sum(only_hfr_diff_msn) > 0
        if ~isnan(msn_hfr_r_squared)
            legend_entries_hfr{end+1} = sprintf('MSN (R²=%.3f, p=%.3f)', msn_hfr_r_squared, msn_hfr_corr_p);
        else
            legend_entries_hfr{end+1} = 'MSN';
        end
        legend_handles_hfr(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    end
    
    if sum(only_hfr_diff_fsi) > 0
        if ~isnan(fsi_hfr_r_squared)
            legend_entries_hfr{end+1} = sprintf('FSI (R²=%.3f, p=%.3f)', fsi_hfr_r_squared, fsi_hfr_corr_p);
        else
            legend_entries_hfr{end+1} = 'FSI';
        end
        legend_handles_hfr(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    end
    
    legend(legend_handles_hfr, legend_entries_hfr, 'Location', 'best', 'FontSize', 12);
end

% Right Column, Row 2: Only LFR significant PPC diff
ax3 = subplot(3,2,4);
hold on;

% MSNs with only LFR significant PPC diff
if sum(only_lfr_diff_msn) > 0
    scatter(ax3, msn_lfr_delta_fr, msn_lfr_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

% FSIs with only LFR significant PPC diff
if sum(only_lfr_diff_fsi) > 0
    scatter(ax3, fsi_lfr_delta_fr, fsi_lfr_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

ax3.XLabel.String = '\Delta Firing Rate';
ax3.YLabel.String = '\Delta Peak Frequency (Hz)';
ax3.FontSize = 16;
ax3.Title.String = sprintf('LFR PPC Diff Only (MSN:%d, FSI:%d)', only_lfr_diff_msn_count, only_lfr_diff_fsi_count);
ax3.YLim = [-100 100];

% Simple legend for this subplot
if sum(only_lfr_diff_msn) > 0 || sum(only_lfr_diff_fsi) > 0
    legend_handles_lfr = [];
    legend_entries_lfr = {};
    
    if sum(only_lfr_diff_msn) > 0
        if ~isnan(msn_lfr_r_squared)
            legend_entries_lfr{end+1} = sprintf('MSN (R²=%.3f, p=%.3f)', msn_lfr_r_squared, msn_lfr_corr_p);
        else
            legend_entries_lfr{end+1} = 'MSN';
        end
        legend_handles_lfr(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    end
    
    if sum(only_lfr_diff_fsi) > 0
        if ~isnan(fsi_lfr_r_squared)
            legend_entries_lfr{end+1} = sprintf('FSI (R²=%.3f, p=%.3f)', fsi_lfr_r_squared, fsi_lfr_corr_p);
        else
            legend_entries_lfr{end+1} = 'FSI';
        end
        legend_handles_lfr(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    end
    
    legend(legend_handles_lfr, legend_entries_lfr, 'Location', 'best', 'FontSize', 12);
end

% Right Column, Row 3: Both significant PPC diffs
ax4 = subplot(3,2,6);
hold on;

% MSNs with both HFR and LFR significant PPC diff
if sum(both_diff_msn) > 0
    scatter(ax4, msn_both_delta_fr, msn_both_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

% FSIs with both HFR and LFR significant PPC diff
if sum(both_diff_fsi) > 0
    scatter(ax4, fsi_both_delta_fr, fsi_both_delta_freq, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

ax4.XLabel.String = '\Delta Firing Rate';
ax4.YLabel.String = '\Delta Peak Frequency (Hz)';
ax4.FontSize = 16;
ax4.Title.String = sprintf('Both PPC Diffs (MSN:%d, FSI:%d)', both_diff_msn_count, both_diff_fsi_count);
ax4.YLim = [-100 100];

% Simple legend for this subplot
if sum(both_diff_msn) > 0 || sum(both_diff_fsi) > 0
    legend_handles_both = [];
    legend_entries_both = {};
    
    if sum(both_diff_msn) > 0
        if ~isnan(msn_both_r_squared)
            legend_entries_both{end+1} = sprintf('MSN (R²=%.3f, p=%.3f)', msn_both_r_squared, msn_both_corr_p);
        else
            legend_entries_both{end+1} = 'MSN';
        end
        legend_handles_both(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', c4, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    end
    
    if sum(both_diff_fsi) > 0
        if ~isnan(fsi_both_r_squared)
            legend_entries_both{end+1} = sprintf('FSI (R²=%.3f, p=%.3f)', fsi_both_r_squared, fsi_both_corr_p);
        else
            legend_entries_both{end+1} = 'FSI';
        end
        legend_handles_both(end+1) = scatter(NaN, NaN, 'SizeData', circle_marker_size, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    end
    
    legend(legend_handles_both, legend_entries_both, 'Location', 'best', 'FontSize', 12);
end

% Adjust subplot spacing
% sgtitle('Main Summary: Delta Firing Rate vs Delta Peak Frequency', 'FontSize', 18, 'FontWeight', 'bold');

%% Create second figure with histogram distributions
fig2 = figure('WindowState', 'maximized');

hist_bins = -105:10:105;
% Top Left: All Eligible Cells - Delta Frequency Distribution
ax5 = subplot(2,2,1);
hold on;

% Calculate delta frequencies for all eligible cells
all_msn_delta_freq = msn_summary.hfr_ppc_peak(clean_msn) - msn_summary.lfr_ppc_peak(clean_msn);
all_fsi_delta_freq = fsi_summary.hfr_ppc_peak(clean_fsi) - fsi_summary.lfr_ppc_peak(clean_fsi);

% Create histograms with transparency
if ~isempty(all_msn_delta_freq)
    histogram(ax5, all_msn_delta_freq, hist_bins, 'FaceColor', c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
end
if ~isempty(all_fsi_delta_freq)
    histogram(ax5, all_fsi_delta_freq, hist_bins, 'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
end

% Add vertical lines for means
if ~isempty(all_msn_delta_freq)
    msn_mean_delta_freq = mean(all_msn_delta_freq);
    xline(ax5, msn_mean_delta_freq, 'Color', c4, 'LineWidth', 2, 'LineStyle', '--');
end
if ~isempty(all_fsi_delta_freq)
    fsi_mean_delta_freq = mean(all_fsi_delta_freq);
    xline(ax5, fsi_mean_delta_freq, 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--');
end

ax5.XLabel.String = '\Delta Peak Frequency (Hz)';
ax5.YLabel.String = 'Proportion';
ax5.FontSize = 16;
ax5.Title.String = 'All Eligible Cells';
ax5.XLim = [-100 100];

% Legend
legend_handles_hist1 = [];
legend_entries_hist1 = {};
if ~isempty(all_msn_delta_freq)
    legend_entries_hist1{end+1} = sprintf('MSN (n=%d)', length(all_msn_delta_freq));
    legend_handles_hist1(end+1) = patch(NaN, NaN, c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if ~isempty(all_fsi_delta_freq)
    legend_entries_hist1{end+1} = sprintf('FSI (n=%d)', length(all_fsi_delta_freq));
    legend_handles_hist1(end+1) = patch(NaN, NaN, 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if ~isempty(legend_entries_hist1)
    legend(legend_handles_hist1, legend_entries_hist1, 'Location', 'best', 'FontSize', 12);
end

% Top Right: Only HFR significant PPC diff - Delta Frequency Distribution
ax6 = subplot(2,2,2);
hold on;

% Calculate delta frequencies for HFR-only cells
if sum(only_hfr_diff_msn) > 0
    msn_hfr_only_delta_freq = msn_summary.hfr_ppc_peak(only_hfr_diff_msn) - msn_summary.lfr_ppc_peak(only_hfr_diff_msn);
    histogram(ax6, msn_hfr_only_delta_freq, hist_bins, 'FaceColor', c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
    xline(ax6, mean(msn_hfr_only_delta_freq), 'Color', c4, 'LineWidth', 2, 'LineStyle', '--');
end

if sum(only_hfr_diff_fsi) > 0
    fsi_hfr_only_delta_freq = fsi_summary.hfr_ppc_peak(only_hfr_diff_fsi) - fsi_summary.lfr_ppc_peak(only_hfr_diff_fsi);
    histogram(ax6, fsi_hfr_only_delta_freq, hist_bins, 'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
    xline(ax6, mean(fsi_hfr_only_delta_freq), 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--');
end

ax6.XLabel.String = '\Delta Peak Frequency (Hz)';
ax6.YLabel.String = 'Proportion';
ax6.FontSize = 16;
ax6.Title.String = sprintf('HFR PPC Diff Only (MSN:%d, FSI:%d)', only_hfr_diff_msn_count, only_hfr_diff_fsi_count);
ax6.XLim = [-100 100];

% Legend
legend_handles_hist2 = [];
legend_entries_hist2 = {};
if sum(only_hfr_diff_msn) > 0
    legend_entries_hist2{end+1} = sprintf('MSN (n=%d)', only_hfr_diff_msn_count);
    legend_handles_hist2(end+1) = patch(NaN, NaN, c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if sum(only_hfr_diff_fsi) > 0
    legend_entries_hist2{end+1} = sprintf('FSI (n=%d)', only_hfr_diff_fsi_count);
    legend_handles_hist2(end+1) = patch(NaN, NaN, 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if ~isempty(legend_entries_hist2)
    legend(legend_handles_hist2, legend_entries_hist2, 'Location', 'best', 'FontSize', 12);
end

% Bottom Left: Only LFR significant PPC diff - Delta Frequency Distribution
ax7 = subplot(2,2,3);
hold on;

% Calculate delta frequencies for LFR-only cells
if sum(only_lfr_diff_msn) > 0
    msn_lfr_only_delta_freq = msn_summary.hfr_ppc_peak(only_lfr_diff_msn) - msn_summary.lfr_ppc_peak(only_lfr_diff_msn);
    histogram(ax7, msn_lfr_only_delta_freq, hist_bins, 'FaceColor', c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
    xline(ax7, mean(msn_lfr_only_delta_freq), 'Color', c4, 'LineWidth', 2, 'LineStyle', '--');
end

if sum(only_lfr_diff_fsi) > 0
    fsi_lfr_only_delta_freq = fsi_summary.hfr_ppc_peak(only_lfr_diff_fsi) - fsi_summary.lfr_ppc_peak(only_lfr_diff_fsi);
    histogram(ax7, fsi_lfr_only_delta_freq, hist_bins, 'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
    xline(ax7, mean(fsi_lfr_only_delta_freq), 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--');
end

ax7.XLabel.String = '\Delta Peak Frequency (Hz)';
ax7.YLabel.String = 'Proportion';
ax7.FontSize = 16;
ax7.Title.String = sprintf('LFR PPC Diff Only (MSN:%d, FSI:%d)', only_lfr_diff_msn_count, only_lfr_diff_fsi_count);
ax7.XLim = [-100 100];

% Legend
legend_handles_hist3 = [];
legend_entries_hist3 = {};
if sum(only_lfr_diff_msn) > 0
    legend_entries_hist3{end+1} = sprintf('MSN (n=%d)', only_lfr_diff_msn_count);
    legend_handles_hist3(end+1) = patch(NaN, NaN, c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if sum(only_lfr_diff_fsi) > 0
    legend_entries_hist3{end+1} = sprintf('FSI (n=%d)', only_lfr_diff_fsi_count);
    legend_handles_hist3(end+1) = patch(NaN, NaN, 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if ~isempty(legend_entries_hist3)
    legend(legend_handles_hist3, legend_entries_hist3, 'Location', 'best', 'FontSize', 12);
end

% Bottom Right: Both significant PPC diffs - Delta Frequency Distribution
ax8 = subplot(2,2,4);
hold on;

% Calculate delta frequencies for both-diff cells
if sum(both_diff_msn) > 0
    msn_both_delta_freq = msn_summary.hfr_ppc_peak(both_diff_msn) - msn_summary.lfr_ppc_peak(both_diff_msn);
    histogram(ax8, msn_both_delta_freq, hist_bins, 'FaceColor', c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
    xline(ax8, mean(msn_both_delta_freq), 'Color', c4, 'LineWidth', 2, 'LineStyle', '--');
end

if sum(both_diff_fsi) > 0
    fsi_both_delta_freq = fsi_summary.hfr_ppc_peak(both_diff_fsi) - fsi_summary.lfr_ppc_peak(both_diff_fsi);
    histogram(ax8, fsi_both_delta_freq, hist_bins, 'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'Normalization', 'probability');
    xline(ax8, mean(fsi_both_delta_freq), 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--');
end

ax8.XLabel.String = '\Delta Peak Frequency (Hz)';
ax8.YLabel.String = 'Proportion';
ax8.FontSize = 16;
ax8.Title.String = sprintf('Both PPC Diffs (MSN:%d, FSI:%d)', both_diff_msn_count, both_diff_fsi_count);
ax8.XLim = [-100 100];

% Legend
legend_handles_hist4 = [];
legend_entries_hist4 = {};
if sum(both_diff_msn) > 0
    legend_entries_hist4{end+1} = sprintf('MSN (n=%d)', both_diff_msn_count);
    legend_handles_hist4(end+1) = patch(NaN, NaN, c4, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if sum(both_diff_fsi) > 0
    legend_entries_hist4{end+1} = sprintf('FSI (n=%d)', both_diff_fsi_count);
    legend_handles_hist4(end+1) = patch(NaN, NaN, 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end
if ~isempty(legend_entries_hist4)
    legend(legend_handles_hist4, legend_entries_hist4, 'Location', 'best', 'FontSize', 12);
end

% Print summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Total MSNs: %d, Eligible MSNs: %d\n', total_msn, eligible_msn);
fprintf('Total FSIs: %d, Eligible FSIs: %d\n', total_fsi, eligible_fsi);
fprintf('\nMSN PPC Differences:\n');
fprintf('  HFR only: %d\n', only_hfr_diff_msn_count);
fprintf('  LFR only: %d\n', only_lfr_diff_msn_count);
fprintf('  Both: %d\n', both_diff_msn_count);
fprintf('\nFSI PPC Differences:\n');
fprintf('  HFR only: %d\n', only_hfr_diff_fsi_count);
fprintf('  LFR only: %d\n', only_lfr_diff_fsi_count);
fprintf('  Both: %d\n', both_diff_fsi_count);
