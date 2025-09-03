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

circle_marker_size = 200;
symbol_marker_size = 30;

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

% Create figure
fig = figure('WindowState', 'maximized');

%% Top Left: Main Summary Plot (All eligible cells with symbols)
ax1 = subplot(2,2,1);
hold on;

% Plot MSNs - First scatter all eligible MSNs as circles
msn_delta_fr = msn_summary.hfr_mean(clean_msn) - msn_summary.lfr_mean(clean_msn);
msn_delta_freq = msn_summary.hfr_ppc_peak(clean_msn) - msn_summary.lfr_ppc_peak(clean_msn);
s1 = scatter(ax1, msn_delta_fr, msn_delta_freq, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.4);

% Plot FSIs - First scatter all eligible FSIs as circles
fsi_delta_fr = fsi_summary.hfr_mean(clean_fsi) - fsi_summary.lfr_mean(clean_fsi);
fsi_delta_freq = fsi_summary.hfr_ppc_peak(clean_fsi) - fsi_summary.lfr_ppc_peak(clean_fsi);
s2 = scatter(ax1, fsi_delta_fr, fsi_delta_freq, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.4);

% Now overlay markers for significant PPC differences
% MSNs with only HFR significant PPC diff
only_hfr_diff_msn = clean_msn & dif_hfr_msn & ~dif_lfr_msn;
only_hfr_diff_msn_data = only_hfr_diff_msn(clean_msn);
if sum(only_hfr_diff_msn_data) > 0
    s3 = scatter(ax1, msn_delta_fr(only_hfr_diff_msn_data), msn_delta_freq(only_hfr_diff_msn_data), symbol_marker_size,'black', 'Marker', '>');
end

% MSNs with only LFR significant PPC diff
only_lfr_diff_msn = clean_msn & ~dif_hfr_msn & dif_lfr_msn;
only_lfr_diff_msn_data = only_lfr_diff_msn(clean_msn);
if sum(only_lfr_diff_msn_data) > 0
    s4 = scatter(ax1, msn_delta_fr(only_lfr_diff_msn_data), msn_delta_freq(only_lfr_diff_msn_data), symbol_marker_size,'black', 'Marker', '<');
end

% MSNs with both HFR and LFR significant PPC diff
both_diff_msn = clean_msn & dif_hfr_msn & dif_lfr_msn;
both_diff_msn_data = both_diff_msn(clean_msn);
if sum(both_diff_msn_data) > 0
    s5 = scatter(ax1, msn_delta_fr(both_diff_msn_data), msn_delta_freq(both_diff_msn_data), symbol_marker_size,'black', 'Marker', 'x');
end

% FSIs with only HFR significant PPC diff
only_hfr_diff_fsi = clean_fsi & dif_hfr_fsi & ~dif_lfr_fsi;
only_hfr_diff_fsi_data = only_hfr_diff_fsi(clean_fsi);
if sum(only_hfr_diff_fsi_data) > 0
    s6 = scatter(ax1, fsi_delta_fr(only_hfr_diff_fsi_data), fsi_delta_freq(only_hfr_diff_fsi_data), symbol_marker_size,'black', 'Marker', '>');
end

% FSIs with only LFR significant PPC diff
only_lfr_diff_fsi = clean_fsi & ~dif_hfr_fsi & dif_lfr_fsi;
only_lfr_diff_fsi_data = only_lfr_diff_fsi(clean_fsi);
if sum(only_lfr_diff_fsi_data) > 0
    s7 = scatter(ax1, fsi_delta_fr(only_lfr_diff_fsi_data), fsi_delta_freq(only_lfr_diff_fsi_data), symbol_marker_size,'black', 'Marker', '<');
end

% FSIs with both HFR and LFR significant PPC diff
both_diff_fsi = clean_fsi & dif_hfr_fsi & dif_lfr_fsi;
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
    legend_entries{end+1} = sprintf('MSN: %d/%d cells', eligible_msn, total_msn);
    legend_handles(end+1) = scatter(NaN, NaN, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.4);
end

if total_fsi > 0
    legend_entries{end+1} = sprintf('FSI: %d/%d cells', eligible_fsi, total_fsi);
    legend_handles(end+1) = scatter(NaN, NaN, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.4);
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

%% Top Right: Only HFR significant PPC diff
ax2 = subplot(2,2,2);
hold on;

% MSNs with only HFR significant PPC diff
if sum(only_hfr_diff_msn) > 0
    msn_hfr_only_delta_fr = msn_summary.hfr_mean(only_hfr_diff_msn) - msn_summary.lfr_mean(only_hfr_diff_msn);
    msn_hfr_only_delta_freq = msn_summary.hfr_ppc_peak(only_hfr_diff_msn) - msn_summary.lfr_ppc_peak(only_hfr_diff_msn);
    scatter(ax2, msn_hfr_only_delta_fr, msn_hfr_only_delta_freq, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.6);
end

% FSIs with only HFR significant PPC diff
if sum(only_hfr_diff_fsi) > 0
    fsi_hfr_only_delta_fr = fsi_summary.hfr_mean(only_hfr_diff_fsi) - fsi_summary.lfr_mean(only_hfr_diff_fsi);
    fsi_hfr_only_delta_freq = fsi_summary.hfr_ppc_peak(only_hfr_diff_fsi) - fsi_summary.lfr_ppc_peak(only_hfr_diff_fsi);
    scatter(ax2, fsi_hfr_only_delta_fr, fsi_hfr_only_delta_freq, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
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
        legend_entries_hfr{end+1} = 'MSN';
        legend_handles_hfr(end+1) = scatter(NaN, NaN, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    if sum(only_hfr_diff_fsi) > 0
        legend_entries_hfr{end+1} = 'FSI';
        legend_handles_hfr(end+1) = scatter(NaN, NaN, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    legend(legend_handles_hfr, legend_entries_hfr, 'Location', 'best', 'FontSize', 12);
end

%% Bottom Left: Only LFR significant PPC diff
ax3 = subplot(2,2,3);
hold on;

% MSNs with only LFR significant PPC diff
if sum(only_lfr_diff_msn) > 0
    msn_lfr_only_delta_fr = msn_summary.hfr_mean(only_lfr_diff_msn) - msn_summary.lfr_mean(only_lfr_diff_msn);
    msn_lfr_only_delta_freq = msn_summary.hfr_ppc_peak(only_lfr_diff_msn) - msn_summary.lfr_ppc_peak(only_lfr_diff_msn);
    scatter(ax3, msn_lfr_only_delta_fr, msn_lfr_only_delta_freq, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.6);
end

% FSIs with only LFR significant PPC diff
if sum(only_lfr_diff_fsi) > 0
    fsi_lfr_only_delta_fr = fsi_summary.hfr_mean(only_lfr_diff_fsi) - fsi_summary.lfr_mean(only_lfr_diff_fsi);
    fsi_lfr_only_delta_freq = fsi_summary.hfr_ppc_peak(only_lfr_diff_fsi) - fsi_summary.lfr_ppc_peak(only_lfr_diff_fsi);
    scatter(ax3, fsi_lfr_only_delta_fr, fsi_lfr_only_delta_freq, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
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
        legend_entries_lfr{end+1} = 'MSN';
        legend_handles_lfr(end+1) = scatter(NaN, NaN, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    if sum(only_lfr_diff_fsi) > 0
        legend_entries_lfr{end+1} = 'FSI';
        legend_handles_lfr(end+1) = scatter(NaN, NaN, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    legend(legend_handles_lfr, legend_entries_lfr, 'Location', 'best', 'FontSize', 12);
end

%% Bottom Right: Both significant PPC diffs
ax4 = subplot(2,2,4);
hold on;

% MSNs with both HFR and LFR significant PPC diff
if sum(both_diff_msn) > 0
    msn_both_delta_fr = msn_summary.hfr_mean(both_diff_msn) - msn_summary.lfr_mean(both_diff_msn);
    msn_both_delta_freq = msn_summary.hfr_ppc_peak(both_diff_msn) - msn_summary.lfr_ppc_peak(both_diff_msn);
    scatter(ax4, msn_both_delta_fr, msn_both_delta_freq, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.6);
end

% FSIs with both HFR and LFR significant PPC diff
if sum(both_diff_fsi) > 0
    fsi_both_delta_fr = fsi_summary.hfr_mean(both_diff_fsi) - fsi_summary.lfr_mean(both_diff_fsi);
    fsi_both_delta_freq = fsi_summary.hfr_ppc_peak(both_diff_fsi) - fsi_summary.lfr_ppc_peak(both_diff_fsi);
    scatter(ax4, fsi_both_delta_fr, fsi_both_delta_freq, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
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
        legend_entries_both{end+1} = 'MSN';
        legend_handles_both(end+1) = scatter(NaN, NaN, circle_marker_size, c4, 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    if sum(both_diff_fsi) > 0
        legend_entries_both{end+1} = 'FSI';
        legend_handles_both(end+1) = scatter(NaN, NaN, circle_marker_size, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    legend(legend_handles_both, legend_entries_both, 'Location', 'best', 'FontSize', 12);
end

% Adjust subplot spacing
% sgtitle('Main Summary: Delta Firing Rate vs Delta Peak Frequency', 'FontSize', 18, 'FontWeight', 'bold');

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
