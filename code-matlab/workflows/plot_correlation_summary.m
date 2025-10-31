%% Load summary files
load('D:\vStr_oscillatory_switch_results\cellFingerPrintsWithCorrelation\msn_summary.mat');
load('D:\vStr_oscillatory_switch_results\cellFingerPrintsWithCorrelation\fsi_summary.mat');

% Firing rate difference threshold (Hz)
fr_diff_threshold = 0;

% Define firing rate difference mask
fr_diff_msn = abs(msn_summary.hfr_mean - msn_summary.lfr_mean) >= fr_diff_threshold;
fr_diff_fsi = abs(fsi_summary.hfr_mean - fsi_summary.lfr_mean) >= fr_diff_threshold;

% Extract correlation data
% MSN data
msn_ppc_corr = msn_summary.ppc_lfr_hfr_corr(fr_diff_msn);
msn_sts_corr = msn_summary.sts_lfr_hfr_corr(fr_diff_msn);
msn_ppc_control_corrs = msn_summary.ppc_control_corrs(fr_diff_msn);
msn_sts_control_corrs = msn_summary.sts_control_corrs(fr_diff_msn);

% FSI data
fsi_ppc_corr = fsi_summary.ppc_lfr_hfr_corr(fr_diff_fsi);
fsi_sts_corr = fsi_summary.sts_lfr_hfr_corr(fr_diff_fsi);
fsi_ppc_control_corrs = fsi_summary.ppc_control_corrs(fr_diff_fsi);
fsi_sts_control_corrs = fsi_summary.sts_control_corrs(fr_diff_fsi);

% Calculate control split statistics
% Calculate std and mean of control correlations for each cell
msn_ppc_control_std = cellfun(@(x) std(x), msn_ppc_control_corrs);
msn_ppc_control_mean = cellfun(@(x) mean(x), msn_ppc_control_corrs);
msn_sts_control_std = cellfun(@(x) std(x), msn_sts_control_corrs);
msn_sts_control_mean = cellfun(@(x) mean(x), msn_sts_control_corrs);

fsi_ppc_control_std = cellfun(@(x) std(x), fsi_ppc_control_corrs);
fsi_ppc_control_mean = cellfun(@(x) mean(x), fsi_ppc_control_corrs);
fsi_sts_control_std = cellfun(@(x) std(x), fsi_sts_control_corrs);
fsi_sts_control_mean = cellfun(@(x) mean(x), fsi_sts_control_corrs);

% Calculate percentiles
% Initialize arrays for percentiles
msn_ppc_percentiles = zeros(length(msn_ppc_corr), 1);
msn_sts_percentiles = zeros(length(msn_sts_corr), 1);
fsi_ppc_percentiles = zeros(length(fsi_ppc_corr), 1);
fsi_sts_percentiles = zeros(length(fsi_sts_corr), 1);

% Calculate percentiles for MSNs
for i = 1:length(msn_ppc_corr)
    % PPC percentile
    control_dist = msn_ppc_control_corrs{i};
    msn_ppc_percentiles(i) = sum(control_dist <= msn_ppc_corr(i)) / length(control_dist) * 100;
    
    % STS percentile
    control_dist = msn_sts_control_corrs{i};
    msn_sts_percentiles(i) = sum(control_dist <= msn_sts_corr(i)) / length(control_dist) * 100;
end

% Calculate percentiles for FSIs
for i = 1:length(fsi_ppc_corr)
    % PPC percentile
    control_dist = fsi_ppc_control_corrs{i};
    fsi_ppc_percentiles(i) = sum(control_dist <= fsi_ppc_corr(i)) / length(control_dist) * 100;
    
    % STS percentile
    control_dist = fsi_sts_control_corrs{i};
    fsi_sts_percentiles(i) = sum(control_dist <= fsi_sts_corr(i)) / length(control_dist) * 100;
end

% Perform KS tests
[~, ks_p_ppc_corr, ~] = kstest2(fsi_ppc_corr, msn_ppc_corr, 'tail', 'smaller');
[~, ks_p_ppc_control_std, ~] = kstest2(fsi_ppc_control_std, msn_ppc_control_std, 'tail', 'larger');
[~, ks_p_ppc_control_mean, ~] = kstest2(fsi_ppc_control_mean, msn_ppc_control_mean, 'tail', 'smaller');

% Calculate differences (HFR-LFR correlation - mean control correlation)
msn_ppc_diff = msn_ppc_corr - msn_ppc_control_mean;
fsi_ppc_diff = fsi_ppc_corr - fsi_ppc_control_mean;

% Perform sign test (H0: HFR-LFR correlation >= mean control correlation)
% Count how many cells have HFR-LFR correlation >= mean control correlation
msn_sign_test = sum(msn_ppc_diff >= 0);
fsi_sign_test = sum(fsi_ppc_diff >= 0);
msn_total = length(msn_ppc_diff);
fsi_total = length(fsi_ppc_diff);

% Calculate p-values using binomial test
msn_p_value = 1 - binocdf(msn_sign_test - 1, msn_total, 0.5);
fsi_p_value = 1 - binocdf(fsi_sign_test - 1, fsi_total, 0.5);


% Colors
msn_color = [0.8500 0.3250 0.0980]; % Orange
fsi_color = [0 0.4470 0.7410]; % Blue

%% ===== FIGURE 1: DIAGNOSTIC PLOTS ===== 
fig = figure('WindowState', 'maximized');
fig.Renderer = 'painters';
fontname(fig, 'Helvetica');

% Row 1, Col 1: Box plot of mean PPC correlation
subplot(2,3,1);
% Create grouping variable for boxplot
group_labels = [repmat({'MSN'}, length(msn_ppc_corr), 1); repmat({'FSI'}, length(fsi_ppc_corr), 1)];
boxplot([msn_ppc_corr; fsi_ppc_corr], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_corr), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_corr), 1) - 0.5);
scatter(msn_x, msn_ppc_corr, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_corr, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC LFR-HFR Correlation');
title('Mean PPC Correlation');
ylim([-1.25 1.25]);
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', ks_p_ppc_corr), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

% Row 1, Col 2: Box plot of PPC correlation percentiles
subplot(2,3,2);
group_labels = [repmat({'MSN'}, length(msn_ppc_percentiles), 1); repmat({'FSI'}, length(fsi_ppc_percentiles), 1)];
boxplot([msn_ppc_percentiles; fsi_ppc_percentiles], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_percentiles), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_percentiles), 1) - 0.5);
scatter(msn_x, msn_ppc_percentiles, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_percentiles, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC Correlation Percentile');
title('PPC Correlation vs Control Distribution');
ylim([0 105]);

% Row 1, Col 3: Box plot of PPC control split std
subplot(2,3,3);
group_labels = [repmat({'MSN'}, length(msn_ppc_control_std), 1); repmat({'FSI'}, length(fsi_ppc_control_std), 1)];
boxplot([msn_ppc_control_std; fsi_ppc_control_std], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_control_std), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_control_std), 1) - 0.5);
scatter(msn_x, msn_ppc_control_std, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_std, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split SD');
title('PPC Control Split Variability');
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', ks_p_ppc_control_std), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');


% Row 2, Col 1: Box plot of mean STS correlation
subplot(2,3,4);
group_labels = [repmat({'MSN'}, length(msn_sts_corr), 1); repmat({'FSI'}, length(fsi_sts_corr), 1)];
boxplot([msn_sts_corr; fsi_sts_corr], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_sts_corr), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_sts_corr), 1) - 0.5);
scatter(msn_x, msn_sts_corr, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_sts_corr, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('STS LFR-HFR Correlation');
title('Mean STS Correlation');
ylim([-1.25 1.25]);

% Row 2, Col 2: Box plot of STS correlation percentiles
subplot(2,3,5);
group_labels = [repmat({'MSN'}, length(msn_sts_percentiles), 1); repmat({'FSI'}, length(fsi_sts_percentiles), 1)];
boxplot([msn_sts_percentiles; fsi_sts_percentiles], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_sts_percentiles), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_sts_percentiles), 1) - 0.5);
scatter(msn_x, msn_sts_percentiles, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_sts_percentiles, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('STS Correlation Percentile');
title('STS Correlation vs Control Distribution');
ylim([0 105]);

% Row 2, Col 3: Box plot of PPC control split mean
subplot(2,3,6);
group_labels = [repmat({'MSN'}, length(msn_ppc_control_mean), 1); repmat({'FSI'}, length(fsi_ppc_control_mean), 1)];
boxplot([msn_ppc_control_mean; fsi_ppc_control_mean], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_control_mean), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_control_mean), 1) - 0.5);
scatter(msn_x, msn_ppc_control_mean, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_mean, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split Mean');
title('PPC Control Split Mean');


%% ===== FIGURE 2: COMBINED LAYOUT =====

fig3 = figure('WindowState', 'maximized');
fig3.Renderer = 'painters';
fontname(fig3, 'Helvetica');

% Top row: First, third, and sixth subplots from the first figure
% Row 1, Col 1: Box plot of mean PPC correlation (from subplot 2,3,1)
subplot(2,3,1);
group_labels = [repmat({'MSN'}, length(msn_ppc_corr), 1); repmat({'FSI'}, length(fsi_ppc_corr), 1)];
boxplot([msn_ppc_corr; fsi_ppc_corr], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_corr), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_corr), 1) - 0.5);
scatter(msn_x, msn_ppc_corr, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_corr, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC LFR-HFR Correlation');
title('Mean PPC Correlation');
ylim([-0.75 1.25]);
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', ks_p_ppc_corr), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

% Row 1, Col 2: Box plot of PPC control split std (from subplot 2,3,3)
subplot(2,3,2);
group_labels = [repmat({'MSN'}, length(msn_ppc_control_std), 1); repmat({'FSI'}, length(fsi_ppc_control_std), 1)];
boxplot([msn_ppc_control_std; fsi_ppc_control_std], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_control_std), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_control_std), 1) - 0.5);
scatter(msn_x, msn_ppc_control_std, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_std, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split SD');
title('PPC Control Split Variability');
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', ks_p_ppc_control_std), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
% ylim([-1.25 1.25]);

% Row 1, Col 3: Box plot of PPC control split mean (from subplot 2,3,6)
subplot(2,3,3);
group_labels = [repmat({'MSN'}, length(msn_ppc_control_mean), 1); repmat({'FSI'}, length(fsi_ppc_control_mean), 1)];
boxplot([msn_ppc_control_mean; fsi_ppc_control_mean], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_control_mean), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_control_mean), 1) - 0.5);
scatter(msn_x, msn_ppc_control_mean, 30, msn_color, 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_mean, 30, fsi_color, 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split Mean');
title('PPC Control Split Mean');
ylim([-0.75 1.25]);

% Add KS test result to this subplot (using pre-calculated value)
text(0.05, 0.85, sprintf('KS Test: p = %.4f', ks_p_ppc_control_mean), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

% Bottom row: One big subplot spanning all 3 columns - the second summary plot
subplot(2,3,[4,5,6]);
hold on;

% Plot MSN data
for i = 1:length(msn_ppc_corr)
    % Connect HFR-LFR correlation to mean control correlation
    h = plot([1, 2], [msn_ppc_corr(i), msn_ppc_control_mean(i)], 'Color', msn_color, 'LineWidth', 0.5);
    h.Color(4) = 0.3; % Set alpha transparency
    % Scatter points
    scatter(1, msn_ppc_corr(i), 50, msn_color, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);
    scatter(2, msn_ppc_control_mean(i), 50, msn_color, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);
end

% Plot FSI data
for i = 1:length(fsi_ppc_corr)
    % Connect HFR-LFR correlation to mean control correlation
    h = plot([3, 4], [fsi_ppc_corr(i), fsi_ppc_control_mean(i)], 'Color', fsi_color, 'LineWidth', 0.5);
    h.Color(4) = 0.3; % Set alpha transparency
    % Scatter points
    scatter(3, fsi_ppc_corr(i), 50, fsi_color, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);
    scatter(4, fsi_ppc_control_mean(i), 50, fsi_color, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);
end

% Add mean lines
plot([1, 2], [mean(msn_ppc_corr), mean(msn_ppc_control_mean)], 'Color', msn_color, 'LineWidth', 3, 'LineStyle', '--');
plot([3, 4], [mean(fsi_ppc_corr), mean(fsi_ppc_control_mean)], 'Color', fsi_color, 'LineWidth', 3, 'LineStyle', '--');

% Customize plot
xlim([0.5, 4.5]);
ylim([-1.25, 1.25]);
xticks([1, 2, 3, 4]);
xticklabels({'MSN HFR-LFR', 'MSN Control Mean', 'FSI HFR-LFR', 'FSI Control Mean'});
ylabel('PPC Correlation');
title('HFR-LFR vs Control Mean PPC Correlation Comparison');
grid on;

% Add legend
legend_handles = [];
legend_handles(end+1) = scatter(NaN, NaN, 50, msn_color, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);
legend_handles(end+1) = scatter(NaN, NaN, 50, fsi_color, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);
legend_handles(end+1) = plot(NaN, NaN, 'Color', msn_color, 'LineWidth', 3, 'LineStyle', '--');
legend_handles(end+1) = plot(NaN, NaN, 'Color', fsi_color, 'LineWidth', 3, 'LineStyle', '--');
legend(legend_handles, {'MSN Individual', 'FSI Individual', 'MSN Mean', 'FSI Mean'}, 'Location', 'best');

% Add p-value text
text(0.05, 0.95, sprintf('MSN Sign Test: %d/%d cells ≥ control mean (p = %.4f)', msn_sign_test, msn_total, msn_p_value), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, 0.85, sprintf('FSI Sign Test: %d/%d cells ≥ control mean (p = %.4f)', fsi_sign_test, fsi_total, fsi_p_value), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

%% Print summary statistics
fprintf('\n=== CORRELATION SUMMARY STATISTICS ===\n');
fprintf('Firing Rate Difference Threshold: %.1f Hz\n', fr_diff_threshold);
fprintf('Total MSNs: %d, Eligible MSNs (after FR filter): %d\n', length(msn_summary.ppc_lfr_hfr_corr), length(msn_ppc_corr));
fprintf('Total FSIs: %d, Eligible FSIs (after FR filter): %d\n', length(fsi_summary.ppc_lfr_hfr_corr), length(fsi_ppc_corr));

fprintf('\n=== PPC Correlation Statistics ===\n');
fprintf('MSN: Mean = %.3f ± %.3f, Median = %.3f\n', mean(msn_ppc_corr), std(msn_ppc_corr), median(msn_ppc_corr));
fprintf('FSI: Mean = %.3f ± %.3f, Median = %.3f\n', mean(fsi_ppc_corr), std(fsi_ppc_corr), median(fsi_ppc_corr));
fprintf('MSN PPC Percentile: Mean = %.1f ± %.1f\n', mean(msn_ppc_percentiles), std(msn_ppc_percentiles));
fprintf('FSI PPC Percentile: Mean = %.1f ± %.1f\n', mean(fsi_ppc_percentiles), std(fsi_ppc_percentiles));
fprintf('MSN PPC Control SD: Mean = %.3f ± %.3f\n', mean(msn_ppc_control_std), std(msn_ppc_control_std));
fprintf('FSI PPC Control SD: Mean = %.3f ± %.3f\n', mean(fsi_ppc_control_std), std(fsi_ppc_control_std));
fprintf('MSN PPC Control Mean: Mean = %.3f ± %.3f\n', mean(msn_ppc_control_mean), std(msn_ppc_control_mean));
fprintf('FSI PPC Control Mean: Mean = %.3f ± %.3f\n', mean(fsi_ppc_control_mean), std(fsi_ppc_control_mean));

fprintf('\n=== STS Correlation Statistics ===\n');
fprintf('MSN: Mean = %.3f ± %.3f, Median = %.3f\n', mean(msn_sts_corr), std(msn_sts_corr), median(msn_sts_corr));
fprintf('FSI: Mean = %.3f ± %.3f, Median = %.3f\n', mean(fsi_sts_corr), std(fsi_sts_corr), median(fsi_sts_corr));
fprintf('MSN STS Percentile: Mean = %.1f ± %.1f\n', mean(msn_sts_percentiles), std(msn_sts_percentiles));
fprintf('FSI STS Percentile: Mean = %.1f ± %.1f\n', mean(fsi_sts_percentiles), std(fsi_sts_percentiles));

fprintf('\n=== KS Test Results ===\n');
fprintf('PPC Correlation (FSI < MSN): p = %.4f\n', ks_p_ppc_corr);
fprintf('PPC Control SD (FSI > MSN): p = %.4f\n', ks_p_ppc_control_std);
fprintf('PPC Control Mean KS Test (FSI < MSN): p = %.4f\n', ks_p_ppc_control_mean);

fprintf('\n=== Sign Test Results (HFR-LFR vs Control Mean) ===\n');
fprintf('MSN: %d/%d cells have HFR-LFR correlation ≥ control mean (p = %.4f)\n', msn_sign_test, msn_total, msn_p_value);
fprintf('FSI: %d/%d cells have HFR-LFR correlation ≥ control mean (p = %.4f)\n', fsi_sign_test, fsi_total, fsi_p_value);

%% ===== FIGURE 3: BOX PLOTS WITH JITTERED DATA =====
% Parameters for third figure
marker_size = 250; % Exposed parameter for marker size
marker_alpha = 0.2;
box_plot_width = 0.75;

% Create third figure
fig4 = figure('WindowState', 'maximized');
fig4.Renderer = 'painters';
fontname(fig4, 'Helvetica');

subplot(2,3,[1,4]);
hold on;

% Create box plot
boxplot([msn_ppc_corr; fsi_ppc_corr], [ones(length(msn_ppc_corr), 1); 2*ones(length(fsi_ppc_corr), 1)], ...
    'Labels', {'MSN', 'FSI'}, 'Colors', [msn_color; fsi_color], 'Widths', box_plot_width);

% Add jittered individual data points
jitter_msn = 1 + 0.6 * (rand(length(msn_ppc_corr), 1) - 0.5);
jitter_fsi = 2 + 0.6 * (rand(length(fsi_ppc_corr), 1) - 0.5);
scatter(jitter_msn, msn_ppc_corr, 'SizeData', marker_size, 'MarkerFaceColor', msn_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
scatter(jitter_fsi, fsi_ppc_corr, 'SizeData', marker_size, 'MarkerFaceColor', fsi_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);

% Line to mark means and draw lines
scatter(1, mean(msn_ppc_corr), 'SizeData', 3*marker_size, 'MarkerEdgeColor', msn_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
scatter(2, mean(fsi_ppc_corr), 'SizeData', 3*marker_size, 'MarkerEdgeColor', fsi_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
plot([1, 2], [mean(msn_ppc_corr), mean(fsi_ppc_corr)], 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--');
ylabel('PPC Correlation');

% Add hypothesis testing lines
y_pos = 1.1;
plot([1, 2], [y_pos, y_pos], 'k--', 'LineWidth', 2);
plot([1, 1], [y_pos-0.025, y_pos+0.025], 'k-', 'LineWidth', 2);
plot([2, 2], [y_pos-0.025, y_pos+0.025], 'k-', 'LineWidth', 2);
text(1.5, y_pos+0.1, sprintf('p = %.4f', ks_p_ppc_corr), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
ylim([-0.65 1.25]);
xlim([0.5 2.5]);
yticks([-0.5:0.25:1])
ax = gca;
ax.TickDir = 'out';
box off;

pause(3);
subplot(2,3,3);
hold on;

% Create box plot
boxplot([msn_ppc_control_mean; fsi_ppc_control_mean], [ones(length(msn_ppc_control_mean), 1); 2*ones(length(fsi_ppc_control_mean), 1)], ...
    'Labels', {'MSN', 'FSI'}, 'Colors', [msn_color; fsi_color], 'Widths', box_plot_width);

% Add jittered individual data points
jitter_msn = 1 + 0.6 * (rand(length(msn_ppc_control_mean), 1) - 0.5);
jitter_fsi = 2 + 0.6 * (rand(length(fsi_ppc_control_mean), 1) - 0.5);
scatter(jitter_msn, msn_ppc_control_mean, 'SizeData', marker_size, 'MarkerFaceColor', msn_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
scatter(jitter_fsi, fsi_ppc_control_mean, 'SizeData', marker_size, 'MarkerFaceColor', fsi_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);

% Line to mark means and draw lines
scatter(1, mean(msn_ppc_control_mean), 'SizeData', 3*marker_size, 'MarkerEdgeColor', msn_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
scatter(2, mean(fsi_ppc_control_mean), 'SizeData', 3*marker_size, 'MarkerEdgeColor', fsi_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
plot([1, 2], [mean(msn_ppc_control_mean), mean(fsi_ppc_control_mean)], 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--');
ylabel('PPC Control-split Correlation Mean');

% Add hypothesis testing lines
y_pos = 1.1;
plot([1, 2], [y_pos, y_pos], 'k--', 'LineWidth', 2);
plot([1, 1], [y_pos-0.05, y_pos+0.05], 'k-', 'LineWidth', 2);
plot([2, 2], [y_pos-0.05, y_pos+0.05], 'k-', 'LineWidth', 2);
text(1.5, y_pos+0.1, sprintf('p = %.4f', ks_p_ppc_control_mean), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
ylim([-0.4 1.25]);
xlim([0.5 2.5]);
yticks([-0.25:0.25:1])
ax = gca;
ax.TickDir = 'out';
box off;

pause(3);
subplot(2,3,6);
hold on;

% Create box plot
boxplot([msn_ppc_control_std; fsi_ppc_control_std], [ones(length(msn_ppc_control_std), 1); 2*ones(length(fsi_ppc_control_std), 1)], ...
    'Labels', {'MSN', 'FSI'}, 'Colors', [msn_color; fsi_color], 'Widths', box_plot_width);

% Add jittered individual data points
jitter_msn = 1 + 0.6 * (rand(length(msn_ppc_control_std), 1) - 0.5);
jitter_fsi = 2 + 0.6 * (rand(length(fsi_ppc_control_std), 1) - 0.5);
scatter(jitter_msn, msn_ppc_control_std, 'SizeData', marker_size, 'MarkerFaceColor', msn_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
scatter(jitter_fsi, fsi_ppc_control_std, 'SizeData', marker_size, 'MarkerFaceColor', fsi_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);

% Line to mark means and draw lines
scatter(1, mean(msn_ppc_control_std), 'SizeData', 3*marker_size, 'MarkerEdgeColor', msn_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
scatter(2, mean(fsi_ppc_control_std), 'SizeData', 3*marker_size, 'MarkerEdgeColor', fsi_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
plot([1, 2], [mean(msn_ppc_control_std), mean(fsi_ppc_control_std)], 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--');

ylabel('PPC Control-split Correlation SD');

% Add hypothesis testing lines
y_pos = 0.35;
plot([1 2], [y_pos, y_pos], 'k--', 'LineWidth', 2);
plot([1, 1], [y_pos-0.025, y_pos+0.02], 'k-', 'LineWidth', 2);
plot([2, 2], [y_pos-0.025, y_pos+0.02], 'k-', 'LineWidth', 2);
text(1.5, y_pos+0.05, sprintf('p = %.4f', ks_p_ppc_control_std), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
ylim([-0.05 0.4]);
xlim([0.5 2.5]);
yticks([0:0.1:0.35]);
ax = gca;
ax.TickDir = 'out';
box off;

pause(3);
subplot(2,3,[2,5]);
hold on

% Plot MSN data
for i = 1:length(msn_ppc_corr)
    % Connect HFR-LFR correlation to mean control correlation
    h = plot([1, 2], [msn_ppc_corr(i), msn_ppc_control_mean(i)], 'Color', msn_color, 'LineWidth', 0.5);
    h.Color(4) = 0.3; % Set alpha transparency
    % Scatter points with larger markers
    scatter(1, msn_ppc_corr(i), 'SizeData', marker_size, 'MarkerFaceColor', msn_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    scatter(2, msn_ppc_control_mean(i), 'SizeData', marker_size, 'MarkerFaceColor', msn_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

% Plot FSI data
for i = 1:length(fsi_ppc_corr)
    % Connect HFR-LFR correlation to mean control correlation
    h = plot([3, 4], [fsi_ppc_corr(i), fsi_ppc_control_mean(i)], 'Color', fsi_color, 'LineWidth', 0.5);
    h.Color(4) = 0.3; % Set alpha transparency
    % Scatter points with larger markers
    scatter(3, fsi_ppc_corr(i), 'SizeData', marker_size, 'MarkerFaceColor',  fsi_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
    scatter(4, fsi_ppc_control_mean(i), 'SizeData', marker_size, 'MarkerFaceColor',  fsi_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
end

% Add mean lines and circles
scatter(1, mean(msn_ppc_corr), 'SizeData', 3*marker_size, 'MarkerEdgeColor', msn_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
scatter(2, mean(msn_ppc_control_mean), 'SizeData', 3*marker_size, 'MarkerEdgeColor', msn_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
scatter(3, mean(fsi_ppc_corr), 'SizeData', 3*marker_size, 'MarkerEdgeColor',  fsi_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);
scatter(4, mean(fsi_ppc_control_mean), 'SizeData', 3*marker_size, 'MarkerEdgeColor',  fsi_color, 'MarkerEdgeAlpha', 2*marker_alpha, 'MarkerFaceAlpha', 0);

plot([1, 2], [mean(msn_ppc_corr), mean(msn_ppc_control_mean)], 'Color', msn_color, 'LineWidth', 3, 'LineStyle', '--');
plot([3, 4], [mean(fsi_ppc_corr), mean(fsi_ppc_control_mean)], 'Color', fsi_color, 'LineWidth', 3, 'LineStyle', '--');

xlabel('Condition');
ylabel('PPC Correlation');
title('HFR-LFR split vs Mean of Control Splits');
set(gca, 'Box', 'off', 'TickDir', 'out');
xticks([1, 2, 3, 4]);
xlim([0.75 4.25])
xticklabels({'MSN HFR-LFR', 'MSN Control', 'FSI HFR-LFR', 'FSI Control'});

% Add legend
legend_handles = [];
legend_entries = {};

% Create legend entries
legend_handles(end+1) = scatter(NaN, NaN, 'SizeData', marker_size, 'MarkerFaceColor',  msn_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
legend_entries{end+1} = sprintf('MSN (%d cells)', length(msn_ppc_corr));

legend_handles(end+1) = scatter(NaN, NaN, 'SizeData', marker_size, 'MarkerFaceColor',  fsi_color, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', 0);
legend_entries{end+1} = sprintf('FSI (%d cells)', length(fsi_ppc_corr));

legend_handles(end+1) = plot(NaN, NaN, 'Color', msn_color, 'LineWidth', 3, 'LineStyle', '--');
legend_entries{end+1} = sprintf('MSN Mean (%d/%d ≥ control, p=%.3f)', msn_sign_test, msn_total, msn_p_value);

legend_handles(end+1) = plot(NaN, NaN, 'Color', fsi_color, 'LineWidth', 3, 'LineStyle', '--');
legend_entries{end+1} = sprintf('FSI Mean (%d/%d ≥ control, p=%.3f)', fsi_sign_test, fsi_total, fsi_p_value);

legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 12);

ax = gca;
ax.TickDir = 'out';
box off;

