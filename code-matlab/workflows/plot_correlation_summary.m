%% Load summary files
load('D:\vStr_oscillatory_switch_results\cellFingerPrintsWithCorrelation\msn_summary.mat');
load('D:\vStr_oscillatory_switch_results\cellFingerPrintsWithCorrelation\fsi_summary.mat');

%% Extract correlation data
% MSN data - these are already numeric, no need for cell2mat
msn_ppc_corr = msn_summary.ppc_lfr_hfr_corr;
msn_sts_corr = msn_summary.sts_lfr_hfr_corr;
msn_ppc_control_corrs = msn_summary.ppc_control_corrs;
msn_sts_control_corrs = msn_summary.sts_control_corrs;

% FSI data - these are already numeric, no need for cell2mat
fsi_ppc_corr = fsi_summary.ppc_lfr_hfr_corr;
fsi_sts_corr = fsi_summary.sts_lfr_hfr_corr;
fsi_ppc_control_corrs = fsi_summary.ppc_control_corrs;
fsi_sts_control_corrs = fsi_summary.sts_control_corrs;

%% Calculate control split statistics
% Calculate std and mean of control correlations for each cell
msn_ppc_control_std = cellfun(@(x) std(x), msn_ppc_control_corrs);
msn_ppc_control_mean = cellfun(@(x) mean(x), msn_ppc_control_corrs);
msn_sts_control_std = cellfun(@(x) std(x), msn_sts_control_corrs);
msn_sts_control_mean = cellfun(@(x) mean(x), msn_sts_control_corrs);

fsi_ppc_control_std = cellfun(@(x) std(x), fsi_ppc_control_corrs);
fsi_ppc_control_mean = cellfun(@(x) mean(x), fsi_ppc_control_corrs);
fsi_sts_control_std = cellfun(@(x) std(x), fsi_sts_control_corrs);
fsi_sts_control_mean = cellfun(@(x) mean(x), fsi_sts_control_corrs);

%% Calculate percentiles
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

%% Perform KS tests
[h1, p1, ~] = kstest2(fsi_ppc_corr, msn_ppc_corr, 'tail', 'smaller');
[h2, p2, ~] = kstest2(fsi_ppc_control_std, msn_ppc_control_std, 'tail', 'larger');

%% Create 2x3 plots 
fig = figure('WindowState', 'maximized');

% Colors
msn_color = [0.8500 0.3250 0.0980]; % Orange
fsi_color = [0 0.4470 0.7410]; % Blue

% Row 1, Col 1: Box plot of mean PPC correlation
subplot(2,3,1);
% Create grouping variable for boxplot
group_labels = [repmat({'MSN'}, length(msn_ppc_corr), 1); repmat({'FSI'}, length(fsi_ppc_corr), 1)];
boxplot([msn_ppc_corr; fsi_ppc_corr], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_corr), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_corr), 1) - 0.5);
scatter(msn_x, msn_ppc_corr, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_corr, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('PPC LFR-HFR Correlation');
title('Mean PPC Correlation');
ylim([-1.25 1.25]);
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', p1), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

% Row 1, Col 2: Box plot of PPC correlation percentiles
subplot(2,3,2);
group_labels = [repmat({'MSN'}, length(msn_ppc_percentiles), 1); repmat({'FSI'}, length(fsi_ppc_percentiles), 1)];
boxplot([msn_ppc_percentiles; fsi_ppc_percentiles], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_percentiles), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_percentiles), 1) - 0.5);
scatter(msn_x, msn_ppc_percentiles, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_percentiles, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
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
scatter(msn_x, msn_ppc_control_std, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_std, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split Std');
title('PPC Control Split Variability');
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', p2), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

% Row 2, Col 1: Box plot of mean STS correlation
subplot(2,3,4);
group_labels = [repmat({'MSN'}, length(msn_sts_corr), 1); repmat({'FSI'}, length(fsi_sts_corr), 1)];
boxplot([msn_sts_corr; fsi_sts_corr], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_sts_corr), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_sts_corr), 1) - 0.5);
scatter(msn_x, msn_sts_corr, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_sts_corr, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
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
scatter(msn_x, msn_sts_percentiles, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_sts_percentiles, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
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
scatter(msn_x, msn_ppc_control_mean, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_mean, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split Mean');
title('PPC Control Split Mean');


%% Create summary figure: Combined layout

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
fig3 = figure('WindowState', 'maximized');

% Top row: First, third, and sixth subplots from the first figure
% Row 1, Col 1: Box plot of mean PPC correlation (from subplot 2,3,1)
subplot(2,3,1);
group_labels = [repmat({'MSN'}, length(msn_ppc_corr), 1); repmat({'FSI'}, length(fsi_ppc_corr), 1)];
boxplot([msn_ppc_corr; fsi_ppc_corr], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_corr), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_corr), 1) - 0.5);
scatter(msn_x, msn_ppc_corr, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_corr, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('PPC LFR-HFR Correlation');
title('Mean PPC Correlation');
ylim([-0.75 1.25]);
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', p1), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

% Row 1, Col 2: Box plot of PPC control split std (from subplot 2,3,3)
subplot(2,3,2);
group_labels = [repmat({'MSN'}, length(msn_ppc_control_std), 1); repmat({'FSI'}, length(fsi_ppc_control_std), 1)];
boxplot([msn_ppc_control_std; fsi_ppc_control_std], group_labels, 'Colors', [msn_color; fsi_color]);
hold on;
% Add jittered scatter plots
msn_x = 0.8 + 0.4 * (rand(length(msn_ppc_control_std), 1) - 0.5);
fsi_x = 1.8 + 0.4 * (rand(length(fsi_ppc_control_std), 1) - 0.5);
scatter(msn_x, msn_ppc_control_std, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_std, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split Std');
title('PPC Control Split Variability');
% Add p-value text
text(0.05, 0.95, sprintf('p = %.4f', p2), 'Units', 'normalized', ...
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
scatter(msn_x, msn_ppc_control_mean, 30, msn_color, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(fsi_x, fsi_ppc_control_mean, 30, fsi_color, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('PPC Control Split Mean');
title('PPC Control Split Mean');
ylim([-0.75 1.25]);

% Bottom row: One big subplot spanning all 3 columns - the second summary plot
subplot(2,3,[4,5,6]);
hold on;

% Plot MSN data
for i = 1:length(msn_ppc_corr)
    % Connect HFR-LFR correlation to mean control correlation
    h = plot([1, 2], [msn_ppc_corr(i), msn_ppc_control_mean(i)], 'Color', msn_color, 'LineWidth', 0.5);
    h.Color(4) = 0.3; % Set alpha transparency
    % Scatter points
    scatter(1, msn_ppc_corr(i), 50, msn_color, 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(2, msn_ppc_control_mean(i), 50, msn_color, 'filled', 'MarkerFaceAlpha', 0.7);
end

% Plot FSI data
for i = 1:length(fsi_ppc_corr)
    % Connect HFR-LFR correlation to mean control correlation
    h = plot([3, 4], [fsi_ppc_corr(i), fsi_ppc_control_mean(i)], 'Color', fsi_color, 'LineWidth', 0.5);
    h.Color(4) = 0.3; % Set alpha transparency
    % Scatter points
    scatter(3, fsi_ppc_corr(i), 50, fsi_color, 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(4, fsi_ppc_control_mean(i), 50, fsi_color, 'filled', 'MarkerFaceAlpha', 0.7);
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
legend_handles(end+1) = scatter(NaN, NaN, 50, msn_color, 'filled', 'MarkerFaceAlpha', 0.7);
legend_handles(end+1) = scatter(NaN, NaN, 50, fsi_color, 'filled', 'MarkerFaceAlpha', 0.7);
legend_handles(end+1) = plot(NaN, NaN, 'Color', msn_color, 'LineWidth', 3, 'LineStyle', '--');
legend_handles(end+1) = plot(NaN, NaN, 'Color', fsi_color, 'LineWidth', 3, 'LineStyle', '--');
legend(legend_handles, {'MSN Individual', 'FSI Individual', 'MSN Mean', 'FSI Mean'}, 'Location', 'best');

% Add p-value text
text(0.05, 0.95, sprintf('MSN Sign Test: %d/%d cells ≥ control mean (p = %.4f)', msn_sign_test, msn_total, msn_p_value), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, 0.85, sprintf('FSI Sign Test: %d/%d cells ≥ control mean (p = %.4f)', fsi_sign_test, fsi_total, fsi_p_value), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');