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

% Updated headers with peak magnitudes
headers = {'label', 'lfr_min', 'lfr_max', 'lfr_mean', 'hfr_min', ...
    'hfr_max', 'hfr_mean', 'lfr_sts_peak', 'lfr_sts_diff', ...
    'hfr_sts_peak', 'hfr_sts_diff', 'lfr_ppc_peak', 'lfr_ppc_diff', ...
    'hfr_ppc_peak', 'hfr_ppc_diff', 'lfr_sts_peak_mag', 'hfr_sts_peak_mag', ...
    'lfr_ppc_peak_mag', 'hfr_ppc_peak_mag', 'lfr_sts_overall_peak', ...
    'hfr_sts_overall_peak', 'lfr_sts_overall_peak_mag', 'hfr_sts_overall_peak_mag'};

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
                headers);
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
%% Old summary viz 
% Scatter plot of MSNs and FSIs that have both significantly phase-locked HFR and LFR peaks
only_hfr_msn = ~isnan(msn_summary.hfr_ppc_peak);
only_lfr_msn = ~isnan(msn_summary.lfr_ppc_peak);
clean_msn = only_hfr_msn & only_lfr_msn;

only_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_peak);
only_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_peak);
clean_fsi = only_hfr_fsi & only_lfr_fsi;

dif_hfr_msn = ~isnan(msn_summary.hfr_ppc_diff);
dif_lfr_msn = ~isnan(msn_summary.lfr_ppc_diff);
dif_msn = dif_hfr_msn | dif_lfr_msn;
dif_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_diff);
dif_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_dn b  iff);
dif_fsi = dif_hfr_fsi | dif_lfr_fsi;

bdif_msn = dif_hfr_msn & dif_lfr_msn;
bdif_fsi = dif_hfr_fsi & dif_lfr_fsi;

fig = figure('WindowState', 'maximized');
ax1 = subplot(2,2,1);
s1 = scatter(ax1,(msn_summary.hfr_mean(clean_msn) - msn_summary.lfr_mean(clean_msn)), ...
    (msn_summary.hfr_ppc_peak(clean_msn) - msn_summary.lfr_ppc_peak(clean_msn)));
s1.Marker = 'o';
s1.MarkerFaceColor = c4;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
hold on

s1 = scatter(ax1, (fsi_summary.hfr_mean(clean_fsi) - fsi_summary.lfr_mean(clean_fsi)), ...
    (fsi_summary.hfr_ppc_peak(clean_fsi) - fsi_summary.lfr_ppc_peak(clean_fsi)));
s1.Marker = 'o';
s1.MarkerFaceColor = 'blue';
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;

% Outline significant MSN diffs
s1 = scatter(ax1,(msn_summary.hfr_mean(dif_msn) - msn_summary.lfr_mean(dif_msn)), ...
    (msn_summary.hfr_ppc_peak(dif_msn) - msn_summary.lfr_ppc_peak(dif_msn)));
s1.Marker = 'o';
s1.MarkerEdgeColor = 'black';
s1.MarkerFaceAlpha = 0;
s1.MarkerEdgeAlpha = 1;
s1.SizeData = 100;

% Outline significant FSI diffs

s1 = scatter(ax1,(fsi_summary.hfr_mean(dif_fsi) - fsi_summary.lfr_mean(dif_fsi)), ...
    (fsi_summary.hfr_ppc_peak(dif_fsi) - fsi_summary.lfr_ppc_peak(dif_fsi)));
s1.Marker = 'o';
s1.MarkerEdgeColor = 'black';
s1.MarkerFaceAlpha = 0;
s1.MarkerEdgeAlpha = 1;
s1.SizeData = 100;

ax1.XLabel.String = '\Delta F.R';
ax1.YLabel.String = '\Delta Freq';
ax1.FontSize = 16;

legend({sprintf('MSNs: %d/%d, diff: %d ', sum(clean_msn), length(clean_msn), sum(dif_msn)), ...
    sprintf('FSIs: %d/%d, diff: %d ', sum(clean_fsi), length(clean_fsi), sum(dif_fsi))}, 'FontSize', 14, 'Location', 'best');
% Scatter plot of MSNs with significant diffs
ax1 = subplot(2,2,2);
hold off;
s1 = scatter(ax1,(msn_summary.hfr_mean(dif_hfr_msn) - msn_summary.lfr_mean(dif_hfr_msn)), ...
    msn_summary.hfr_ppc_diff(dif_hfr_msn));
s1.Marker = 'o';
s1.MarkerFaceColor = c2;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
hold on
s1 = scatter(ax1,(msn_summary.hfr_mean(dif_lfr_msn) - msn_summary.lfr_mean(dif_lfr_msn)), ...
    msn_summary.lfr_ppc_diff(dif_lfr_msn));
s1.Marker = 'o';
s1.MarkerFaceColor = c1;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
    
ax1.XLabel.String = '\Delta F.R';
ax1.YLabel.String = 'Freq. of sig difference';
ax1.FontSize = 16;
ax1.Title.String = 'Sig diff MSNs';

legend({sprintf('sig HFR: %d/%d', sum(dif_hfr_msn), sum(clean_msn)), ...
    sprintf('sig LFR: %d/%d', sum(dif_lfr_msn), sum(clean_msn))}, 'FontSize', 14, 'Location', 'best')

% Scatter plot of FSIs with significant diffs
ax1 = subplot(2,2,3);
hold off;
s1 = scatter(ax1,(fsi_summary.hfr_mean(dif_hfr_fsi) - fsi_summary.lfr_mean(dif_hfr_fsi)), ...
    fsi_summary.hfr_ppc_diff(dif_hfr_fsi));
s1.Marker = 'o';
s1.MarkerFaceColor = c2;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
hold on
s1 = scatter(ax1,(fsi_summary.hfr_mean(dif_lfr_fsi) - fsi_summary.lfr_mean(dif_lfr_fsi)), ...
    fsi_summary.lfr_ppc_diff(dif_lfr_fsi));
s1.Marker = 'o';
s1.MarkerFaceColor = c1;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
    
ax1.XLabel.String = '\Delta F.R';
ax1.YLabel.String = 'Freq. of sig difference';
ax1.FontSize = 16;
ax1.Title.String = 'Sig diff FSIs';

legend({sprintf('sig HFR: %d/%d', sum(dif_hfr_fsi), sum(clean_fsi)), ...
    sprintf('sig LFR: %d/%d', sum(dif_lfr_fsi), sum(clean_fsi))}, 'FontSize', 14, 'Location', 'best')

% Scatter plot of both sig
ax1 = subplot(2,2,4);
hold off
s1 = scatter(ax1,(msn_summary.hfr_mean(bdif_msn) - msn_summary.lfr_mean(bdif_msn)), ...
    (msn_summary.hfr_ppc_diff(bdif_msn) - msn_summary.lfr_ppc_diff(bdif_msn)));
s1.Marker = 'o';
s1.MarkerFaceColor = c4;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
hold on

s1 = scatter(ax1, (fsi_summary.hfr_mean(bdif_fsi) - fsi_summary.lfr_mean(bdif_fsi)), ...
    (fsi_summary.hfr_ppc_diff(bdif_fsi) - fsi_summary.lfr_ppc_diff(bdif_fsi)));
s1.Marker = 'o';
s1.MarkerFaceColor = 'blue';
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;

ax1.XLabel.String = '\Delta F.R';
ax1.YLabel.String = '\Delta Freq Sig';
ax1.FontSize = 16;
ax1.YLim = [-100 100];

legend({sprintf('MSNs: %d', sum(bdif_msn)), sprintf('FSIs: %d', sum(bdif_fsi))}, ...
    'FontSize', 14, 'Location', 'best');

%% Apply masks for various cases
% Define colors (same scheme as before)
c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
c2 = [26/255 255/255 26/255]; % Green for HFR
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange

% Define masks
only_hfr_msn = ~isnan(msn_summary.hfr_ppc_peak);
only_lfr_msn = ~isnan(msn_summary.lfr_ppc_peak);
clean_msn = only_hfr_msn & only_lfr_msn;

only_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_peak);
only_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_peak);
clean_fsi = only_hfr_fsi & only_lfr_fsi;

dif_hfr_msn = ~isnan(msn_summary.hfr_ppc_diff);
dif_lfr_msn = ~isnan(msn_summary.lfr_ppc_diff);
dif_hfr_fsi = ~isnan(fsi_summary.hfr_ppc_diff);
dif_lfr_fsi = ~isnan(fsi_summary.lfr_ppc_diff);

%% New summary viz 1
% Scatter plot of both sig

fig = figure('WindowState', 'maximized');
ax1 = gca;
s1 = scatter(ax1,(msn_summary.hfr_mean(bdif_msn) - msn_summary.lfr_mean(bdif_msn)), ...
    (msn_summary.hfr_ppc_diff(bdif_msn) - msn_summary.lfr_ppc_diff(bdif_msn)));
s1.Marker = 'o';
s1.MarkerFaceColor = c4;
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;
hold on

s1 = scatter(ax1, (fsi_summary.hfr_mean(bdif_fsi) - fsi_summary.lfr_mean(bdif_fsi)), ...
    (fsi_summary.hfr_ppc_diff(bdif_fsi) - fsi_summary.lfr_ppc_diff(bdif_fsi)));
s1.Marker = 'o';
s1.MarkerFaceColor = 'blue';
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeAlpha = 0;
s1.SizeData = 100;

ax1.XLabel.String = '\Delta F.R';
ax1.YLabel.String = '\Delta Freq Sig';
ax1.FontSize = 16;
ax1.YLim = [-100 100];

legend({sprintf('MSNs: %d', sum(bdif_msn)), sprintf('FSIs: %d', sum(bdif_fsi))}, ...
    'FontSize', 14, 'Location', 'best');

%% New summary viz 2
% Define bin edges for 10 bins from 0-100 (each bin is 10 Hz wide)
bin_edges = 0:10:100;

% Create figure
fig = figure('WindowState', 'maximized');

% Left subplot - MSNs
subplot(1,2,1);
hold on;
h1 = histogram(msn_summary.hfr_ppc_diff(dif_hfr_msn), bin_edges, 'FaceColor', c2, 'FaceAlpha', 0.25, 'EdgeColor', c2);
h2 = histogram(msn_summary.lfr_ppc_diff(dif_lfr_msn), bin_edges, 'FaceColor', c1, 'FaceAlpha', 0.25, 'EdgeColor', c1);
xlabel('Frequency (Hz)');
ylabel('Count');
title('MSN PPC Diff Frequencies');
legend({sprintf('HFR: n=%d', sum(dif_hfr_msn)), sprintf('LFR: n=%d', sum(dif_lfr_msn))}, 'Location', 'best');
% grid on;

% Right subplot - FSIs
subplot(1,2,2);
hold on;
h3 = histogram(fsi_summary.hfr_ppc_diff(dif_hfr_fsi), bin_edges, 'FaceColor', c2, 'FaceAlpha', 0.25, 'EdgeColor', c2);
h4 = histogram(fsi_summary.lfr_ppc_diff(dif_lfr_fsi), bin_edges, 'FaceColor', c1, 'FaceAlpha', 0.25, 'EdgeColor', c1);
xlabel('Frequency (Hz)');
ylabel('Count');
title('FSI PPC Diff Frequencies');
legend({sprintf('HFR: n=%d', sum(dif_hfr_fsi)), sprintf('LFR: n=%d', sum(dif_lfr_fsi))}, 'Location', 'best');
% grid on;

%% New summary viz 3 (no thresholding)
% Define colors (same scheme as before)
c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
c2 = [26/255 255/255 26/255]; % Green for HFR
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange

% Calculate the contrast index and mean frequency for MSNs
msn_hfr_mag = msn_summary.hfr_sts_overall_peak_mag;
msn_lfr_mag = msn_summary.lfr_sts_overall_peak_mag;
msn_hfr_freq = msn_summary.hfr_sts_overall_peak;
msn_lfr_freq = msn_summary.lfr_sts_overall_peak;

% Contrast index: (HFR - LFR) / (HFR + LFR)
msn_contrast = (msn_hfr_mag - msn_lfr_mag) ./ (msn_hfr_mag + msn_lfr_mag);
% Mean frequency
msn_mean_freq = (msn_hfr_freq + msn_lfr_freq) / 2;
msn_mean_freq = (msn_hfr_freq - msn_lfr_freq);

% Calculate the contrast index and mean frequency for FSIs
fsi_hfr_mag = fsi_summary.hfr_sts_overall_peak_mag;
fsi_lfr_mag = fsi_summary.lfr_sts_overall_peak_mag;
fsi_hfr_freq = fsi_summary.hfr_sts_overall_peak;
fsi_lfr_freq = fsi_summary.lfr_sts_overall_peak;

% Contrast index: (HFR - LFR) / (HFR + LFR)
fsi_contrast = (fsi_hfr_mag - fsi_lfr_mag) ./ (fsi_hfr_mag + fsi_lfr_mag);
% Mean frequency
fsi_mean_freq = (fsi_hfr_freq + fsi_lfr_freq) / 2;
fsi_mean_freq = (fsi_hfr_freq - fsi_lfr_freq);

% Remove NaN values for plotting
msn_valid = ~isnan(msn_contrast) & ~isnan(msn_mean_freq);
fsi_valid = ~isnan(fsi_contrast) & ~isnan(fsi_mean_freq);

% Create figure
fig = figure('WindowState', 'maximized');

% Left subplot - MSNs
subplot(1,2,1);
scatter(msn_mean_freq(msn_valid),msn_contrast(msn_valid), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('STS Magnitude Contrast ((HFR-LFR)/(HFR+LFR))');
xlabel('Peak Frequency Diff (Hz)');
title(sprintf('MSN STS Peak Analysis (n=%d)', sum(msn_valid)));
ylim([-1 1]);
xlim([-100,100]);

% Right subplot - FSIs
subplot(1,2,2);
scatter(fsi_mean_freq(fsi_valid), fsi_contrast(fsi_valid), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
ylabel('STS Magnitude Contrast ((HFR-LFR)/(HFR+LFR))');
xlabel('Peak Frequency Diff (Hz)');
title(sprintf('FSI STS Peak Analysis (n=%d)', sum(fsi_valid)));
ylim([-1 1]);
xlim([-100,100]);
%% New summary thresholding viz 4 (With thresholding)

% Define colors (same scheme as before)
c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
c2 = [26/255 255/255 26/255]; % Green for HFR
c3 = [0.7 0.7 0.7]; % Gray
c4 = [0.8500 0.3250 0.0980]; % Orange

% Calculate the contrast index and mean frequency for MSNs (thresholded peaks)
msn_hfr_mag = msn_summary.hfr_sts_peak_mag;
msn_lfr_mag = msn_summary.lfr_sts_peak_mag;
msn_hfr_freq = msn_summary.hfr_sts_peak;
msn_lfr_freq = msn_summary.lfr_sts_peak;

% Contrast index: (HFR - LFR) / (HFR + LFR)
msn_contrast = (msn_hfr_mag - msn_lfr_mag) ./ (msn_hfr_mag + msn_lfr_mag);
% Mean frequency
msn_mean_freq = (msn_hfr_freq + msn_lfr_freq) / 2;
msn_mean_freq = (msn_hfr_freq + msn_lfr_freq);

% Calculate the contrast index and mean frequency for FSIs (thresholded peaks)
fsi_hfr_mag = fsi_summary.hfr_sts_peak_mag;
fsi_lfr_mag = fsi_summary.lfr_sts_peak_mag;
fsi_hfr_freq = fsi_summary.hfr_sts_peak;
fsi_lfr_freq = fsi_summary.lfr_sts_peak;

% Contrast index: (HFR - LFR) / (HFR + LFR)
fsi_contrast = (fsi_hfr_mag - fsi_lfr_mag) ./ (fsi_hfr_mag + fsi_lfr_mag);
% Mean frequency
fsi_mean_freq = (fsi_hfr_freq + fsi_lfr_freq) / 2;
fsi_mean_freq = (fsi_hfr_freq - fsi_lfr_freq);

% Only include cases where BOTH HFR and LFR thresholded peaks exist
msn_valid = ~isnan(msn_hfr_mag) & ~isnan(msn_lfr_mag) & ~isnan(msn_hfr_freq) & ~isnan(msn_lfr_freq);
fsi_valid = ~isnan(fsi_hfr_mag) & ~isnan(fsi_lfr_mag) & ~isnan(fsi_hfr_freq) & ~isnan(fsi_lfr_freq);

% Create figure
fig = figure('WindowState', 'maximized');

% Left subplot - MSNs
subplot(1,2,1);
scatter(msn_mean_freq(msn_valid), msn_contrast(msn_valid), 100, c4, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Peak Frequency Diff (Hz)');
ylabel('STS Magnitude Contrast ((HFR-LFR)/(HFR+LFR))');
title(sprintf('MSN STS Thresholded Peak Analysis (n=%d)', sum(msn_valid)));
ylim([-1 1]);
xlim([-100,100]);

% Right subplot - FSIs
subplot(1,2,2);
scatter(fsi_mean_freq(fsi_valid), fsi_contrast(fsi_valid), 100, 'blue', 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Peak Frequency Diff (Hz)');
ylabel('STS Magnitude Contrast ((HFR-LFR)/(HFR+LFR))');
title(sprintf('FSI STS Thresholded Peak Analysis (n=%d)', sum(fsi_valid)));
ylim([-1 1]);
xlim([-100,100]);
%% Other functions
% Modular function to process each cell type
function [summary_table, clean_count] = process_cell_type(cell_res, cell_labels, cell_type, colors, params, headers)
    % Extract colors
    c1 = colors.c1; c2 = colors.c2; c3 = colors.c3; c4 = colors.c4;
    
    % Extract parameters
    min_freq = params.min_freq;
    pl_thresh = params.pl_thresh;
    diff_thresh = params.diff_thresh;
    peak_freq_tol_win = params.peak_freq_tol_win;
    odir = params.odir;
    
    % Initialize results
    summary_results = cell(0, 23);  % Pre-allocate with correct number of columns
    clean_count = 0;
    
    for iC = 1:length(cell_labels)
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
            [sts_peaks, sts_peak_mags, sts_overall_peaks, sts_overall_peak_mags] = plot_sts_and_find_peaks(cell_res, iC, c1, c2, c3, x1, this_freqs, pl_thresh, cell_type);
            
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
            
            % Compile results row
            this_row = {cell_labels{iC}, firing_rate_stats.lfr_min, firing_rate_stats.lfr_max, ...
                firing_rate_stats.lfr_mean, firing_rate_stats.hfr_min, firing_rate_stats.hfr_max, firing_rate_stats.hfr_mean, ...
                sts_peaks.lfr, sts_diff_peaks.lfr, sts_peaks.hfr, sts_diff_peaks.hfr, ...
                ppc_peaks.lfr, ppc_diff_peaks.lfr, ppc_peaks.hfr, ppc_diff_peaks.hfr, ...
                sts_peak_mags.lfr, sts_peak_mags.hfr, ppc_peak_mags.lfr, ppc_peak_mags.hfr, ...
                sts_overall_peaks.lfr, sts_overall_peaks.hfr, sts_overall_peak_mags.lfr, sts_overall_peak_mags.hfr};
            
            % Debug: Check row length
            fprintf('Row %d length: %d\n', clean_count, length(this_row));
            summary_results(end+1, :) = this_row;
        end
    end
    
    % Convert to table
    if ~isempty(summary_results)
        % Debug: Check dimensions
        fprintf('Summary results size: %d x %d, Headers length: %d\n', size(summary_results,1), size(summary_results,2), length(headers));
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

% Helper function: Plot STS and find peaks
function [peaks, peak_mags, overall_peaks, overall_peak_mags] = plot_sts_and_find_peaks(cell_res, iC, c1, c2, c3, x1, this_freqs, pl_thresh, cell_type)
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
    
    % Initialize overall peak outputs
    overall_peaks.lfr = NaN; overall_peaks.hfr = NaN;
    overall_peak_mags.lfr = NaN; overall_peak_mags.hfr = NaN;
    
    % Find LFR overall peak (maximum regardless of threshold)
    [~, peak_idx] = max(lfr_sts(x1:end));
    overall_peaks.lfr = this_freqs(peak_idx);
    overall_peak_mags.lfr = lfr_sts(x1-1+peak_idx);
    xline(ax2, overall_peaks.lfr, 'Color', c1, 'LineStyle', ':', 'LineWidth', 1.5);
    % Add horizontal line with dotted style for overall peaks
    yline(ax2, overall_peak_mags.lfr, 'Color', c1, 'LineStyle', ':', 'Alpha', 0.7, 'LineWidth', 1.5);
    
    % Find HFR overall peak (maximum regardless of threshold)
    [~, peak_idx] = max(hfr_sts(x1:end));
    overall_peaks.hfr = this_freqs(peak_idx);
    overall_peak_mags.hfr = hfr_sts(x1-1+peak_idx);
    xline(ax2, overall_peaks.hfr, 'Color', c2, 'LineStyle', ':', 'LineWidth', 1.5);
    % Add horizontal line with dotted style for overall peaks
    yline(ax2, overall_peak_mags.hfr, 'Color', c2, 'LineStyle', ':', 'Alpha', 0.7, 'LineWidth', 1.5);
    
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
