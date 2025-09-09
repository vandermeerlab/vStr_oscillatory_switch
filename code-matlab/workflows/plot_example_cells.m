%% Example Cell Plot Script
% cell names to plot
cell_names = {'R117-2007-06-01-TT07_1', ...
    'R132-2007-10-27-TT05_3', ... 
    'R117-2007-06-08-TT03_3', ...
    'R119-2007-06-29-TT06_5'}; 
% Specify the directory where results are stored
results_dir = 'D:\vStr_oscillatory_switch_results\temp3\'; % Change this to your results directory

% Output directory
output_dir = 'E:\Dartmouth College Dropbox\Manish Mohapatra\Figures\Papers\vStr_oscillatory_switch\version1\ExampleCells\'; % Change this to your output directory


% Define colors
c1 = [75/255 0/255 146/255];  % Violet/Purple for LFR
c2 = [26/255 255/255 26/255]; % Green for HFR
c3 = [0.7 0.7 0.7]; % Gray for All Trials

% Plotting parameters
line_width = 1;
tick_font_size = 20;

% Convert single cell name to cell array if needed
if ischar(cell_names)
    cell_names = {cell_names};
end

% First, check which cells exist and count them
valid_cells = {};
valid_cell_info = {};

for cell_idx = 1:length(cell_names)
    current_cell = cell_names{cell_idx};
    
    % Extract session info from cell name (format: Rat-Date-TT-Cell)
    cell_parts = strsplit(current_cell, '-');
    rat_name = cell_parts{1};
    date_part = strjoin(cell_parts(2:4),'-');
    tt_part = cell_parts{5};
    % Extract the cell number from the TT part (e.g., "TT07_1" -> "TT07")
    tt_cell_parts = strsplit(tt_part, '_');
    tt_name = tt_cell_parts{1};
     
    % Find the corresponding session file
    search_pattern = strcat(rat_name, '-', date_part, '_','ft_spec.mat');
    session_files = dir(fullfile(results_dir, search_pattern));
    
    if isempty(session_files)
        warning('No session file found for cell %s. Skipping...', current_cell);
        continue;
    end
    
    % Load the session data
    session_data = load(fullfile(results_dir, session_files(1).name));
    od = session_data.od;
    
    % Find the cell index in the session data
    cell_found = false;
    cell_index = [];
    
    % Check in MSN results
    if isfield(od, 'msn_res')
        msn_labels = od.label(od.cell_type == 1);
        msn_labels_clean = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        
        for i = 1:length(msn_labels_clean)
            if strcmp(msn_labels_clean{i}, current_cell)
                cell_index = i;
                cell_type = 'MSN';
                cell_res = od.msn_res;
                cell_found = true;
                break;
            end
        end
    end
    
    % Check in FSI results if not found in MSN
    if ~cell_found && isfield(od, 'fsi_res')
        fsi_labels = od.label(od.cell_type == 2);
        fsi_labels_clean = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        
        for i = 1:length(fsi_labels_clean)
            if strcmp(fsi_labels_clean{i}, current_cell)
                cell_index = i;
                cell_type = 'FSI';
                cell_res = od.fsi_res;
                cell_found = true;
                break;
            end
        end
    end
    
    if cell_found
        valid_cells{end+1} = current_cell;
        valid_cell_info{end+1} = struct('session_data', session_data, 'cell_index', cell_index, ...
            'cell_type', cell_type, 'cell_res', cell_res, 'session_file', session_files(1).name);
        fprintf('Found cell: %s (%s)\n', current_cell, cell_type);
    else
        warning('Cell %s not found in session data. Skipping...', current_cell);
    end
end

% Create figure with appropriate number of subplots
n_cells = length(valid_cells);
fig = figure('WindowState', 'maximized');

% Store axes handles for easy access
sta_axes = cell(1, n_cells);
ppc_axes = cell(1, n_cells);

% Process each valid cell
for cell_idx = 1:n_cells
    current_cell = valid_cells{cell_idx};
    cell_info = valid_cell_info{cell_idx};
    
    % Plot STA
    ax1 = subplot(2, n_cells, cell_idx);
    sta_axes{cell_idx} = ax1; % Store axis handle
    hold on;
    
    % Get STA data
    sta_time = cell_info.cell_res.near_spec{cell_info.cell_index}.sta_time;
    lfr_sta = cell_info.cell_res.near_lfr_spec{cell_info.cell_index}.sta_vals;
    hfr_sta = cell_info.cell_res.near_hfr_spec{cell_info.cell_index}.sta_vals;
    all_sta = cell_info.cell_res.near_spec{cell_info.cell_index}.sta_vals;
    
    % Plot STA traces
    plot(ax1, sta_time, lfr_sta, 'Color', c1, 'LineWidth', line_width);
    plot(ax1, sta_time, hfr_sta, 'Color', c2, 'LineWidth', line_width);
    plot(ax1, sta_time, all_sta, 'Color', c3, 'LineWidth', line_width);
    
    % Format STA plot
    ax1.Box = 'off';
    ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
    ax1.Title.String = sprintf('STA');
    ax1.XLabel.String = 'Time (sec)';
    ax1.Title.FontSize = 16;
    ax1.XAxis.FontSize = tick_font_size;
    ax1.YAxis.FontSize = tick_font_size;
    ax1.TickDir = 'out';
    ax1.YAxis.Exponent = 0; 
    
    % Plot PPC
    ax2 = subplot(2, n_cells, n_cells + cell_idx);
    ppc_axes{cell_idx} = ax2; % Store axis handle
    hold on;
    
    % Get PPC data based on cell type
    if strcmp(cell_info.cell_type, 'FSI')
        lfr_ppc = cell_info.cell_res.near_lfr_spec{cell_info.cell_index}.subsampled_ppc;
        hfr_ppc = cell_info.cell_res.near_hfr_spec{cell_info.cell_index}.subsampled_ppc;
        all_ppc = cell_info.cell_res.near_spec{cell_info.cell_index}.subsampled_ppc;
        freqs = cell_info.cell_res.near_spec{cell_info.cell_index}.freqs;
    else % MSN
        lfr_ppc = cell_info.cell_res.near_lfr_spec{cell_info.cell_index}.ppc';
        hfr_ppc = cell_info.cell_res.near_hfr_spec{cell_info.cell_index}.ppc';
        all_ppc = cell_info.cell_res.near_spec{cell_info.cell_index}.ppc';
        freqs = cell_info.cell_res.near_spec{cell_info.cell_index}.freqs;
    end
    
    % Plot PPC traces
    plot(ax2, freqs, lfr_ppc, 'Color', c1, 'LineWidth', line_width);
    plot(ax2, freqs, hfr_ppc, 'Color', c2, 'LineWidth', line_width);
    plot(ax2, freqs, all_ppc, 'Color', c3, 'LineWidth', line_width);
    
    % Find and plot PPC peaks for HFR and LFR trials
    % Find LFR peak
    [~, lfr_peak_idx] = max(lfr_ppc);
    lfr_peak_freq = freqs(lfr_peak_idx);
    xline(ax2, lfr_peak_freq, 'Color', c1, 'LineStyle', '--', 'LineWidth', line_width);
    
    % Find HFR peak
    [~, hfr_peak_idx] = max(hfr_ppc);
    hfr_peak_freq = freqs(hfr_peak_idx);
    xline(ax2, hfr_peak_freq, 'Color', c2, 'LineStyle', '--', 'LineWidth', line_width);
    
    % Format PPC plot
    ax2.Box = 'off';
    ax2.XLim = [2 100]; % Focus on frequencies above 2 Hz
    ax2.XTick = [2 10 25 50 75 100];
    ax2.Title.String = 'PPC';
    ax2.XLabel.String = 'Frequency (Hz)';
    ax2.YLabel.String = 'PPC';
    ax2.Title.FontSize = 16;
    ax2.XAxis.FontSize = tick_font_size;
    ax2.YAxis.FontSize = tick_font_size;
    ax2.TickDir = 'out';
    ax2.YAxis.Exponent = 0;
    
    % Add legend for PPC (only for first subplot)
    if cell_idx == 1
        legend(ax2, {'LFR', 'HFR', 'All Trials'}, ...
            'Location', 'best', 'FontSize', 12);
    end
end

%% Set Y-limits manually (uncomment and modify as needed)
% Example:
% % Set STA Y-limits (top row)
% for i = 1:n_cells
%     sta_axes(i).YLim = [-0.1 0.1]; % Adjust these values
% end
% 
% % Set PPC Y-limits (bottom row)
% for i = 1:n_cells
%     ppc_axes{i}.YLim = [-0.02 0.05]; % Adjust these values
%     ppc_axes{i}.YTick = [0:0.01:0.05]; % Set Y-tick labels
% end
for i = 1:n_cells
    sta_axes{i}.YTick = []; % Adjust these values
end


%% Save figure
save_filename = 'example_cells_combined.png';
print(fig, '-dpng', '-r300', fullfile(output_dir, save_filename));

fprintf('\nGenerated combined plot with %d cells\n', n_cells);
fprintf('Use sta_axes and ppc_axes arrays to set Y-limits manually\n');

