clear;
clc;
close all;

% After running the code "PSD_analysis_updated.m" this code is useful for
% plotting each organoid's Overall Average PSD Across Frequency Bands
% Note: in this code, you manually choose which organoid to analyze.


%% Add paths 

% Define the root path and the organoid and run folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/Power Spectral Density/loop'; % contains "PSD_results.mat" which contains "overall_avg"
organoid_folders = arrayfun(@(x) sprintf('organoid%d', x), 1:9, 'UniformOutput', false);
addpath(genpath('/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/Power Spectral Density/loop'));

SAVE_FLAG = 1;
base_save_directory = fullfile(root_path, 'Power Spectral Density', 'loop');

% Data structure to store the averages
overall_avgs = struct();


% MANUALLY CHOOSE AN ORGANOID
% AND adjust " run_folders" and "custom_x_labels" for plotting
chosen_organoid = 1;
run_folders = {'run', 'run_1', 'run_2', 'run_3', 'run_4', 'run_5'}; % all runs
% run_folders = {'run', 'run_1', 'run_2', 'run_4', 'run_5'}; % for org4 and org9
% run_folders = {'run', 'run_1', 'run_2', 'run_5'}; % for org6


% Extract overall_avg_1, overall_avg_2, and overall_avg_3 values.
for org_idx = chosen_organoid % 1:length(organoid_folders)
    if org_idx == 3 || org_idx == 7 || org_idx == 8 % Skip bad recordings
        continue;
    end

    organoid_folder = organoid_folders{org_idx};

    for run_idx = 1:length(run_folders)
        run_folder = run_folders{run_idx};
        
        % Update the paths for the current organoid folder and run data file
        data_path = fullfile(root_path, organoid_folder, run_folder);
        addpath(genpath(data_path));

        % Print the constructed file path for debugging
        data_file = fullfile(data_path, 'PSD_results.mat');
        disp(data_file); % Print the path

        % Check if the file exists before attempting to load it
        if isfile(data_file)
            data = load(data_file, 'overall_avg_1', 'overall_avg_2', 'overall_avg_3');
            disp(data); % Print the loaded data for verification
            overall_avgs(org_idx).(run_folder) = [data.overall_avg_1, data.overall_avg_2, data.overall_avg_3];
        else
            warning('File not found: %s\n', data_file);
        end

    end
end

%% Saving "overall_avg" as a .CSV file

if SAVE_FLAG
    % Flatten the data structure for CSV export
    csv_data = [];
    row_names = {};  % To keep track of row labels   


    for org_idx = chosen_organoid 
        for run_idx = 1:length(run_folders)
            run_folder = run_folders{run_idx};
            
            if isfield(overall_avgs(org_idx), run_folder)
                row = overall_avgs(org_idx).(run_folder);
            else
                row = [NaN, NaN, NaN];  % Use NaN for missing data
            end
    
            csv_data = [csv_data; row];  % Append the row to the data
            row_names{end+1} = run_folder;  % Keep track of the run name
        end
    end
    
    % Convert to table for easy CSV writing
    csv_table = array2table(csv_data, 'VariableNames', {'Overall_Avg_1_(0.1-0.4Hz)', 'Overall_Avg_2_(0.5-4Hz)', 'Overall_Avg_3_(5-8Hz)'});
    csv_table.Run = row_names';  % Add run names as a new column
    
    % Move 'Run' column to the first position
    csv_table = csv_table(:, [end, 1:end-1]);
    
    % Define the base directory and create "save_results" folder
    base_save_directory = fullfile(root_path, 'save_results');
    
    % Check if the directory exists, if not, create it
    if ~exist(base_save_directory, 'dir')
        mkdir(base_save_directory);  % Create the directory
    end
    
    % Generate the filename including the org_idx value
    filename = ['overall_avgs_org_' num2str(org_idx) '.csv'];
    
    % Define the CSV file path with the dynamic filename
    csv_file_path = fullfile(base_save_directory, filename);
    
    % Write to CSV
    writetable(csv_table, csv_file_path);
    
    disp(['Data saved to CSV file: ' csv_file_path]);
end


%% Plot Overall Average PSD Across Runs (GROUPED BY CONDITIONS)

org_idx = chosen_organoid;

% Define your custom labels here
custom_x_labels = {'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Opto. 1 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL'};
% custom_x_labels = {'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL'}; % for org4 and org9
% custom_x_labels = {'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 200 \muL'}; % for org6
% Ensure that the number of custom labels matches the number of runs

% Prepare data for plotting
num_runs = length(run_folders);
num_avgs = 3; % Since there are 3 overall_avg values
colors = ['b', 'r', 'g']; % Define colors for each set of averages

% Initialize the data for the grouped bar plot
grouped_plot_data = zeros(num_runs, num_avgs);

% Loop to populate the data matrix
for i = 1:num_runs
    for j = 1:num_avgs
        % Check if data exists for the run, otherwise use NaN
        if isfield(overall_avgs(chosen_organoid), run_folders{i})
            grouped_plot_data(i, j) = overall_avgs(chosen_organoid).(run_folders{i})(j);
        else
            grouped_plot_data(i, j) = NaN;
        end
    end
end

% Plotting with logarithmic y-axis
fig = figure;
set(fig, 'Position', [100, 100, 800, 500]); % Position and size: [left, bottom, width, height] in pixels
hold on;
bar_handle = bar(grouped_plot_data, 'grouped', 'EdgeColor', 'none');
set(bar_handle, {'FaceColor'}, num2cell(colors, 1)'); % Set colors for each avg


% Set x-axis labels for each run
set(gca, 'xtick', 1:num_runs, 'xticklabel', custom_x_labels);
% set(gca, 'xtick', 1:num_runs, 'xticklabel', run_folders); % Set x-axis labels for each run
set(gca, 'YScale', 'log'); % Set the y-axis to logarithmic scale
xtickangle(45);
ylabel('Power (Log Scale)');
title(sprintf('Organoid %d \n Overall Average PSD Across Frequency Bands', org_idx));

% Add legend to clarify which color corresponds to which overall_avg
legend(bar_handle, {'0.1 - 0.4 Hz', '0.5 - 4 Hz', '5 - 8 Hz'});

% Save the plot if required
if SAVE_FLAG
    save_file = fullfile(base_save_directory, sprintf('organoid%d_overall_avg_PSD_1.png', org_idx));
    saveas(fig, save_file);
end


%% Plot Overall Average PSD Across Runs (GROUPED BY FREQUENCIES) - customized labels for x axis

org_idx = chosen_organoid;

% Define your custom x-axis labels
custom_x_labels = {'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Opto. 1 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL', ...
                   'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Opto. 1 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL', ...
                   'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Opto. 1 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL'}; 

% custom_x_labels = {'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL', ...
%                    'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL', ...
%                    'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 100 \muL', 'Ket. 200 \muL'}; % for org4 and org9

% custom_x_labels = {'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 200 \muL', ...
%                    'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 200 \muL', ...
%                    'Spont. Act.', 'Opto. 1 ms, 10 s', 'Opto. 50 ms, 20s', 'Ket. 200 \muL'}; % for org6

% Ensure the number of labels matches the number of bars
if length(custom_x_labels) ~= num_runs * num_avgs
    error('The number of custom x-axis labels does not match the number of bars.');
end

% Define your custom labels for the legend
custom_legend_labels = {'0.1 - 0.4 Hz', '0.5 - 4 Hz', '5 - 8 Hz'};
% Ensure that the number of custom labels matches the number of averages

% Prepare data for plotting
num_runs = length(run_folders);
num_avgs = 3; % Since there are 3 overall_avg values
colors = ['b', 'r', 'g']; % Define colors for each set of averages

% Initialize plot_data and x_labels
plot_data = zeros(num_avgs, num_runs); % Data for each group
x_labels = cell(1, num_runs * num_avgs);

% Loop to aggregate data and labels
for j = 1:num_avgs
    for i = 1:num_runs
        % Check if data exists for the run, otherwise use NaN
        if isfield(overall_avgs(chosen_organoid), run_folders{i})
            plot_data(j, i) = overall_avgs(chosen_organoid).(run_folders{i})(j);
        else
            plot_data(j, i) = NaN;
        end

        % Create label for each bar
        % x_labels{(j-1)*num_runs + i} = sprintf('%s (avg_%d)', run_folders{i}, j);
    end
end

% Plotting with logarithmic y-axis
fig = figure;
set(fig, 'Position', [100, 100, 800, 500]); % Position and size: [left, bottom, width, height] in pixels
hold on;

% Plot each group with a different color
for j = 1:num_avgs
    bar_handle(j) = bar(((1:num_runs) + (j-1)*num_runs), plot_data(j, :), colors(j));
end

% Set x-axis and other plot properties
ax = gca;
ax.XTick = 1:(num_runs*num_avgs); % Set a tick for each bar

% Set your custom labels for the x-axis
ax.XTickLabel = custom_x_labels; % Use your custom labels

set(ax, 'YScale', 'log'); % Set the y-axis to logarithmic scale
xtickangle(45);
ylabel('Power (Log Scale)');
title(sprintf('Organoid %d \n Overall Average PSD Across Frequency Bands', org_idx));

% Create the legend
legend(bar_handle, custom_legend_labels, 'Location', 'best');

hold off; % Release the hold on the current axes

% Save the plot if required
if SAVE_FLAG
    save_file = fullfile(base_save_directory, sprintf('organoid%d_overall_avg_PSD_2.png', org_idx));
    saveas(fig, save_file);
end

close all;

%% Plot Overall Average PSD Across Runs (GROUPED BY FREQUENCIES) - non customized labels for x axis

% % Define your custom labels for the legend
% custom_legend_labels = {'0.1 - 0.4 Hz', '0.5 - 4 Hz', '5 - 8 Hz'};
% % Ensure that the number of custom labels matches the number of averages
% 
% % Prepare data for plotting
% num_runs = length(run_folders);
% num_avgs = 3; % Since there are 3 overall_avg values
% colors = ['b', 'r', 'g']; % Define colors for each set of averages
% 
% % Initialize plot_data and x_labels
% plot_data = zeros(num_avgs, num_runs); % Data for each group
% x_labels = cell(1, num_runs * num_avgs);
% 
% % Loop to aggregate data and labels
% for j = 1:num_avgs
%     for i = 1:num_runs
%         % Check if data exists for the run, otherwise use NaN
%         if isfield(overall_avgs(chosen_organoid), run_folders{i})
%             plot_data(j, i) = overall_avgs(chosen_organoid).(run_folders{i})(j);
%         else
%             plot_data(j, i) = NaN;
%         end
% 
%         % Create label for each bar
%         x_labels{(j-1)*num_runs + i} = sprintf('%s (avg_%d)', run_folders{i}, j);
%     end
% end
% 
% % Plotting with logarithmic y-axis
% fig = figure;
% set(fig, 'Position', [100, 100, 800, 500]); % Position and size: [left, bottom, width, height] in pixels
% hold on;
% 
% % Plot each group with a different color
% for j = 1:num_avgs
%     bar_handle(j) = bar(((1:num_runs) + (j-1)*num_runs), plot_data(j, :), colors(j));
% end
% 
% % Set x-axis and other plot properties
% ax = gca;
% ax.XTick = 1:(num_runs*num_avgs); % Set a tick for each bar
% ax.XTickLabel = x_labels; % Set labels for each tick
% set(ax, 'YScale', 'log'); % Set the y-axis to logarithmic scale
% xtickangle(45);
% ylabel('Power (Log Scale)');
% title(sprintf('Organoid %d \n Overall Average PSD Across Frequency Bands', org_idx));
% 
% % Create the legend
% legend(bar_handle, custom_legend_labels, 'Location', 'best');
% 
% hold off; % Release the hold on the current axes
% 
% % Save the plot if required
% if SAVE_FLAG
%     save_file = fullfile(base_save_directory, sprintf('organoid%d_overall_avg_PSD_2.png', org_idx));
%     saveas(fig, save_file);
% end


%% Prepare data for plotting - outdated plot, the other ones are better

% plot_data = [];
% x_labels = {};
% 
% % Loop through each overall_avg
% for j = 1:3
%     % Loop through each run
%     for i = 1:length(run_folders)
%         avg_label = sprintf('avg_%d', j);
% 
%         % Check if data exists for the run, otherwise use NaN
%         if isfield(overall_avgs(4), run_folders{i})  % Assuming data for Organoid 4 is in overall_avgs(4)
%             plot_data(end+1) = overall_avgs(4).(run_folders{i})(j);
%         else
%             plot_data(end+1) = NaN;
%         end
% 
%         x_labels{end+1} = sprintf('%s (%s)', run_folders{i}, avg_label);
%     end
% end
% 
% % Plotting with logarithmic y-axis
% 
% bar_handle = bar(plot_data, 'grouped');
% set(gca, 'xticklabel', x_labels);
% set(gca, 'YScale', 'log'); % Set the y-axis to logarithmic scale
% xtickangle(45);
% ylabel('Power (Log Scale)');
% title('Overall Average PSD for Organoid 4 Across Runs'); 
