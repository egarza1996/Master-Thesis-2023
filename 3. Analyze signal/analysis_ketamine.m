clear;
clc;
close all;

% Code for ketamine (condition 4 and 5)

%% Add paths and confirm the run that you want to analyze - only things you need to modify

% Define the root path and the organoid folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/LOOP/exp1-2_11-25_07_23/'; % contains the results file we want to analyze
organoid_folders = arrayfun(@(x) sprintf('organoid%d', x), 1:9, 'UniformOutput', false);

addpath(genpath('/home/silviu/MATLAB Add-Ons/Collections/suplabel')); % contains suplabel function
addpath(genpath('/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Spikes-and-Deconvolution/loop')); 

% Set SAVE_FLAG to 1 if you want to save variables, otherwise set to 0
SAVE_FLAG = 1;

% Define saving path
base_save_directory = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Spikes-and-Deconvolution/loop';

% Loop through each organoid folder
for org_idx = 1:length(organoid_folders) % org3 is faulty
    if org_idx == 7 || org_idx == 8 % these two organoid recordings are very bad
        continue;
    end
    tic; % check for elapsed time of running the code
    current_folder = organoid_folders{org_idx};
    
    % Update the paths for the current organoid folder and "run_4" or "run_5" data file
    data_path = fullfile(root_path, current_folder, 'run_4');
    addpath(genpath(data_path)); % contains 'analysis_results.hdf5' file

    % Create the subfolder for the current organoid and run
    current_save_directory = fullfile(base_save_directory, current_folder, 'run_4');
    
    % Check if the folder exists, if not, create it
    if ~exist(current_save_directory, 'dir')
        mkdir(current_save_directory);
    end

%% Using the following signals: C, F_dff_dec, and S_dff

C = h5read('analysis_results.hdf5', '/estimates/C'); % raw fluorescence signal 
h5disp('analysis_results.hdf5', '/estimates/C');

F_dff_dec = h5read('analysis_results.hdf5', '/estimates/F_dff_dec'); % deconvolved fluorescence signal, representing an estimate of the underlying neural activity
h5disp('analysis_results.hdf5', '/estimates/F_dff_dec');

S_dff = h5read('analysis_results.hdf5', '/estimates/S_dff'); % inferred spiking activity based on the fluorescence signal
h5disp('analysis_results.hdf5', '/estimates/S_dff');


%% Select 3 traces for visualization (for random selection of traces see section below)

% Specify the neurons for visualization (can be any number of traces)
% selected_indices = [30, 50, 70]; % Example: [111, 126, 153, 157, 160]
% 
% num_traces = length(selected_indices);

% OR: Randomly select 3 neurons for visualization
num_neurons = size(S_dff, 2);  % Assuming S_dff has columns as traces/neurons
selected_indices = randperm(num_neurons, 3);  % Randomly select 3 traces

num_traces = length(selected_indices);

%% Time Series Plots for S_dff and F_dff_dec
% Helps in visualizing the temporal dynamics of the signals.

figure;

for i = 1:num_traces
    idx = selected_indices(i);
    
    % S_dff plot
    subplot(2, num_traces, i);
    plot(S_dff(:, idx), 'b');
    title(['Neuron #', num2str(idx)]);
    ylabel('Activity');
    xlabel('Frame Number');
    xlim([0 12600]); % ketamine (12600)
    
    % Add title for the top row
    if i == ceil(num_traces/2)
        text(0.5, 1.2, 'Inferred Spiking Activity for Selected Neurons', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    end
end

for i = 1:num_traces
    idx = selected_indices(i);
    
    % F_dff_dec plot
    subplot(2, num_traces, i + num_traces);
    plot(F_dff_dec(:, idx), 'r');
    title(['Neuron #', num2str(idx)]);
    ylabel('Activity');
    xlabel('Frame Number');
    xlim([0 12600]); % ketamine (12600)
    
    % Add title for the bottom row
    if i == ceil(num_traces/2)
        text(0.5, 1.2, 'Deconvolved Fluorescence Signal for Selected Neurons', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    end
end

%% Overlay Plots for C, S_dff, and F_dff_dec
% This provides a comparative visualization to see how the deconvolved signal and inferred spiking activity relate to the original fluorescence signal.

figure;

% Add a main title to the figure
sgtitle('Comparative Visualization of C, S_{dff} , and F_{dff-dec}');

for i = 1:num_traces
    idx = selected_indices(i);
    
    % Raw fluorescence
    subplot(num_traces, 1, i);
    
    % Left y-axis for Raw Fluorescence (C)
    yyaxis left;
    plot(1:size(C, 1), C(:, idx), 'g');
    ylabel('RFS');
    ylim([min(C(:, idx)) max(C(:, idx))]);
    xlim([0 12600]); % CHANGE depending on cond1 (6600) or cond2/3 (12600)
    set(gca, 'ycolor', [0 1 0]);  % Set y-axis color to green
    
    % Right y-axis for S_dff and F_dff_dec
    yyaxis right;
    plot(1:size(C, 1), S_dff(:, idx), 'b');
    hold on;
    plot(1:size(C, 1), F_dff_dec(:, idx), 'r');
    ylabel('ISA / DFS');
    ylim([min(S_dff(:, idx)) max(S_dff(:, idx))]); % Adjusting the y-axis limits based on S_dff data
    xlim([0 12600]); % CHANGE depending on cond1 (6600) or cond2/3 (12600)
    
    title(['Neuron #', num2str(idx)]);
    legend('Raw Fluorescence Signal (RFS)', 'Inferred Spiking Activity (ISA)', 'Deconvolved Fluorescence Signal (DFS)');
    xlabel('Frame Number');
end

%% Raster Plot for all neurons
% Quick overview of the spiking patterns of multiple neurons. Particularly for Sdff, where each row represents a neuron, and the dots indicate the times of inferred spikes.

% Extract number of frames and neurons
[num_frames, num_neurons] = size(S_dff);
frames_of_interest = 1:12600;  % For cond1 (6600) and for cond2/3 (12600)

% Define ketamine administration point (at 1 minute)
ketamine_administration_point = 1200; % at 60 seconds, given 20 frames per second

% For raster plots and other visualizations
figure;
hold on;

% Iterate over each neuron to plot spikes
for neuron = 1:num_neurons
    spike_frames = find(S_dff(:, neuron)); % Find frames with spikes for the current neuron
    plot(spike_frames, repmat(neuron, length(spike_frames), 1), '.', 'MarkerSize', 10);
end

% Add vertical line for ketamine administration
line([ketamine_administration_point ketamine_administration_point], [0, num_neurons + 1], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1.5); % Red dashed line for clarity

xlim([0, frames_of_interest(end)]);
ylim([0, num_neurons + 1]);
xlabel('Frame Number');
ylabel('Neuron Number');
title('Raster Plot of Inferred Spiking Activity for All Neurons');


%% Raster Plot for selected neurons

num_neurons = size(S_dff, 2);  % Total number of neurons

frames_of_interest = 1:12600;  % Adjust based on your data

% Define ketamine administration point (at 1 minute)
ketamine_administration_point = 1200; % at 60 seconds, given 20 frames per second

% For raster plots and other visualizations
figure;
hold on;

% Iterate over the selected neurons to plot spikes
for idx = 1:num_traces
    neuron = selected_indices(idx);
    spike_frames = find(S_dff(:, neuron)); % Find frames with spikes for the current neuron
    plot(spike_frames, repmat(idx, length(spike_frames), 1), '.', 'MarkerSize', 10);
end

% Add vertical line for ketamine administration
line([ketamine_administration_point ketamine_administration_point], [0, num_traces + 1], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1.5); % Red dashed line

xlim([0, 12600]);
ylim([0, num_traces + 1]);
xlabel('Frame Number');
ylabel('Selected Neuron');
title('Raster Plot of Inferred Spiking Activity for Selected Neurons');
set(gca, 'YTick', 1:num_traces, 'YTickLabel', arrayfun(@num2str, selected_indices, 'UniformOutput', false)); % Label y-axis with neuron IDs

%% Event Rate Plot for All Neurons after Ketamine Administration
% Calculate and plot the event rate (spike rate) over time for the inferred spiking activity of all neurons in 1-minute windows after ketamine administration.

% % Define constants based on dataset details
% frame_rate = 20; % 20 frames per second
% window_length = frame_rate * 60; % 60 seconds in frames (1 minute)
% total_frames = 12600; % this remains unchanged, total frames of interest
% 
% % Define ketamine administration point (at 1 minute)
% ketamine_administration_point = 1200; % at 60 seconds, given 20 frames per second
% 
% % Define 10 windows of 1 minute each after the ketamine administration
% all_windows = arrayfun(@(s) s:min(s+window_length-1, total_frames), ketamine_administration_point:window_length:ketamine_administration_point+10*window_length, 'UniformOutput', false);
% 
% % Number of neurons (all of them)
% num_neurons = size(S_dff, 2);
% 
% % Calculate and plot event rate for each neuron
% figure;
% for neuron = 1:num_neurons
%     subplot(num_neurons, 1, neuron);
%     hold on;
% 
%     for w = 1:length(all_windows)
%         window_frames = all_windows{w};
% 
%         % Calculate event rate for the current window
%         event_rate = sum(S_dff(window_frames, neuron)) / length(window_frames);
% 
%         % Plotting - the x-axis will be the central frame of the window
%         central_frame = mean(window_frames);
%         plot(central_frame, event_rate, 'o');
%     end
% 
%     % Add vertical line for ketamine administration
%     line([ketamine_administration_point ketamine_administration_point], ylim, 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1.5); % Red dashed line for clarity
% 
%     xlabel('Frame Number');
%     ylabel('Event Rate');
%     title(['Neuron ', num2str(neuron)]);
%     xlim([0, total_frames]);
% end

%% Event Rate Plot for Selected Neurons after Ketamine Administration
% Calculate and plot the event rate (spike rate) over time for the inferred spiking activity of selected neurons in 1-minute windows after ketamine administration.

% Define constants based on dataset details
frame_rate = 20; % 20 frames per second
window_length = frame_rate * 60; % 60 seconds in frames (1 minute)
total_frames = 12600; % this remains unchanged, total frames of interest

% Define ketamine administration point (at 1 minute)
ketamine_administration_point = 1200; % at 60 seconds, given 20 frames per second

% Define 10 windows of 1 minute each after the ketamine administration
all_windows = arrayfun(@(s) s:min(s+window_length-1, total_frames), ketamine_administration_point:window_length:ketamine_administration_point+10*window_length, 'UniformOutput', false);

% Calculate and plot event rate for each selected neuron
figure;

% Add a main title to the figure
sgtitle('Event Rate Plot for Selected Neurons');

for idx = 1:num_traces
    neuron = selected_indices(idx);
    subplot(num_traces, 1, idx);
    hold on;
    
    for w = 1:length(all_windows)
        window_frames = all_windows{w};
        
        % Calculate event rate for the current window
        event_rate = sum(S_dff(window_frames, neuron)) / length(window_frames);
        
        % Plotting - the x-axis will be the central frame of the window
        central_frame = mean(window_frames);
        plot(central_frame, event_rate, 'o');
    end
    
    % Add vertical line for ketamine administration
    line([ketamine_administration_point ketamine_administration_point], ylim, 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1.5); % Red dashed line for clarity
    
    xlabel('Frame Number');
    ylabel('Event Rate');
    title(['Neuron ', num2str(neuron)]);
    xlim([0, total_frames]);
end

%% Alternative view: Event Rate Plot for selected neurons (with gray dotted lines every 30 seconds)
    % (Bar Representation) This view provides a continuous representation of event rates over the entire duration of interest, giving the appearance of lines or bars.
    % (Dot Representation) Wheras the previous section provides a discrete representation of event rates for specific windows, resulting in dots plotted at the central frame of each window.

% Define constants based on dataset details
frame_rate = 20; % 20 frames per second
window_length = frame_rate * 60; % 60 seconds in frames (1 minute)
ketamine_administration_point = 1200; % at 60 seconds, given 20 frames per second

% Define the last start point for a full window within total_frames
last_start_point = total_frames - window_length;

% Define 10 windows of 1 minute each after the ketamine administration
% Ensure that we do not exceed the last_start_point
all_windows = arrayfun(@(s) s:min(s+window_length-1, total_frames), ketamine_administration_point:window_length:last_start_point, 'UniformOutput', false);

% Calculate and plot event rate for each selected neuron
figure;

% Add a main title to the figure
sgtitle('Event Rate Plot for Selected Neurons');

for idx = 1:num_traces
    neuron = selected_indices(idx);
    subplot(num_traces, 1, idx);
    hold on;
        
    event_rate = zeros(total_frames, 1); % Initialize the event rate array with zeros
    for w = 1:length(all_windows)
        window_frames = all_windows{w};
        start_frame = window_frames(1);
        end_frame = window_frames(end);

        % Calculate the event rate for the current window
        event_rate(start_frame:end_frame) = sum(S_dff(start_frame:end_frame, neuron)) / window_length;

    end
    
    % Plotting the event rate over time as a continuous line
    plot(event_rate);
    
    % Add gray dotted vertical lines every 60 seconds
    for frame = window_length:window_length:total_frames
        line([frame frame], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1.5);
    end

    % Add vertical line for ketamine administration
    line([ketamine_administration_point ketamine_administration_point], ylim, 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1.5); % Red dashed line for clarity

    xlabel('Frame Number');
    ylabel('Event Rate');
    title(['Neuron ', num2str(neuron)]);
    xlim([0, total_frames]);
end

%% Saving variables and figures

if SAVE_FLAG
    % Gather all variables from the workspace into a structure
    allVariables = whos;
    workspaceData = struct();
    
    for idx = 1:length(allVariables)
        varName = allVariables(idx).name;
        workspaceData.(varName) = eval(varName);
    end
    
    % Save the variables to the specified file
    variable_file_name = 'spikes_deconvolution_data.mat';
    full_path = fullfile(current_save_directory, variable_file_name);
    save(full_path, 'workspaceData');

    % Save all open figures
    % all_figs = findall(0, 'Type', 'figure');
    % for idx = 1:length(all_figs)
    %     fig = all_figs(idx);
    %     fig_file_name = sprintf('figure_%d.png', fig.Number); % You can adjust the format, e.g., .fig, .jpg, etc.
    %     full_fig_path = fullfile(current_save_directory, fig_file_name);
    %     saveas(fig, full_fig_path);
    % end

    % Save all open figures
    all_figs = findall(0, 'Type', 'figure');
    for idx = 1:length(all_figs)
        fig = all_figs(idx);
        set(fig, 'WindowState', 'maximized'); % Maximize the figure window
        drawnow; % Update figure window
    
        fig_file_name = sprintf('figure_%d.png', fig.Number); % You can adjust the format, e.g., .fig, .jpg, etc.
        full_fig_path = fullfile(current_save_directory, fig_file_name);
        saveas(fig, full_fig_path);
    
        % set(fig, 'WindowState', 'normal'); % Optional: Return the figure window to its original state
    end

end


%% Loop ends and starts over
    % At the end of the loop iteration, remove the added path to keep the workspace clean
    rmpath(genpath(data_path));
    close all;
end
