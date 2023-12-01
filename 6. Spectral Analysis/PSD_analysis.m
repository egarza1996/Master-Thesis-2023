clear;
clc;
close all;

% Calculate the Power Spectral Density (PSD) for all organoids using their complete deconvolved fluorescence signal.
% This code used the function "calc_PSD.m".

%% Add paths 

% Define the root path and the organoid and run folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/loop-new'; % contains "clustering_data_new.mat" which contains "F_dff_dec"
organoid_folders = arrayfun(@(x) sprintf('organoid%d', x), 1:9, 'UniformOutput', false);
run_folders = {'run', 'run_1', 'run_2', 'run_3', 'run_4', 'run_5'}; 
addpath(genpath('/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/Power Spectral Density/loop'));


SAVE_FLAG = 1;
base_save_directory = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/Power Spectral Density/loop';


% Parameters for calc_PSD function
NFFT = 200;  % Number of frequency bins
    % We chose 200 because we want our freq. res. to be 0.1 Hz and since our
    % frame rate is 20, that means 20/0.1 = 200
sr = 20;     % Sample rate (Hz)
window = 5;  % Window length in seconds (e.g. 2 seconds)
overlap = 2.5; % Overlap in seconds (e.g. 1 sec would be 50% of the window size, if window size is 2 sec)


% Loop through each organoid and run folder
for org_idx = 1:length(organoid_folders)
    if org_idx == 3 || org_idx == 7 || org_idx == 8 % Skip bad recordings
        continue;
    end
    
    tic; % Start timer
    organoid_folder = organoid_folders{org_idx};
    
    for run_idx = 1:length(run_folders)
        run_folder = run_folders{run_idx};

        % Update the paths for the current organoid folder and run data file
        data_path = fullfile(root_path, organoid_folder, run_folder);
        addpath(genpath(data_path));

        % Print the constructed file path for debugging
        data_file = fullfile(data_path, 'clustering_data_new.mat');
        disp(data_file); % Print the path

        % Check if the file exists before attempting to load it
        if isfile(data_file)
            load(data_file, 'F_dff_dec'); % Load the file if it exists

            % Get the number of rows in F_dff_dec
            num_rows = size(F_dff_dec, 1);

            % Slicing logic based on run_idx
            if strcmp(run_folder, 'run')
                % Slice for "run" if the number of rows is sufficient
                F_dff_dec = F_dff_dec(1:min(num_rows, 6599), :); % 5min 30sec
            elseif strcmp(run_folder, 'run_1')
                % Slice for "run_1" if the number of rows is sufficient
                F_dff_dec = F_dff_dec(600:min(num_rows, 6599), :); % 5min (discard first 30sec of baseline activity)
            elseif strcmp(run_folder, 'run_2')
                % Slice for "run_2" if the number of rows is sufficient
                F_dff_dec = F_dff_dec(600:min(num_rows, 12599), :); % 10min (discard first 30sec of baseline activity)
            elseif strcmp(run_folder, 'run_3')
                % Slice for "run_3" if the number of rows is sufficient
                F_dff_dec = F_dff_dec(600:min(num_rows, 12599), :); % 10min (discard first 30sec of baseline activity)
            elseif strcmp(run_folder, 'run_4')
                % Slice for "run_4" if the number of rows is sufficient
                F_dff_dec = F_dff_dec(1200:min(num_rows, 12599), :); % 9.5min (discard first 1min of baseline activity)            
            else
                % Slice for "run_5" if the number of rows is sufficient
                F_dff_dec = F_dff_dec(1200:min(num_rows, 12599), :); % 9.5min (discard first 1min of baseline activity)
            end

            % Create the subfolder for the current organoid and run
            current_save_directory = fullfile(base_save_directory, organoid_folder, run_folder);

            % Check and create the folder if it doesn't exist
            if ~exist(current_save_directory, 'dir')
                mkdir(current_save_directory);
            end


            %% Call calc_PSD function
            [DSA, frequencies] = calc_PSD(F_dff_dec, NFFT, sr, window, overlap);

            % DSA: A 3D matrix where each element represents the power spectral density (PSD) of a specific frequency, for a specific channel (neuron), during a specific time window.
            % frequencies: A vector indicating the frequencies corresponding to each element in the PSD.


            %% Plots

            % Plot Average PSD Across All neurons
            % To get an overall sense of the frequency content across all neurons
            avg_psd = mean(mean(DSA, 3), 2);  % Average across neurons and windows
            figure;
            plot(frequencies, avg_psd);
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('Average PSD Across Neurons');

            % Heatmap of PSD Across Neurons
            % Useful for identifying patterns across multiple neurons
            avg_psd_neurons = mean(DSA, 3);  % Average across windows
            figure;
            imagesc(avg_psd_neurons);
            xlabel('Neuron #');
            ylabel('Frequency Bin');
            title('PSD Across Neurons');
            colorbar;


            %% %% Average across frequency bands
            % 1) Identify the indices in the frequency vector that correspond to the desired range, e.g. 0.1 to 0.4 Hz.
            % 2) Extract the relevant portions of the DSA matrix that correspond to these frequencies, i.e. DSA(indices, :, :)
            % 3) Average the extracted data across the time dimension.

            %% For 0.1 to 0.4 Hz
            freq_min_1 = 0.1;
            freq_max_1 = 0.4;
            
            indices_1 = find(frequencies >= freq_min_1 & frequencies <= freq_max_1);
            average_dsa_1 = mean(DSA(indices_1, :, :), 3); % Averaging over the third dimension (time)
                % Represents the average PSD values across time for each frequency bin and neuron within the specified frequency range (0.1 to 0.4 Hz)
            
            % Average across frequency bins (rows of average_dsa)
            avg_across_freq_bins_1 = mean(average_dsa_1, 1); % Result is a row vector
            % Reshape to get a column vector
            avg_across_freq_bins_1 = avg_across_freq_bins_1(:);
                % This will give you the average PSD for each channel, averaged over both time and the specified frequency range.

            % Average across neurons (columns of average_dsa)
            avg_across_neurons_1 = mean(average_dsa_1, 2); % Result is a column vector
                 % This will give you the average PSD for each frequency bin, averaged over all neurons and time windows.
            
            % Overall average across all elements of average_dsa
            overall_avg_1 = mean(average_dsa_1(:));
                % This will give you a single number representing the average PSD across all frequency bins, neurons, and time windows in the specified frequency range.
            
            
            %% For 0.5 to 4 Hz
            freq_min_2 = 0.5;
            freq_max_2 = 4;
            
            indices_2 = find(frequencies >= freq_min_2 & frequencies <= freq_max_2);
            average_dsa_2 = mean(DSA(indices_2, :, :), 3); % Averaging over the third dimension (time)
                % Represents the average PSD values across time for each frequency bin and channel within the specified frequency range (0.1 to 0.4 Hz)
               
            % Average across frequency bins (rows of average_dsa)
            avg_across_freq_bins_2 = mean(average_dsa_2, 1); % Result is a row vector
            % Reshape to get a column vector
            avg_across_freq_bins_2 = avg_across_freq_bins_2(:);
                % This will give you the average PSD for each channel, averaged over both time and the specified frequency range.

            % Average across neurons (columns of average_dsa)
            avg_across_neurons_2 = mean(average_dsa_2, 2); % Result is a column vector
                 % This will give you the average PSD for each frequency bin, averaged over all neurons and time windows.
            
            % Overall average across all elements of average_dsa
            overall_avg_2 = mean(average_dsa_2(:));
                % This will give you a single number representing the average PSD across all frequency bins, neurons, and time windows in the specified frequency range.
            
            %% For 5 to 8 Hz
            freq_min_3 = 5;
            freq_max_3 = 8;
            
            indices_3 = find(frequencies >= freq_min_3 & frequencies <= freq_max_3);
            average_dsa_3 = mean(DSA(indices_3, :, :), 3); % Averaging over the third dimension (time)
                % Represents the average PSD values across time for each frequency bin and channel within the specified frequency range (0.1 to 0.4 Hz)
            
            % Average across frequency bins (rows of average_dsa)
            avg_across_freq_bins_3 = mean(average_dsa_3, 1); % Result is a row vector
            % Reshape to get a column vector
            avg_across_freq_bins_3 = avg_across_freq_bins_3(:);
                % This will give you the average PSD for each channel, averaged over both time and the specified frequency range.

            % Average across neurons (columns of average_dsa)
            avg_across_neurons_3 = mean(average_dsa_3, 2); % Result is a column vector
                 % This will give you the average PSD for each frequency bin, averaged over all neurons and time windows.
            
            % Overall average across all elements of average_dsa
            overall_avg_3 = mean(average_dsa_3(:));
                % This will give you a single number representing the average PSD across all frequency bins, neurons, and time windows in the specified frequency range.


            %% Plots

            % Compare the average PSD across frequency bins for different conditions or time periods.


            % Create a scatter plot with logarithmic y-axis
            fig = figure;
            set(fig, 'Position', [100, 100, 800, 500]); % Position and size: [left, bottom, width, height] in pixels
            hold on;
            hold on; % Allows multiple plots in the same figure
            
            % Plotting each dataset with different colors and markers
            scatter(1:length(avg_across_freq_bins_1), avg_across_freq_bins_1, 'filled', 'MarkerFaceColor', 'r');
            scatter(1:length(avg_across_freq_bins_2), avg_across_freq_bins_2, 'filled', 'MarkerFaceColor', 'g');
            scatter(1:length(avg_across_freq_bins_3), avg_across_freq_bins_3, 'filled', 'MarkerFaceColor', 'b');
            
            % Set the Y-axis to logarithmic scale
            set(gca, 'YScale', 'log');
            
            xlabel('Neurons');
            ylabel('Average PSD (Log Scale)');
            title('Comparison of Average PSD Across Frequency Bands');
            legend('0.1 - 0.4 Hz', '0.5 - 4 Hz', '5 - 8 Hz');
            grid on;
            hold off; % Release the figure for other plots

            %% Save results if SAVE_FLAG is set
            if SAVE_FLAG
                save_file = fullfile(current_save_directory, 'PSD_results.mat');
                save(save_file, 'DSA', 'frequencies', ...
                    "indices_1", "average_dsa_1", "avg_across_freq_bins_1", "avg_across_neurons_1", "overall_avg_1", ...
                    "indices_2", "average_dsa_2", "avg_across_freq_bins_2", "avg_across_neurons_2", "overall_avg_2", ...
                    "indices_3", "average_dsa_3", "avg_across_freq_bins_3", "avg_across_neurons_3", "overall_avg_3");

                % Save all open figures
                all_figs = findall(0, 'Type', 'figure');
                for idx = 1:length(all_figs)
                    fig = all_figs(idx);
                    fig_file_name = sprintf('figure_new_%d.png', fig.Number); % You can adjust the format, e.g., .fig, .jpg, etc.
                    full_fig_path = fullfile(current_save_directory, fig_file_name);
                    saveas(fig, full_fig_path);
                end
            end
            close all;

        else
            fprintf('File not found: %s\n', data_file);
            % Handle the case where the file does not exist
        end

    end
    elapsed_time = toc;
    fprintf('Elapsed time for processing organoid %d is %.6f seconds.\n', org_idx, elapsed_time);

    % At the end of the loop iteration, remove the added path to keep the workspace clean
    rmpath(genpath(data_path));

end
