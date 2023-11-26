clear;
clc;
close all;

% K-means algorithm that clusters neurons based on their temporal activity patterns. 
% Here, the K-means++ clustering is performed on the correlation matrix to categorize neurons based on their similarity in activity during the entire signal or specific windows of interest. 

% Code for condition 1 (1msec stimulation every 10sec)

%% Add paths and define F -- only things you need to modify

% Define the root path and the organoid folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/LOOP/exp1-2_11-25_07_23/';
organoid_folders = arrayfun(@(x) sprintf('organoid%d', x), 1:9, 'UniformOutput', false);
addpath(genpath('/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/signal-correlation-master')); % contains meta_k_means_cond1.m function
addpath(genpath('/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/loop-new')); 

% Set SAVE_FLAG to 1 if you want to save variables, otherwise set to 0
SAVE_FLAG = 1;

% Define saving path
base_save_directory = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Clustering-Correlation-Distance/loop-new';

% Loop through each organoid folder
for org_idx = 6 % 6:length(organoid_folders)
    if org_idx == 3 ||org_idx == 7 || org_idx == 8 % these organoid recordings are very bad
        continue;
    end
    tic; % check for elapsed time of running the code
    current_folder = organoid_folders{org_idx};
    
    % Update the paths for the current organoid folder and "run_1" data file
    data_path = fullfile(root_path, current_folder, 'run_1');
    addpath(genpath(data_path)); % contains 'analysis_results.hdf5' file

    % Create the subfolder for the current organoid and run
    current_save_directory = fullfile(base_save_directory, current_folder, 'run_1');
    
    % Check if the folder exists, if not, create it
    if ~exist(current_save_directory, 'dir')
        mkdir(current_save_directory);
    end

%% Using the deconvolved fluorescence signal (F_dff_dec) since it will likely provide a richer representation of the oscillatory dynamics in the post-stimulation windows.
% As opposed to the Inferred Spiking Activity (S_dff), which represents a binary or discrete estimate of when a neuron is firing.
analysis_results_file = fullfile(data_path, 'analysis_results.hdf5');
F_dff_dec = h5read(analysis_results_file, '/estimates/F_dff_dec');
h5disp(analysis_results_file, '/estimates/F_dff_dec');

%% Need to remove columns/neurons with a lot of zeroes before running k-means, otherwise it fails
    % (usually only about 2 to 6 columns get removed)
    % Needed to incorporate this in order for org3 to run (otherwise it did not run)... Might as well do it for all runs to keep it consistent

% Number of columns before removal
cols_before = size(F_dff_dec, 2);

% Threshold for consecutive zeros
threshold = 1200; % if 60sec of continuous inactivity (i.e. zeroes) are found, the whole column gets discarded
% 1200 frames / 20 fr/sec = 50sec

% Create a logical vector to track columns with consecutive zeros exceeding the threshold
remove_columns = false(1, cols_before);

for col = 1:cols_before
    % Identify start and stop of zero sequences
    zero_diff = diff([0; (F_dff_dec(:, col) == 0); 0]);
    starts = find(zero_diff == 1);
    stops = find(zero_diff == -1) - 1;

    % Check for consecutive zero sequences that exceed the threshold
    zero_lengths = stops - starts + 1;
    if any(zero_lengths >= threshold)
        remove_columns(col) = true;
    end
end

% Remove columns that exceed the threshold
F_dff_dec(:, remove_columns) = [];

% Display how many columns were removed
cols_after = size(F_dff_dec, 2);
fprintf('Removed %d columns. Remaining columns: %d\n', cols_before - cols_after, cols_after);


%% Initialize constants and define our windows of interest (post-stimulation windows)
frame_rate = 20; % 20 frames per second
num_neurons = size(F_dff_dec, 2); % includes ALL neurons
stimulation_start = 600; % stimulation starts at 30 seconds
stimulation_interval = 200; % CHANGE: for cond1 every 10 sec (200) and for cond2/3 every 20 sec (400)
total_frames = 6000; % CHANGE: For cond1 (6600) and for cond2/3 (12600)
% usnure why the duration of organoid 1 is only 5min so I will change total_frames to 6000 instead of the usual 6600.

% Define stimulation points within the range of interest
stimulation_points = stimulation_start:stimulation_interval:total_frames;

% Define windows of interest: 15-second window starting 5 seconds after each stimulation
% windows_15s = arrayfun(@(s) s+frame_rate*5:min(s+frame_rate*20-1, total_frames), stimulation_points, 'UniformOutput', false);

% Adapted for cond1:
% Define windows of interest: 5-second window starting 5 seconds after each stimulation
windows_5s = arrayfun(@(s) s+frame_rate*5:min(s+frame_rate*10-1, total_frames), stimulation_points, 'UniformOutput', false);

% Filter out empty windows
all_windows = windows_5s(~cellfun('isempty', windows_5s));

% Segment the F_dff_dec data using the defined windows of interest
segmented_data = cellfun(@(window) F_dff_dec(window, :), all_windows, 'UniformOutput', false);

% Combine these segments to create a new matrix with the desired format
F_dff_dec_segmented = vertcat(segmented_data{:});

%% Meta k means
[allclusters, centroidcorr, dendmem, dunnsinitial] = meta_k_means_cond1(F_dff_dec_segmented', 'correlation');

% Input:
    % eventsmat: each row is expected to be a a neuron and each column a frame.

% Return:
    % allclusters: A cell array containing the final clusters and their centroids.
    % centroidcorr: Correlation coefficients between cluster centroids.
    % dendmem: A matrix indicating graded membership of dendrites to clusters.
    % dunnsinitial: A measure of clustering validity, possibly the Dunn's index.

%% Elapsed time for meta_k_means
elapsed_time = toc;
fprintf('Elapsed time for meta_k_means calculation is %.6f seconds.\n', elapsed_time);

%% Fixing the neuron clustering repetition issue
% The meta_k_means function sometimes assigns one neuron to several
% clusters. We want the neuron to only be assigned to one cluster.

% Instead of running "meta_k_means" again, just load the data that we had previously calculated
% load('clustering_data.mat');

% Number of neurons and clusters
num_neurons = size(F_dff_dec_segmented, 2); % Columns represent neurons
num_clusters = size(allclusters, 1); % Number of clusters

% Compute centroids for each cluster
centroids = cell(num_clusters, 1);
for cluster_index = 1:num_clusters
    cluster_neurons = allclusters{cluster_index, 1}; % Neuron IDs in the cluster
    cluster_signals = F_dff_dec_segmented(:, cluster_neurons); % Extracting signals of these neurons
    centroids{cluster_index} = mean(cluster_signals, 2)'; % Calculating the centroid
end

% Initialize a cell array for unique cluster assignments
unique_clusters = cell(num_clusters, 1);
for i = 1:num_clusters
    unique_clusters{i} = [];
end

% Create a map of all cluster assignments
cluster_assignments = cell(num_neurons, 1);
for cluster_index = 1:num_clusters
    cluster_neurons = allclusters{cluster_index, 1};
    for i = 1:length(cluster_neurons)
        neuron = cluster_neurons(i);
        cluster_assignments{neuron} = [cluster_assignments{neuron}, cluster_index];
    end
end

% Assign neurons to clusters based on the highest correlation with centroid
for neuron = 1:num_neurons
    assigned_clusters = cluster_assignments{neuron};
    if length(assigned_clusters) > 1
        % Compute correlation with each centroid
        neuron_signal = F_dff_dec_segmented(:, neuron)';
        correlations = zeros(1, length(assigned_clusters));
        for idx = 1:length(assigned_clusters)
            cluster_index = assigned_clusters(idx);
            centroid = centroids{cluster_index};
            correlations(idx) = corr(neuron_signal', centroid');
        end
        
        % Find the highest correlation
        [~, max_index] = max(correlations);
        chosen_cluster = assigned_clusters(max_index);
        
        % Assign neuron to the cluster with the highest correlation
        unique_clusters{chosen_cluster} = [unique_clusters{chosen_cluster}, neuron];
    elseif length(assigned_clusters) == 1
        % If the neuron is only in one cluster, keep that assignment
        unique_clusters{assigned_clusters} = [unique_clusters{assigned_clusters}, neuron];
    end
end
% unique_clusters now contains neurons assigned to a single cluster based on the highest correlation criterion


% Initialize allclusters_new as a cell array with the same number of clusters as unique_clusters
allclusters_new = cell(num_clusters, 2);

% Populate allclusters_new with the neuron IDs and calculate centroids
for i = 1:num_clusters
    if ~isempty(unique_clusters{i})
        % Add neuron IDs to the first column
        allclusters_new{i, 1} = unique_clusters{i};

        % Extract the signals for neurons in the current cluster
        cluster_neurons_signals = F_dff_dec_segmented(:, unique_clusters{i});

        % Calculate the centroid for the current cluster
        % Mean across columns (neurons), resulting in a vector for time frames (rows)
        centroid = mean(cluster_neurons_signals, 2)'; % Transpose to make it 1x(number of frames)

        % Add the centroid to the second column
        allclusters_new{i, 2} = centroid;
    else
        % If there are no neurons in a cluster, assign empty arrays
        allclusters_new{i, 1} = [];
        allclusters_new{i, 2} = [];
    end
end
% allclusters_new now resembles the structure of allclusters, with updated unique cluster assignments and centroids


% Print a message saying which neurons, if any, got moved from Cluster X to Cluster Y

% Flatten original cluster assignments for easy comparison
original_cluster_assignments = cell(1, num_neurons);
for cluster_index = 1:num_clusters
    cluster_neurons = allclusters{cluster_index, 1};
    for i = 1:length(cluster_neurons)
        neuron = cluster_neurons(i);
        original_cluster_assignments{neuron} = [original_cluster_assignments{neuron}, cluster_index];
    end
end

% Compare new assignments in unique_clusters with original assignments
neuron_changes = {};
for neuron = 1:num_neurons
    original_assigned_clusters = original_cluster_assignments{neuron};
    new_assigned_cluster = find(cellfun(@(x) ismember(neuron, x), unique_clusters));
    
    % Generate message if there is a change in the cluster assignment
    if ~isequal(original_assigned_clusters, new_assigned_cluster)
        original_clusters_str = sprintf('Cluster %d, ', original_assigned_clusters);
        new_cluster_str = sprintf('Cluster %d', new_assigned_cluster);
        change_message = sprintf('Neuron %d was initially assigned to %s. It is now assigned to %s.', neuron, original_clusters_str(1:end-2), new_cluster_str);
        neuron_changes{end+1} = change_message;
    end
end

% Display the results
if isempty(neuron_changes)
    disp('No neurons were changed from one cluster to another.');
else
    disp('Details of neurons that changed clusters:');
    for i = 1:length(neuron_changes)
        disp(neuron_changes{i});
    end
end

%% Plot of neuron temporal activity patterns (inspired by plots in Dombeck, 2008)
% This section organizes the neurons according to their clusters, 
% visually separates the clusters with black lines, and labels each cluster. 
% This makes it easier to visually assess and interpret the patterns of activity within and between clusters.

% Reorder the data based on the clusters
new_order = [];
for i = 1:length(allclusters_new)
    new_order = [new_order; allclusters_new{i, 1}'];
end

% Load the F_dff_dec_segmented data and reorder it based on the new order
F_reordered = F_dff_dec_segmented(:, new_order);

% Plot the heatmap
figure;
imagesc(F_reordered');
title('Clustered Heatmap of Neurons');
xlabel('Frame #');

% Primary y-axis for neuron index
ylabel('Neuron #', 'Color', 'k');
yyaxis left;

% Set the colors appropriately for neuron index
set(gca, 'ycolor', 'k');
yyaxis right;
axis tight;

cluster_starts = [0; cumsum(cellfun('length', allclusters_new(:, 1)))];
cluster_middles = cluster_starts(1:end-1) + 0.5 * cellfun('length', allclusters_new(:, 1));
yticks(cluster_middles);
cluster_labels = strcat('Cluster ', arrayfun(@num2str, 1:length(allclusters_new), 'UniformOutput', false));
yticklabels(cluster_labels);
set(gca, 'ycolor', 'k');
set(gca, 'YDir','reverse'); % Reverse the y-axis order

% Add vertical black lines indicating cluster boundaries
for start = cluster_starts'
    line([0, size(F_reordered, 1)], [start, start], 'Color', 'k', 'LineWidth', 1);
end

h = colorbar;
ylabel(h, 'ΔF/F');

%% Plot neuron-neuron correlation matrices according to CLUSTERING

% Compute the neuron-to-neuron correlation matrix from F_reordered
correlation_reordered = corr(F_reordered);

% Define cluster boundaries based on allclusters
cluster_starts = [1; cumsum(cellfun('length', allclusters_new(:, 1))) + 1];
cluster_ends = cluster_starts(2:end) - 1;

% Plot the reordered neuron-to-neuron correlation matrix
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
imagesc(correlation_reordered, [-1, 1]);
colormap('jet');
title('Neuron-to-Neuron Correlation Ordered by Clusters');

% Primary y-axis for neuron count
ylabel('Neuron #', 'Color', 'k');
yyaxis left;

% Set the colors appropriately for neuron count
set(gca, 'ycolor', 'k');
yyaxis right;
axis tight;

cluster_middles = (cluster_starts(1:end-1) + cluster_ends) ./ 2;
yticks(cluster_middles);
cluster_labels = strcat('Cluster ', arrayfun(@num2str, 1:length(allclusters_new), 'UniformOutput', false));
yticklabels(cluster_labels);
set(gca, 'ycolor', 'k');
set(gca, 'YDir','reverse'); % Reverse the y-axis order

% Adjust x-axis ticks to start at 0 and increase in counts of 10
% xticks(10:10:length(new_order));
% xticklabels(10:10:length(new_order));
xticks(50:50:length(new_order));
xticklabels(50:50:length(new_order));
xlabel('Neuron #', 'Color', 'k');

% Add vertical and horizontal black lines indicating cluster boundaries
% Skip the first line to avoid overlap with the axes
hold on;
for idx = 2:length(cluster_starts)
    line([cluster_starts(idx), cluster_starts(idx)], [0, length(new_order)], 'Color', 'k', 'LineWidth', 1);
    line([0, length(new_order)], [cluster_starts(idx), cluster_starts(idx)], 'Color', 'k', 'LineWidth', 1);
end
hold off;

h = colorbar;
ylabel(h, 'Correlation Coefficient');

%% Plot neuron-neuron correlation matrices RANDOMLY

% Randomly reorder the neurons
random_order = randperm(size(F_reordered, 2));
F_randomly_ordered = F_reordered(:, random_order);

% Compute the neuron-to-neuron correlation matrix from F_randomly_ordered
correlation_random = corr(F_randomly_ordered);

% Plot the randomly reordered neuron-to-neuron correlation matrix
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
imagesc(correlation_random, [-1, 1]);
colormap('jet');
title('Neuron-to-Neuron Correlation Ordered Randomly');

% Primary y-axis for neuron count
ylabel('Neuron #', 'Color', 'k');
yyaxis left;

% Set the colors appropriately for neuron count
set(gca, 'ycolor', 'k');
yyaxis right;
axis tight;
set(gca, 'YTick', []); % Remove yticks on the right y-axis

% Adjust x-axis ticks to start at 10 and increase in counts of 10
% xticks(10:10:length(random_order));
% xticklabels(10:10:length(random_order));
xticks(50:50:length(random_order));
xticklabels(50:50:length(random_order));
xlabel('Neuron #', 'Color', 'k');

h = colorbar;
ylabel(h, 'Correlation Coefficient');

%% Saving variables and figures

if SAVE_FLAG
    % Save the variables to the specified file
    variable_file_name = 'clustering_data_new.mat';
    full_path = fullfile(current_save_directory, variable_file_name);
    save(full_path, 'allclusters_new', 'neuron_changes', 'remove_columns', 'F_dff_dec', 'F_dff_dec_segmented', 'full_path', 'allclusters', 'centroidcorr', 'dendmem', 'dunnsinitial');
    
    % Save all open figures
    all_figs = findall(0, 'Type', 'figure');
    for idx = 1:length(all_figs)
        fig = all_figs(idx);
        fig_file_name = sprintf('figure_new_%d.png', fig.Number); % You can adjust the format, e.g., .fig, .jpg, etc.
        full_fig_path = fullfile(current_save_directory, fig_file_name);
        saveas(fig, full_fig_path);
    end
end

%% Loop ends and starts over
    % At the end of the loop iteration, remove the added path to keep the workspace clean
    rmpath(genpath(data_path));
    close all;
end
