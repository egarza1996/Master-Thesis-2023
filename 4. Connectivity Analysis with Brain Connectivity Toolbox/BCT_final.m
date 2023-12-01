clear;
clc;
close all;

% This code merges "Pre_BCT_1.m" and "BCT_2.m" and creates a loop across
% all organoid folders and saves important figures and variables.

% Analyzes the WHOLE signal. Check "BCT_final_segmented" for analysis of
% only the post-stimulation-windows for conditions 1, 2, and 3.

%% Part 1: Correlation matrices
% 0. Add paths and define root/saving paths
% 1. Loaded CaImAn's output file.
% 2. Computed the pairwise Pearson correlation coefficients.
% 3. Visualized the full correlation matrix.
% 4. Applied 70% proportional thresholding to the correlation matrix.
% 5. Visualized the thresholded correlation matrix.
% continue with Part 2

%% Add paths and define F -- only things you need to modify

% Define the root path and the organoid folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/LOOP/exp1-2_11-25_07_23/';
organoid_folders = arrayfun(@(x) sprintf('organoid%d', x), 1:9, 'UniformOutput', false);
 
addpath('/mnt/data/erika-organoid-data/CODE/CURRENT-FOLDER/PreBCT-and-BCT/loop');
addpath('/home/silviu/BCT/2019_03_03_BCT'); % contains BCT toolbox
addpath('/home/silviu/schemaball'); % contains schemaball

% Set SAVE_FLAG to 1 if you want to save variables, otherwise set to 0
SAVE_FLAG = 1;

% Define saving path
base_save_directory = '/mnt/data/erika-organoid-data/CODE/CURRENT-FOLDER/PreBCT-and-BCT/loop';

% Define the different runs you are interested in
% runs = {'run', 'run_1', 'run_2', 'run_3', 'run_4', 'run_5'};
% runs = {'run', 'run_4', 'run_5'};
runs = {'run_5'};
% runs = {'run_1', 'run_2', 'run_3'};

% Loop through each organoid folder
for org_idx = 1:length(organoid_folders)
    if org_idx == 3 || org_idx == 7 || org_idx == 8 % these organoid recordings are very bad
        continue;
    end
    current_folder = organoid_folders{org_idx};
    
    % Loop through each run
    for run_idx = 1:length(runs)
        % Loop through the desired analysis runs
        desired_run = runs{run_idx};

        % Update the paths for the current organoid folder and run data file
        data_path = fullfile(root_path, current_folder, desired_run);
        addpath(genpath(data_path)); % contains 'analysis_results.hdf5' file
        
        % Create the subfolder for the current organoid and run
        current_save_directory = fullfile(base_save_directory, current_folder, desired_run);
        
        % Check if the folder exists, if not, create it
        if ~exist(current_save_directory, 'dir')
            mkdir(current_save_directory);
        end


        %% Step 1: Load CaImAn's output file "analysis_results.hdf5"
        % Load C (raw fluorscence traces)
        % C = h5read('analysis_results.hdf5', '/estimates/C');
        
        % OR Load F_dff_dec % deconvolved fluorescence signal, representing an estimate of the underlying neural activity
        F_dff_dec = h5read('analysis_results.hdf5', '/estimates/F_dff_dec');
        h5disp('analysis_results.hdf5', '/estimates/F_dff_dec');
        
        % Get the number of rows in F_dff_dec
        num_rows = size(F_dff_dec, 1);
        
        % Slicing logic based on run ID
        if strcmp(runs, 'run')
            % Slice for "run" if the number of rows is sufficient
            F_dff_dec = F_dff_dec(1:min(num_rows, 6600), :); % 5min 30sec
        elseif strcmp(runs, 'run_1')
            % Slice for "run_1" if the number of rows is sufficient
            F_dff_dec = F_dff_dec(600:min(num_rows, 6600), :); % 5min (discard first 30sec of baseline activity)
        elseif strcmp(runs, 'run_2')
            % Slice for "run_2" if the number of rows is sufficient
            F_dff_dec = F_dff_dec(600:min(num_rows, 12600), :); % 10min (discard first 30sec of baseline activity)
        elseif strcmp(runs, 'run_3')
            % Slice for "run_3" if the number of rows is sufficient
            F_dff_dec = F_dff_dec(600:min(num_rows, 12600), :); % 10min (discard first 30sec of baseline activity)
        elseif strcmp(runs, 'run_4')
            % Slice for "run_4" if the number of rows is sufficient
            F_dff_dec = F_dff_dec(1200:min(num_rows, 12600), :); % 9.5min (discard first 1min of baseline activity)            
        else
            % Slice for "run_5" if the number of rows is sufficient
            F_dff_dec = F_dff_dec(1200:min(num_rows, 12600), :); % 9.5min (discard first 1min of baseline activity)
        end
        
        
        %% Step 2: Compute the pairwise Pearson correlation coefficients
            % C corresponds to neurons (columns) and C' to time points (rows)
        correlation_matrix = corr(F_dff_dec);
        
        %% Remove possible NaN values in the matrix
        % Initialize an array to hold the indices of rows/columns to be removed
        indices_to_remove = [];
        
        % Loop through each row
        for i = 1:size(correlation_matrix, 1)
            if all(isnan(correlation_matrix(i, :)))
                indices_to_remove(end + 1) = i;
            end
        end
        
        % Remove the identified rows and corresponding columns
        correlation_matrix(indices_to_remove, :) = [];
        correlation_matrix(:, indices_to_remove) = [];
        
        %% Step 3: Visualize the Full Correlation Matrix
        figure;
        imagesc(correlation_matrix);
        colorbar;
        clim([-1 1]);
        title('Full Correlation Matrix');
        xlabel('Neurons');
        ylabel('Neurons');
        
        %% Step 4: Apply 70% Proportional Thresholding using the "threshold_proportional.m" function from BCT
        % i.e. preserve the top 70% strongest connections in the matrix and remove the remaining 30%
        
        thresholded_matrix = threshold_proportional(correlation_matrix, 0.70);
        
        %% Step 5: Visualize the Thresholded Correlation Matrix
        figure;
        imagesc(thresholded_matrix);
        colorbar;
        clim([-1 1]);
        title('Proportionally Thresholded Correlation Matrix (by 70%)');
        xlabel('Neurons');
        ylabel('Neurons');
        
        
        %% Part 2: After getting the correlation matrices, this code runs a Connectivity Analysis with BCT functions:
        % 1) Degree Distribution: Number of connections or edges the node has to other nodes (neurons).
        % 2) Clustering Coefficient: Measures the degree to which nodes in a graph tend to cluster together. It gives an idea of the cliquishness of a typical neighborhood in the network.
        % 3) Path Length: Average shortest path length between all pairs of nodes.
        % 4) Modularity: Detects communities in networks, where nodes in the same community are more densely connected to each other than to nodes in other communities. Modularity is a statistic that quantifies the degree to which the network may be subdivided into clearly delineated groups or modules.
        % 5) Centrality Measures: These can identify the most influential nodes in a network. (Betweenness Centrality: A measure of centrality in a graph based on shortest paths.)
        % 6) Small-Worldness: A network is considered small-world if it has a clustering coefficient significantly higher than random networks, and a similar characteristic path length to random networks. To determine small-worldness, you'd typically compare the clustering coefficient and path length of your network to those of random networks.
        % 7) Extra: Schemaball
        
        % Choose matrix (i.e. full correlation matrix or thresholded correlation matrix)
        % matrix = thresholded_matrix;
        
        % Abosulute value (added on 10.11.23)
        matrix = abs(thresholded_matrix);
        
            % QUESTION: since there are a lot of negative values in the correlation
            % matrix (and just a few in the thresholded) should I just do the absolute
            % value from the start? i.e.: matrix = abs(correlation_matrix);
            % Yes, you can take the abs from the start, but make sure to mention
            % this in the discussion section of thesis.
        
        %% 1) Degree Distribution
        degree_dist = degrees_und(matrix);
        figure;
        histogram(degree_dist);
        title('Degree Distribution');
        xlabel('Degree');
        ylabel('Number of Nodes');
        
        %% 2) Clustering Coefficient
        clustering_coeff = clustering_coef_wu_sign(matrix);
        figure;
        bar(clustering_coeff);
        title('Clustering Coefficient');
        xlabel('Node');
        ylabel('Coefficient');
        
        %% 3) Path Length
        D = distance_wei(1./abs(matrix)); % Convert weights to distances
            % Note:  Since we are dealing with a signed correlation matrix (with both positive and negative values), 
            % one common approach is to take the absolute value of the matrix before computing distances.
        path_length = mean(mean(D));
        disp(['Average Path Length: ', num2str(path_length)]);
        
        %% 4) Modularity
        [community_structure, ~] = community_louvain(matrix, [], [], 'negative_sym');
            % 'negative_sym': This method treats the network as undirected (the connections between nodes do not have a direction) with signed weights. 
        figure;
        imagesc(community_structure);
        title('Community Structure');
        xlabel('Node');
        ylabel('Node');
        colorbar;
        
        %% 5) Betweenness Centrality
        betweenness_centrality = betweenness_wei(1./abs(matrix)); % Convert weights to distances
            % While running this section on the spont. data (i.e. run 0) had to add abs(matrix) in order for this fn to run
        figure;
        bar(betweenness_centrality);
        title('Betweenness Centrality');
        xlabel('Node');
        ylabel('Centrality');
        
        %% 6) Small-Worldness 
        % Takes a very long time to run on full datasets.
        num_random_networks = 20;
        clustering_coeff_random = zeros(1, num_random_networks);
        path_length_random = zeros(1, num_random_networks);
        for i = 1:num_random_networks
            random_network = randmio_und(abs(matrix), 5);
            clustering_coeff_random(i) = mean(clustering_coef_wu(random_network));
            D_random = distance_wei(1./random_network);
            path_length_random(i) = mean(mean(D_random));
        end
        sigma = (mean(clustering_coeff) / mean(clustering_coeff_random)) / (path_length / mean(path_length_random));
        disp(['Small-Worldness: ', num2str(sigma)]);
        
        %% 7) Extra: Schemaball
        % Ideally, in optogenetic stimulation we want to see changes 
        % in the network dynamics from "small world" to "fully connected"
        
        % schemaball(correlation_matrix);
        
        schemaball(thresholded_matrix);
        
        
        %% Saving variables and figures
        
        if SAVE_FLAG
            % Save the variables to the specified file
            variable_file_name = 'BCT.mat';
            full_path = fullfile(current_save_directory, variable_file_name);
            save(full_path, 'F_dff_dec', 'indices_to_remove', 'correlation_matrix', 'thresholded_matrix', 'degree_dist', 'clustering_coeff', 'path_length', 'community_structure', 'betweenness_centrality', 'sigma');
            
            % Save all open figures
            all_figs = findall(0, 'Type', 'figure');
            for idx = 1:length(all_figs)
                fig = all_figs(idx);
                fig_file_name = sprintf('figure_%d.png', fig.Number); % You can adjust the format, e.g., .fig, .jpg, etc.
                full_fig_path = fullfile(current_save_directory, fig_file_name);
                saveas(fig, full_fig_path);
            end
        end
        
        
        %% Loop ends and starts over
            % At the end of the loop iteration, remove the added path to keep the workspace clean
            rmpath(genpath(data_path));
            close all;
    end
end
