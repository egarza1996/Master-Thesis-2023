clear;
clc;
close all;

%% Creating a .CSV file from the previously calculated BCT parameters for all organoids and all runs

% In this script:
%     You set the parameter_name variable to the name of the parameter you want to extract.
%     The script checks if the parameter is one of those that require the mean to be calculated.
%     It extracts either the mean value or the single value, as appropriate.
%     It saves the results to a CSV file named after the parameter.


% Define the root path and the organoid folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/PreBCT-and-BCT/loop';
organoid_folders = {'organoid1', 'organoid2', 'organoid4', 'organoid5', 'organoid6', 'organoid9'};
runs = {'run', 'run_1', 'run_1_segmented', 'run_2', 'run_2_segmented', 'run_3', 'run_3_segmented', 'run_4', 'run_5'};

% Parameter to extract from the BCT.mat file
parameter_name = 'sigma';  % Change this to the parameter you want to extract
    % "BCT.mat" contains the following parameters:
    % 'F_dff_dec', 'indices_to_remove', 'correlation_matrix', 'thresholded_matrix', ...
    % 'degree_dist', 'clustering_coeff', 'path_length', 'community_structure', 'betweenness_centrality', 'sigma');

% Parameters that need mean calculation
mean_parameters = {'degree_dist', 'clustering_coeff', 'community_structure', 'betweenness_centrality'};

% Initialize a table to store the parameter values
num_organoids = length(organoid_folders);
num_runs = length(runs);
extract = array2table(zeros(num_organoids, num_runs), 'VariableNames', runs, 'RowNames', organoid_folders);

% Loop through each organoid folder
for i = 1:num_organoids
    organoid = organoid_folders{i};
    
    % Loop through each run
    for j = 1:num_runs
        run = runs{j};
        
        % Construct the full path to the BCT.mat file
        mat_file_path = fullfile(root_path, organoid, run, 'BCT.mat');
        
        % Check if the BCT.mat file exists
        if exist(mat_file_path, 'file')
            data = load(mat_file_path);
            
            % Check if the parameter is a field in the loaded data
            if isfield(data, parameter_name)
                % Determine if the parameter is a single value or needs mean calculation
                if any(strcmp(mean_parameters, parameter_name))
                    % Calculate the mean value of the parameter
                    value = mean(data.(parameter_name)(:));
                elseif any(strcmp("F_dff_dec", parameter_name))
                    % Extract the total number of neurons from "F_dff_dec"
                    value = data.(parameter_name);
                    value = size(value);
                    value = value(1,2);
                else
                    % Extract the single value parameter directly
                    value = data.(parameter_name);
                end
                
                % Store the value of the parameter in the table
                extract{organoid, run} = value;
            else
                warning('%s not found in %s', parameter_name, mat_file_path);
            end
        else
            warning('%s does not exist.', mat_file_path);
        end
    end
end

% Save the table to a CSV file with the specified parameter name
csv_file_path = fullfile(root_path, [parameter_name '_values.csv']);
writetable(extract, csv_file_path, 'WriteRowNames', true);


%% Optional: calculating the proportional degree distribution
% Number of connections the neuron has to other neurons
% divided by the the number of total connections (i.e. neurons^2) 
% minus the number of repeated neurons.
% (Like a correlation matrix, where the diagonal values are not considered.)

% % Load data from CSV files
% degreeDist = readmatrix('degree_dist_values.csv');
% neurons = readmatrix('F_dff_dec_values.csv'); % contains the total number of neurons
% 
% % Check if the dimensions of the two matrices are the same
% if all(size(degreeDist) == size(neurons))
%     % Perform element-wise division
%     result = degreeDist ./ ((neurons).^2 - neurons) ;
% 
%     % Handle any division-by-zero by replacing Inf and NaN with a defined value
%     % Here, we choose to replace them with zero
%     result(isinf(result) | isnan(result)) = 0;
% 
%     % Append the file name to the directory path
%     outputFilePath = fullfile(root_path, 'proportional_degree_dist_values.csv');
% 
%     % Save the result to the new CSV file at the specified path
%     writematrix(result, outputFilePath);
% else
%     error('The dimensions of the two CSV files do not match.');
% end
