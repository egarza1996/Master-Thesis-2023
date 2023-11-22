clear;
clc;
close all;

% Code for conditions 2/3 (50msec every 20sec, and 1msec every 20sec)
% Note: Because this code saves a lot of videos, sometimes it is best to run
% sections manually just to make sure that they get saved properly.

%% Add paths and confirm the run that you want to analyze - only things you need to modify

% Define the root path and the organoid folders
root_path = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/LOOP/exp1-2_11-25_07_23/'; % contains the results file we want to analyze
organoid_folders = arrayfun(@(x) sprintf('organoid%d', x), 1:9, 'UniformOutput', false);

addpath(genpath('/home/silviu/matlab_code')); % contains shadedErrorBar.m function
addpath(genpath('/mnt/data/erika-organoid-data/CODE/CURRENT-FOLDER')); % contains plotNeurons.m function
addpath(genpath('/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Average/loop'));

% Set SAVE_FLAG to 1 if you want to save variables, otherwise set to 0
SAVE_FLAG = 1;

% Define saving path
base_save_directory = '/home/silviu/erika-organoid-data/CODE/CURRENT-FOLDER/Average/loop';

% Loop through each organoid folder
for org_idx = 1 % 1:length(organoid_folders)
    if org_idx == 3 || org_idx == 7 || org_idx == 8 % these organoid recordings are very bad
        continue;
    end
    tic; % check for elapsed time of running the code
    current_folder = organoid_folders{org_idx};
    
    % Update the paths for the current organoid folder and "run_2" data file
    data_path = fullfile(root_path, current_folder, 'run_2');
    addpath(genpath(data_path)); % contains 'analysis_results.hdf5' file

    % Create the subfolder for the current organoid and run
    current_save_directory = fullfile(base_save_directory, current_folder, 'run_2');
    
    % Check if the folder exists, if not, create it
    if ~exist(current_save_directory, 'dir')
        mkdir(current_save_directory);
    end

%% Using the following signal: F_dff

F_dff = h5read('analysis_results.hdf5', '/estimates/F_dff'); % normalized "delta F over F" values
h5disp('analysis_results.hdf5', '/estimates/F_dff');

%% Create .avi file with all neuron traces 
% You can then open it in Fiji to visualize each frame

% Set folder and filename where you want to save the video
name = 'AllNeuronTraces.avi';
fullpath = fullfile(current_save_directory, name);
vidfile = VideoWriter(fullpath, 'Motion JPEG AVI');  %'Uncompressed AVI'

open(vidfile);

% figure('Position', [10 10 1200 600])
figure('Position', get(0, 'ScreenSize')) % makes the figure enter full screen, which is better for saving the video

fontSize = 14;  % Increase this value if you want bigger text
for i=1:size(F_dff,2)
    plot(F_dff(:,i)*100);  % Multiply y values by 100
    title(strcat('Neuron #', num2str(i)), 'FontSize', fontSize + 2);
    xlabel('Frame Number', 'FontSize', fontSize); 
    ylabel('% ΔF/F', 'FontSize', fontSize);
    
    % Set x-axis limit
    xlim([1, size(F_dff,1)]);
    
    % Set font size for tick marks
    ax = gca;
    ax.FontSize = fontSize;

    set(gcf,'color', 'w');
    drawnow;
    y = getframe(gcf);
    writeVideo(vidfile, y)
    clf
end

close(vidfile);
close;


%% Compute the average neuron values across stimulation cycles

AllTiffFiles = dir('*.tif');
Fall = cell(1, length(AllTiffFiles)); % Preallocate the cell array
for i_tf = 1:length(AllTiffFiles)
    Fall{i_tf} = loadtiff(AllTiffFiles(i_tf).name);   
end

    % f = cat(3,Fall{1},Fall{2},Fall{3},Fall{4},Fall{5}); % In the case that there are 5 TIF files
f = cat(3, Fall{:}); % Automatically concatenate all the TIF files
f = double(f); % Convert to double

neurons_stim_cycles = {}; % initialize variable -- represents the averaged values of a neuron across stimulus cycles.

stim_start = 600:400:11600; % cond2 and cond3

% frame rate is 20 frames/sec
% stim starts at 30sec for all conditions, so we start at frame 600 (i.e. 30*20)
% middle number indicates how often there was stimulation:
    % For cond1 every 10sec, so frame is 200 (i.e. 10*20)
    % For cond2/3 every 20sec, so frame is 400 (i.e. 20*20)
% last number indicates until which frame it should go:
        % For cond1, duration was 5min30sec (i.e. 330sec*20 = 6600 frames)
        % For cond2/3, duration was 10min30sec (i.e. 630sec*20 = 12600 frames)
% IMPORTANT: we want to analyze the window corresponding to 1 second before stimulation, and 10sec or 20sec after stimulation
    % For con1, use 5800 (=6600-30*20-10*20)
    % For cond2/3, use: 11600 (=12600-30*20-20*20)

av_neuron = zeros(420,28); % cond2/3

% For cond1: we want to analyze the window corresponding to 1 second
% (i.e. 20 frames) before stimulation, and 10sec (i.e. 200 frames) after
% stimulation, which gives a total "window" length of 220 frames.
% Therefore, 5800/220=27.36. This means the signal will be divided into 27 sections.

% For cond2/3: we want to analyze the window corresponding to 1 second
% (i.e. 20 frames) before stimulation, and 20sec (i.e. 400 frames) after
% stimulation, which gives a total "window" length of 420 frames.
% Therefore, 11600/420=27.61. This means the signal will be divided into 28 sections.

for i = 1:size(F_dff,2)
    f = F_dff(:,i); % analyze neuron by neuron
    k = 0;

  for m = stim_start(1:end)
      k = k+1;    
      
      display(strcat('processing......... Neuron # ',num2str(i), '... Stimulation Cycle # ' ,num2str(k)))
     
      av_neuron(:,k) = f(m-20:m+399); % cond2/3: 1sec (20frames) before stim, and 20sec after (400frames)  

  end

    neurons_stim_cycles{i} = av_neuron;

end


%% Plot average with stdev of NEURON #50 's average  
figure;
% 19.95 for cond2/3 (corresponds to stimulus cycle time length)
shadedErrorBar(-1:0.05:19.95, neurons_stim_cycles{50}'*100, {@mean,@std}); % TIME and also multiply by 100 because of % ΔF/F
hold on;

% 20sec for cond2/3
xlim([-1 20]); 
% ylim([ymin ymax]);
vline(0, 'r'); % indicates optogenetic stimulation at t=0sec

fontSize = 14;  % Increase this value if you want bigger text
title("Neuron #50 - Average with Standard Deviation", 'FontSize', fontSize + 2);
xlabel('Time (sec)', 'FontSize', fontSize);
ylabel('% ΔF/F', 'FontSize', fontSize);

% Set font size for tick marks
ax = gca;
ax.FontSize = fontSize;

% Another way of viewing this with frames instead of time on the x axis...
    % 420 for cond2/3
% figure;
% shadedErrorBar(1:420, neurons_stim_cycles{50}'*100, {@mean,@std}); % FRAMES
% hold on;
% xlim([0 420]); % ylim([ymin ymax]);
% vline(20, 'r'); % indicates optogenetic stimulation at t=1sec which is equal to frame #20
% title("Neuron #50 - Average with Standard Deviation");
% xlabel('Frame Number');
% ylabel('% ΔF/F');


%% Plot average with stdev of ALL NEURONS averages
figure;
% 19.95 for cond2/3 (corresponds to stimulus cycle time length)
shadedErrorBar(-1:0.05:19.95, av_neuron'*100, {@mean,@std}); % TIME and also multiply by 100 because of % ΔF/F
hold on;
    
% 20sec for cond2/3
xlim([-1 20]); 
% ylim([ymin ymax]);
vline(0, 'r'); % indicates optogenetic stimulation at t=0sec

fontSize = 14;  % Increase this value if you want bigger text
title("All Neurons - Average with Standard Deviation", 'FontSize', fontSize + 2);
xlabel('Time (sec)', 'FontSize', fontSize);
ylabel('% ΔF/F', 'FontSize', fontSize);

% Set font size for tick marks
ax = gca;
ax.FontSize = fontSize;

% Another way of viewing this...
    % 420 for cond2/3
% figure;
% shadedErrorBar(1:420, av_neuron', {@mean,@std}); % FRAMES
% hold on;
% xlim([0 420]); % ylim([ymin ymax]);
% vline(20, 'r'); % indicates optogenetic stimulation at t=1sec which is equal to frame #20
% title("Neuron #50 - Average with Standard Deviation");
% xlabel('Frame Number');
% ylabel('% ΔF/F');


%% Plot 10 neurons 
    % See function (made by me): plotNeurons(av_neuron, startNeuron, endNeuron, figureTitle, xLimit)
        
% xLimit is 420 for cond2/3    
% yLimits sometimes also need to be changed but do that directly in the function if needed

% Plot neurons 1 through 5
figure;
plotNeurons(av_neuron*100, 1, 5, "Average Neuron Values Across 5 Stimulation Cycles (SC)", 420);

% Plot neurons 10 through 20
figure;
plotNeurons(av_neuron*100, 10, 20, "Average Neuron Values Across 10 Stimulation Cycles (SC)", 420);


%% Averaged neuronal activity movie (basically averaging stimulations across the whole movie recording, thereby summarizing it into a couple of secs)

% Identify all the tif files of that organoid
file_pattern = '*.tif';
AllTiffFiles = dir(fullfile(data_path, file_pattern));

Fall = cell(1, length(AllTiffFiles)); % Preallocate the cell array
for i_tf = 1:length(AllTiffFiles)
    Fall{i_tf} = loadtiff(AllTiffFiles(i_tf).name);   
end

f = cat(3, Fall{:}); % Automatically concatenate all the TIF files
f = double(f); % Convert to double


av_movie = zeros(size(f,1),size(f,2), 420); % cond2/3


k = 0;
  for m = stim_start(1:end)
      k = k+1;
      
      display(strcat('processing...........',num2str(k)))


      temp = f(:,:,m-20:m+399); % cond2/3

  
    av_movie = av_movie + temp;
   

  end
    
av_movie = av_movie./k;
% FYI: The av_movie variable is a 3D matrix with dimensions XxYxZ, where XxY are the pixel resoultion and Z is the number of frames in the movie.

% Use implay function to see the video or use code below to save video
% implay(av_movie); % remember to change parameter in Tools -- Colormap -- Specify range -- 0 to 2200 (or whatever # in that range works)


%% Save "Averaged neuronal activity movie" video as an .avi file

hFig = figure; % Create a new figure
colormap_range = [0, 2200]; % Specify the range for the colormap
hImage = imagesc(av_movie(:, :, 1)); % Display the first frame
axis off; % Turn off the axis labels
colormap('gray'); % Set the colormap (you can choose the one you prefer)
clim(colormap_range); % Set the colormap range
title('Averaged neuronal activity'); % Set the title

% Set folder and filename where you want to save the video
filename = 'AveragedNeuronalActivity';
fullpath = fullfile(current_save_directory, filename);
video = VideoWriter(fullpath);
open(video);

for i = 1:size(av_movie, 3) % Loop over each frame
    set(hImage, 'CData', av_movie(:, :, i)); % Update the frame
    drawnow; % Force MATLAB to render the frame
    frame = getframe(hFig); % Capture the figure window as a frame
    writeVideo(video, frame); % Write the frame to the video
end

% Close the video writer object
close(video);
close(hFig);


%% Create .avi file of the video above with its corresponding signal
% (sometimes you need to run this section several times until it works, not sure why)

name = 'AveragedNeuronalActivityWithSignal';
fullpath = fullfile(current_save_directory, name);
vidfile = VideoWriter(fullpath, 'Motion JPEG AVI');  %'Uncompressed AVI'

% Stimulation Cycle duration
time = -1:0.05:19.95; % cond2/3

roi = squeeze(mean(mean(av_movie,1),2));

open(vidfile);
figure('Position', [10 10 1200 600])
% figure('Position', get(0, 'ScreenSize')) % makes the figure enter full screen, which is better for saving the video

for i = 1:size(av_movie,3)
    
subplot(4,20,1:60);
     imagesc(av_movie(:,:,i), [0 2200]);  axis image; axis off;colormap(gray); %% 3D data goes here
                        % Change 2200 for better visualization, if needed
    title(['Frame #', num2str(i)]);
   
 subplot(4,20,62:79);  

 plot(time, roi, 'LineWidth', 1);  vline(time(i),'k');
 
 vline(0,'r');
 % line([0, 0], [520, 560], 'Color', 'r') % sometimes this line works better than the line above

 % ylim([520 560]); % sometimes you need to manually
 xlim([-1 20]);


 xlabel('Time (sec)'); 
 ylabel('(arbitrary units)'); % arbitrary units

% set(gcf,'color', 'w');
drawnow;
y = getframe(gcf);
writeVideo(vidfile, y)
clf
end

close(vidfile);


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
    variable_file_name = 'average_data.mat';
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
