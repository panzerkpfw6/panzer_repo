% Clear workspace and close figures
clear;
close all;

% Parameters
dx = 25.0; % Inline spacing (m)
dy = 25.0; % Crossline spacing (m)
filename = '../../data/shots.txt'; % Path to shots.txt (adjust if needed)
picture_dir = 'output_plots'; % Directory to save the plot (adjust as needed)
picture_dir='/media/plotnips/sdd1/Dropbox/Apps/Overleaf/2025_stencil_rtm/figures/rtm_salt/'

% Create output directory if it doesn't exist
if ~exist(picture_dir, 'dir')
    mkdir(picture_dir);
end

% Read the shots.txt file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open shots.txt. Check the file path.');
end

% Initialize arrays to store coordinates
num_shots = 101; % From shot ID 0 to 100
isx = zeros(num_shots, 1);
isy = zeros(num_shots, 1);
shot_ids = zeros(num_shots, 1);

% Parse the file
for i = 1:num_shots
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    % Extract values using regular expressions
    tokens = regexp(line, 'ir=\s*(\d+),\s*shot id=\s*(\d+),\s*isx=(\d+),\s*isy=(\d+),\s*isz=(\d+),\s*isz_rcv=(\d+)', 'tokens');
    if ~isempty(tokens)
        shot_ids(i) = str2double(tokens{1}{2});
        isx(i) = str2double(tokens{1}{3});
        isy(i) = str2double(tokens{1}{4});
    end
end
fclose(fid);

% Convert indices to physical coordinates (in meters)
x = isx * dx;
y = isy * dy;

% Convert to kilometers for plotting
x_km = x / 1000;
y_km = y / 1000;

% Create the figure with a white background
f1 = figure('units', 'pixels', 'position', [2003 303 1200 600], 'Color', [1 1 1]);

% Plot the shot locations as a scatter plot
scatter(x_km, y_km, 50, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
ax1 = gca();

% Set title and labels with LaTeX interpreter
title('\textbf{Shot Locations}', 'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 19);
xlabel('$\textbf{X, km}$', 'Interpreter', 'latex', 'FontWeight', 'bold');
ylabel('$\textbf{Y, km}$', 'Interpreter', 'latex', 'FontWeight', 'bold');

% Set font size and grid
set(ax1, 'FontSize', 19);
grid on;

% Set axes background to white
set(ax1, 'Color', [1 1 1]);

% Make the axis tight
axis tight;

% Minimize white space around the axes
set(ax1, 'LooseInset', get(ax1, 'TightInset') + [0.02 0.02 0.02 0.02]);

% Adjust figure for saving
set(f1, 'PaperPositionMode', 'auto');
set(f1, 'InvertHardcopy', 'off');

% Save the figure
% print(f1, fullfile(picture_dir, 'shot_locations'), '-dpng', '-r400');
exportgraphics(f1, fullfile(picture_dir, 'acquisition.png'), 'Resolution', 400); %, 'BackgroundColor', 'white'
