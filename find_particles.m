%% Initialize settings
clear; clc; close all;

% File path settings
image_dir = 'F:\1202\Aunc-4ms-50ma-75.6fps-+-0-1.4v_files\';  % Image sequence directory
all_init_path = 'F:\1202\all spots.csv';  % First frame particle positions
disp_path = 'F:\1202\one spot.csv';       % Reference particle displacement
plot_dir = 'F:\1202\electrode result\';   % Image save directory

% Parameter configuration
search_radius = 2;        % Search radius (pixels), recommended 2-5
int_radius = 7;           % Intensity extraction radius (pixels), recommended 1-2
min_intensity = 10;       % Valid signal threshold, adjust based on image

%% Create image save directory (if not exists)
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
    fprintf('Created image save directory: %s\n', plot_dir);
end

%% Data loading
fprintf('Loading particle position data...\n');
all_init = readtable(all_init_path);
all_init.Properties.VariableNames{'POSITION_X'} = 'X';
all_init.Properties.VariableNames{'POSITION_Y'} = 'Y';

fprintf('Loading reference particle displacement data...\n');
one_spot = readtable(disp_path);
one_spot.Properties.VariableNames{'POSITION_X'} = 'X';
one_spot.Properties.VariableNames{'POSITION_Y'} = 'Y';

% Calculate displacement relative to first frame
delta_X = one_spot.X - one_spot.X(1);
delta_Y = one_spot.Y - one_spot.Y(1);

% Create displacement data table
disp_data = table(delta_X, delta_Y);
num_frames = height(disp_data);
num_particles = height(all_init);
fprintf('Detected %d particles, %d frames\n', num_particles, num_frames);

% Get image sequence information
image_files = dir(fullfile(image_dir, '*.tif'));
if numel(image_files) ~= num_frames
    error('Frame count mismatch! Displacement data:%d, Images:%d', num_frames, numel(image_files));
end

%% Create result storage structure
particle_data = struct(...
    'ID', num2cell(1:num_particles), ...
    'X', cell(1, num_particles), ...
    'Y', cell(1, num_particles), ...
    'Intensity', cell(1, num_particles), ...
    'Frame', cell(1, num_particles));

% Initialize first frame data
fprintf('Initializing data structure...\n');
for k = 1:num_particles
    particle_data(k).X = all_init.X(k) * ones(num_frames, 1);
    particle_data(k).Y = all_init.Y(k) * ones(num_frames, 1);
    particle_data(k).Intensity = zeros(num_frames, 1);
    particle_data(k).Frame = (1:num_frames)';
end

%% Main processing loop
fprintf('Starting image sequence processing...\n');
progress_bar = waitbar(0, 'Processing image sequence...', 'Name', 'Processing Progress');

for frame = 1:num_frames
    % Read current frame image
    img = imread(fullfile(image_dir, image_files(frame).name));
    
    % Get global displacement
    delta_X = disp_data.delta_X(frame);
    delta_Y = disp_data.delta_Y(frame);
    
    % Process all particles
    for k = 1:num_particles
        % Calculate predicted position
        pred_x = all_init.X(k) + delta_X;
        pred_y = all_init.Y(k) + delta_Y;
        
        % Local window search
        [x, y, intensity] = local_max_search(img, pred_x, pred_y, search_radius, int_radius, min_intensity);
        
        % Store results
        particle_data(k).X(frame) = x;
        particle_data(k).Y(frame) = y;
        particle_data(k).Intensity(frame) = intensity;
    end
    
    % Update progress bar
    progress = frame/num_frames;
    waitbar(progress, progress_bar, sprintf('Processing progress: %d/%d (%.1f%%)', frame, num_frames, progress*100));
end

% Close progress bar
close(progress_bar);
fprintf('Image processing completed!\n');

%% Result output
fprintf('Saving tracking results...\n');

% Preallocate result table
total_rows = num_particles * num_frames;
results = table(...
    zeros(total_rows, 1), ...  % ParticleID
    zeros(total_rows, 1), ...  % Frame
    NaN(total_rows, 1), ...    % X
    NaN(total_rows, 1), ...    % Y
    zeros(total_rows, 1), ...  % Mean_Intensity
    'VariableNames', {'ParticleID', 'Frame', 'X', 'Y', 'Mean_Intensity'});

% Fill preallocated table
for k = 1:num_particles
    start_idx = (k-1)*num_frames + 1;
    end_idx = k*num_frames;
    
    results.ParticleID(start_idx:end_idx) = k;
    results.Frame(start_idx:end_idx) = particle_data(k).Frame;
    results.X(start_idx:end_idx) = particle_data(k).X;
    results.Y(start_idx:end_idx) = particle_data(k).Y;
    results.Mean_Intensity(start_idx:end_idx) = round(particle_data(k).Intensity, 2);
end

% Save results
results_file = 'single_particle_tracking_results.csv';
writetable(results, results_file);
save('particle_tracking_data.mat', 'particle_data');
fprintf('Results saved to: %s\n', results_file);

%% Plot and save intensity curves
fprintf('Starting intensity curve plotting...\n');
plot_progress_bar = waitbar(0, 'Plotting intensity curves...', 'Name', 'Plotting Progress');

% Create unified figure template
fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600], 'Color', 'white');
curve_color = [0, 0.45, 0.74]; % MATLAB default blue

for k = 1:num_particles
    % Prepare data
    frames = particle_data(k).Frame;
    intensities = particle_data(k).Intensity;
    
    % Create new figure
    clf(fig);
    plot(frames, intensities, 'Color', curve_color, 'LineWidth', 1.8);
    
    % Set figure properties
    title(sprintf('Particle %d Average Intensity Curve', k), 'FontSize', 16, 'FontName', 'Arial');
    xlabel('Frame Number', 'FontSize', 14, 'FontName', 'Arial');
    ylabel('Average Intensity (a.u.)', 'FontSize', 14, 'FontName', 'Arial');
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'GridLineStyle', '--', ...
        'GridAlpha', 0.3, 'FontName', 'Arial');
    xlim([1, num_frames]);
    
    % Auto-adjust Y-axis range
    valid_intensities = intensities(~isnan(intensities));
    if ~isempty(valid_intensities)
        y_min = min(valid_intensities);
        y_max = max(valid_intensities);
        if y_min == y_max
            y_upper = max(1, y_max+1);
        else
            margin = 0.1 * (y_max - y_min);
            y_upper = y_max + margin;
        end
        ylim([0, y_upper]);
    end
    box on;
    
    % Save image
    plot_filename = fullfile(plot_dir, sprintf('Particle_%04d.png', k));
    exportgraphics(fig, plot_filename, 'Resolution', 300);
    
    % Update plotting progress
    waitbar(k/num_particles, plot_progress_bar, sprintf('Plotting progress: %d/%d (%.1f%%)', k, num_particles, k/num_particles*100));
end

% Close progress bar and figure
close(plot_progress_bar);
close(fig);
fprintf('Intensity curve plotting completed!\n');
fprintf('Curve plots saved to: %s\n', plot_dir);

%% Local maximum search function
function [x, y, intensity] = local_max_search(img, pred_x, pred_y, radius, int_radius, min_int)
    % Get search window
    x_min = max(1, floor(pred_x) - radius);
    x_max = min(size(img, 2), ceil(pred_x) + radius);
    y_min = max(1, floor(pred_y) - radius);
    y_max = min(size(img, 1), ceil(pred_y) + radius);
    
    % Check window validity
    if x_min > x_max || y_min > y_max
        x = NaN;
        y = NaN;
        intensity = 0;
        return;
    end
    
    % Extract search region
    search_roi = img(y_min:y_max, x_min:x_max);
    
    % Find local maximum
    [max_val, max_idx] = max(search_roi(:));
    
    % Check signal validity
    if max_val < min_int
        x = NaN;
        y = NaN;
        intensity = 0;
        return;
    end
    
    % Calculate local maximum position
    [rel_y, rel_x] = ind2sub(size(search_roi), max_idx);
    x = x_min + rel_x - 1;
    y = y_min + rel_y - 1;
    
    % Extract average intensity
    intensity = extract_intensity(img, x, y, int_radius);
end

%% Intensity extraction function (circular region average intensity)
function avg_intensity = extract_intensity(img, x, y, radius)
    % If radius < 1, return center point intensity
    if radius < 1
        avg_intensity = double(img(round(y), round(x)));
        return;
    end
    
    % Calculate center coordinates
    center_x = round(x);
    center_y = round(y);
    
    % Determine bounding box
    x_min = max(1, center_x - radius);
    x_max = min(size(img, 2), center_x + radius);
    y_min = max(1, center_y - radius);
    y_max = min(size(img, 1), center_y + radius);
    
    % Initialize total intensity and pixel count
    total_intensity = 0;
    pixel_count = 0;
    
    % Traverse pixels within bounding box
    for xi = x_min:x_max
        for yi = y_min:y_max
            dist = sqrt((xi - center_x)^2 + (yi - center_y)^2);
            if dist <= radius
                total_intensity = total_intensity + double(img(yi, xi));
                pixel_count = pixel_count + 1;
            end
        end
    end
    
    % Calculate average intensity
    if pixel_count > 0
        avg_intensity = total_intensity / pixel_count;
    else
        avg_intensity = 0;
    end
end