function processSPRiDataFromImages()
    image_folder = uigetdir('', 'Select folder containing image sequence');
    if isequal(image_folder, 0)
        disp('User cancelled folder selection');
        return;
    end
    
    fprintf('Processing folder: %s\n', image_folder);
    
    output_folder = uigetdir('', 'Select folder to save results');
    if isequal(output_folder, 0)
        output_folder = fullfile(image_folder, 'SPRi_Analysis_Results');
    end
    
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    image_files = dir(fullfile(image_folder, '*.tif*'));
    if isempty(image_files)
        image_files = dir(fullfile(image_folder, '*.png'));
    end
    if isempty(image_files)
        image_files = dir(fullfile(image_folder, '*.jpg'));
    end
    
    if isempty(image_files)
        errordlg('No image files found (supported formats: tif, tiff, png, jpg)', 'File Error');
        return;
    end
    
    [~, idx] = sort({image_files.name});
    image_files = image_files(idx);
    
    num_frames = length(image_files);
    fprintf('Found %d image files\n', num_frames);
    
    first_image = imread(fullfile(image_folder, image_files(1).name));
    [height, width] = size(first_image);
    fprintf('Image dimensions: %d x %d pixels\n', height, width);
    
    image_sequence = zeros(height, width, num_frames, 'double');
    
    fprintf('Reading image sequence...\n');
    for i = 1:num_frames
        if mod(i, 100) == 0
            fprintf('Read %d/%d images\n', i, num_frames);
        end
        img = imread(fullfile(image_folder, image_files(i).name));
        image_sequence(:, :, i) = double(img);
    end
    fprintf('Image sequence reading complete\n');
    
    total_pixels = height * width;
    pixel_counter = 0;
    
    peak_data_table = table();
    
    for row = 1:height
        for col = 1:width
            pixel_id = sprintf('R%d_C%d', row-1, col-1);
            fprintf('Processing pixel %s (%d/%d)...\n', pixel_id, (row-1)*width + col, total_pixels);
            
            A = squeeze(image_sequence(row, col, :));
            
            frame_start = 200;
            frame_end = 2700;
            
            if num_frames < frame_end
                fprintf('Insufficient frames: only %d available\n', num_frames);
                continue;
            end
            
            if frame_end - frame_start + 1 < 100
                fprintf('Pixel %s has insufficient frames (%d), skipping.\n', pixel_id, frame_end - frame_start + 1);
                continue;
            end
            
            A = A(frame_start:frame_end);
            A = double(A);
            t = (1:length(A))';
            
            [p, y_fit, r_squared] = calculatePolyFitWithR2(t, A, 10);
            
            [success, cv_data, peak_data] = processSinglePixel(t, A, pixel_id, output_folder, r_squared);
            
            if success
                pixel_counter = pixel_counter + 1;
                if ~isempty(peak_data)
                    if isempty(peak_data_table)
                        peak_data_table = peak_data;
                    else
                        peak_data_table = [peak_data_table; peak_data];
                    end
                end
                fprintf('Pixel %s: processed successfully | R² = %.4f\n', pixel_id, r_squared);
            end
            
            if mod((row-1)*width + col, 50) == 0
                drawnow;
                pause(0.1);
            end
        end
    end
    
    if ~isempty(peak_data_table)
        peak_file = fullfile(output_folder, 'All_Pixels_Peak_Data.xlsx');
        writetable(peak_data_table, peak_file);
        fprintf('Peak data saved to: %s\n', peak_file);
    end
    
    msg = sprintf('Processing complete! Successfully processed %d/%d pixels\nResults saved to: %s', ...
                  pixel_counter, total_pixels, output_folder);
    msgbox(msg, 'Complete');
    
    if isfolder(output_folder)
        winopen(output_folder);
    end
end

function [p, y_fit, r_squared] = calculatePolyFitWithR2(t, A, poly_order)
    t = double(t);
    A = double(A);
    
    p = polyfit(t, A, poly_order);
    y_fit = polyval(p, t);
    
    y_mean = mean(A);
    ss_tot = sum((A - y_mean).^2);
    ss_res = sum((A - y_fit).^2);
    r_squared = 1 - (ss_res / ss_tot);
end

function [success, cv_data, peak_data] = processSinglePixel(t, A, pixel_id, output_folder, r_squared)
    success = false;
    cv_data = [];
    peak_data = table();
    
    try
        [p, y_fit, ~] = calculatePolyFitWithR2(t, A, 10);
        
        t_cycle = t;
        y_cycle = y_fit;
        
        y_diff = diff(y_cycle);
        
        t_kernel = t_cycle(1:end-1) + eps;
        kernel = 1./sqrt(pi * t_kernel);
        
        y_deconv = y_diff(:) ./ kernel(:);
        
        total_points = length(y_deconv);
        half_len = floor(total_points/2);
        
        y_forward = y_deconv(1:half_len);
        y_reverse = flip(y_deconv(end-half_len+1:end));
        
        potential = linspace(0, 1.4, half_len);
        
        [Ipa, Epa] = findpeaks(-y_forward, potential);
        [Ipc, Epc] = findpeaks(y_reverse, potential);
        
        if ~isempty(Ipa)
            [max_Ipa, ~] = max(Ipa);
        else
            max_Ipa = NaN;
        end
        
        if ~isempty(Ipc)
            [max_Ipc, ~] = max(Ipc);
        else
            max_Ipc = NaN;
        end
        
        peak_data = table({pixel_id}, max_Ipa, max_Ipc, ...
            'VariableNames', {'Pixel_ID', 'Ipa', 'Ipc'});
        
        fig = figure('Position', [100, 100, 800, 600], 'Visible', 'off');
        
        plot(potential, -y_forward, 'b', 'LineWidth', 2, 'DisplayName', 'Forward scan')
        hold on
        plot(potential, -y_reverse, 'r', 'LineWidth', 2, 'DisplayName', 'Reverse scan')
        
        xlabel('Potential (V)')
        ylabel('Current (a.u.)')
        title(sprintf('Pixel %s CV Curve', pixel_id))
        legend('Location', 'best')
        grid on
        
        png_file = fullfile(output_folder, sprintf('Pixel_%s_CV_Curve.png', pixel_id));
        fig_file = fullfile(output_folder, sprintf('Pixel_%s_CV_Curve.fig', pixel_id));
        
        saveas(fig, png_file);
        saveas(fig, fig_file);
        close(fig);
        
        cv_data = struct(...
            'pixel_id', pixel_id, ...
            'forward_scan', y_forward, ...
            'reverse_scan', y_reverse, ...
            'potential_forward', potential, ...
            'potential_reverse', potential, ...
            'deconv_signal', y_deconv, ...
            'r_squared', r_squared);
        
        success = true;
        
    catch ME
        fprintf('Error processing pixel %s: %s\n', pixel_id, ME.message);
        success = false;
    end
end