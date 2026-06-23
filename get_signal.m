%{ 

Get_signal.m
Created on 06/11/26 by Azra Karatan, Wang Lab

Gets the ΔF/F signal from the drawn cortex mask and plot changes in fluorescence over time relative to a baseline

%}

clear
clc
close all

%% Load TIFF movie

filename = '0020.tif'; 

info = imfinfo(filename); % imfinfo function returns information about the graphics file
N = numel(info); % numel returns the number of elements, n, in array A
N_TEST = 8; %TEST, CHANGE LATER TO N
Ny = info(1).Height; % 500
Nx = info(1).Width; % 500

movie = zeros(Ny, Nx, N_TEST, 'uint16'); % uint16 = 6-bit unsigned integer arrays (variables in MATLAB of data type (class) uint16 are stored as 2-byte (16-bit) unsigned integers)

for k = 1:N_TEST %TEST, change back to N
    movie(:,:,k) = imread(filename,k);% read the movie one frame at a time, becomes 500 × 500 × 2400
end 

fprintf('Loaded movie: %d \n', size(movie));

% Separate channels - [Green = hemodynamic(odds), Blue = GCamp(evens)]
movie_green = movie(:,:,1:2:end); %start from the 1st, get every other until the end
movie_blue  = movie(:,:,2:2:end); 


%% Apply Cortex Mask to Every Frame
load('brain_mask.mat')

movie_blue = double(movie_blue);
movie_green = double(movie_green);

for k = 1:size(movie_blue,3) %for each frame, set the outside pixels NaN
    frame_blue = movie_blue(:,:,k);
    frame_blue(~brain_mask) = NaN;
    movie_blue(:,:,k) = frame_blue;

    frame_green = movie_green(:,:,k);
    frame_green(~brain_mask) = NaN;
    movie_green(:,:,k) = frame_green;
end

disp("applied mask")

%sanity checks to make sure mask is applied correct and to every 1200 frame:
nNaN_frame4 = sum(isnan(movie_blue(:,:,4)),'all');
fprintf('Number of NaN pixels in movie_blue frame 4: %d\n', nNaN_frame4);
check = movie_blue(: ,: , 1);
openvar('check')


%% Background Removal
% Preprocessing params (matches Analyze_one_WF_run hard-coded values)
background_blue  = 1650; 
background_green = 1600;

% Remove background F in the absence of GCAMP
function [movie_blue, movie_green] = removeBackground(movie_blue, movie_green, background_blue, background_green)
    movie_blue = movie_blue - background_blue;
    movie_blue(movie_blue < 0) = 1; %making sure everything is above 0 so the /F0 division doesn't crash
    movie_green = movie_green - background_green;
    movie_green(movie_green < 0) = 1;
end 

[movie_blue, movie_green] = removeBackground(movie_blue, movie_green, background_blue, background_green);
disp("background removed")

%% Calculate the Baseline F0  
function F0 = f0_calculation(movie, method, percentile, window_size)

    [Ny, Nx, Nt] = size(movie); % 500x500x1200 for us
    %F0(y,x,t) will be the baseline fluorescence for pixel (y,x) at time t

    disp("begin calculating f0")

    switch method
        case 'sliding_window'
            % For each frame, F0 was estimated from a 10 s sliding window (200 frames) 
            % as the mean of the lowest 50% of fluorescence values in that window.

            %flatten movie to [nPixels x Nt]
            movie_2d = reshape(movie, [], Nt); %[all pixels x time] -- movie_2d becomes = 250000 rows (one for each pixel) x 1200 columns (one for each frame), easier to work with
            fprintf('size of movie_2d: %d\n', size(movie_2d))

            % identify pixels inside the mask
            valid_pixels = ~isnan(movie_2d(:,1));

            % keep only valid cortex pixels
            movie_valid = movie_2d(valid_pixels, :);
            fprintf('size of movie_valid: %d\n', size(movie_valid))

            % allocate space for valid pixels only
            F0_valid = NaN(size(movie_valid));

            % For an even window size (e.g. 200), split around current frame: so 99 before + current frame + 100 after = total window size 200
            left_half_size  = floor((window_size - 1)/2); %floor rounds each element of X to the nearest integer less than or equal to that element.
            right_half_size = ceil((window_size - 1)/2); %ceil rounds each element of X to the nearest integer greater than or equal to that element.
        
            disp("just before the loop")

            for t = 1:Nt %loop over every frame/timepoint (loops through every column)
                 %define window bounds around current frame
                 t_start = max(1, t - left_half_size);%start of the window
                 t_end   = min(Nt, t + right_half_size);%end of the window
                      
                 %extract all pixel traces within this time window: size =[nPixels x window_length] so 200 of pixel1 through this window
                 window_vals = movie_valid(:, t_start:t_end);
              
                 %sort each pixel's values across time within the window
                 window_vals_sorted = sort(window_vals, 2, 'ascend'); %2 means dimension 2 (so it sorts along the columns)
                           
                 %only keep the lower half
                 n_vals = size(window_vals_sorted, 2); %total # of items in the list
                 n_keep = max(1, floor(percentile/100 * n_vals)); %# of values to keep (half of the list) 
                 lowest_vals = window_vals_sorted(:, 1:n_keep); %store the first half
                       
                 % mean of those lowest values = F0 for frame t
                 F0_valid(:, t) = mean(lowest_vals, 2, 'omitnan');
            end
       
            disp("loop is done")
            % reshape back to movie form
        
            F0_2d = NaN(size(movie_2d));
            F0_2d(valid_pixels, :) = F0_valid;

            F0 = reshape(F0_2d, Ny, Nx, Nt); %F0(y,x,t) = baseline fluorescence for pixel (y,x) at frame t, for pxeil-wise df/f calculation
            disp("F0 are calculated by sliding window method and movie reshaped back to original size succesfully")
       
        otherwise
            error('error: unknown F0 calculation method.')
     end
end

F0 = f0_calculation(movie_blue, 'sliding_window', 50, 2); %TEST -- CHANGE TO 200 FOR ACTUAL RUN.  

%% Convert to dF/F 
% df_f = (movie_blue - F0) ./ F0;

%{ 
-- eva's version for dff conversion:

function [df_blue, df_green] = convertToDFF(movie_blue, movie_green, Frame_stimOn, brain_mask, Nimages, Nstim) %said to take out frame_stin0n
    baseline_blue = squeeze(mean(mean(movie_blue(:,:,1:Frame_stimOn-1,:),3),4)) .* brain_mask;
    baseline_blue = repmat(baseline_blue,[1,1,Nimages,Nstim]);
    
    baseline_green = squeeze(mean(mean(movie_green(:,:,1:Frame_stimOn-1,:),3),4)) .* brain_mask;
    baseline_green = repmat(baseline_green,[1,1,Nimages,Nstim]);
    
    df_blue = movie_blue ./ baseline_blue -1;
    df_green = movie_green ./ baseline_green -1;

end

[df_blue, df_green] = convertToDFF(movie_blue, movie_green, Frame_stimOn, brain_mask, Nimages, Nstim); %--same with Nstim? and is Nimages here the same as my N variable?
% [df_blue_linear, df_green_linear] = convertToDFF_linear(movie_blue, movie_green, Frame_stimOn, brain_mask, Nimages, Nstim);

-----

Lowpass filter 
    function [df_blue, df_green] = applyLowPassFilter(df_blue, df_green, FrameRate_PCO, cutoff, Ny, Nx, Nstim, Nimages)
    % Ensure that Ny, Nx, Nstim, and Nimages are defined before calling this function

    % Reshape the matrices
    df_blue = reshape(df_blue,[Ny,Nx,Nstim*Nimages]);
    df_green = reshape(df_green,[Ny,Nx,Nstim*Nimages]);
    
    % Apply low-pass filter
    df_blue = Low_pass_TC(df_blue,FrameRate_PCO,cutoff);
    df_green = Low_pass_TC(df_green,FrameRate_PCO,cutoff);
end

[df_blue, df_green] = applyLowPassFilter(df_blue, df_green, FrameRate_PCO, parameters(4), Ny, Nx, Nstim, Nimages);

% Regress and subtract
function [df_movie,df_green_reshaped] = regressAndSubtract(df_blue, df_green, Nstim, Nimages, Ny, Nx)
    % Calculate the mean of df_blue and df_green
    mean_df_blue = repmat(mean(df_blue, 3), [1, 1, Nstim * Nimages]);
    mean_df_green = repmat(mean(df_green, 3), [1, 1, Nstim * Nimages]);

    % Calculate slope using regression
    slope = sum((df_blue - mean_df_blue) .* (df_green - mean_df_green), 3) ...
            ./ sum((df_green - mean_df_green).^2, 3);
    slope = repmat(slope, [1, 1, Nimages * Nstim]);

    % Subtract rescaled green from blue
    df_movie = df_blue - df_green .* slope;

    % Reshape df_movie to original dimensions
    df_movie = reshape(df_movie, [Ny, Nx, Nimages, Nstim]);
    df_green_reshaped = reshape(df_green, [Ny, Nx, Nimages, Nstim]);
end

[df_movie,df_green_reshaped]  = regressAndSubtract(df_blue, df_green, Nstim, Nimages, Ny, Nx);
%df_movie = HemoCorrection(df_blue,df_green,movie_green(:,:,1:Frame_stimOn-1,:),5)



%% Compute dF/F movie

df_f = zeros(size(movie_blue)); %create an empty array of this size (to be filled in in the next step)

for k = 1:size(movie_blue,3) %if size is 500 x 500 x 1200, this gets the 1200 which is the 3rd item - so 1 frame at a time for 1:1200
    df_f(:,:,k) = (movie_blue(:,:,k) - F0) ./ F0; %calculates df_f for every frame
end % this "./" syntax means element-wise division, so every pixel will get divided by its matching pixel

disp('dF/F movie computed succesfully.') %resulting df_f movie contains percent fluorescence change over time for each pixel




%% Compute mean cortex signal

mean_signal = zeros(1,size(df_f,3)); % create empty space to store the mean signal

for k = 1:size(df_f,3) %for 1200: processes one frame at a time
    frame = df_f(:,:,k); 
    mean_signal(k) = mean(frame(:),'omitnan'); %converts the pixels of the frame into a single linear list because it's easier to calculate, then computes the mean of this frame
end %result gives you the average ΔF/F across the whole cortex at this frame

%% Plot the dF/F signal

figure
plot(mean_signal,'LineWidth',2)

xlabel('Frame Number')
ylabel('\DeltaF/F')
title('Mean Cortex dF/F Signal')



%% Load TIFF movie
%% Setup the Mask and Compute the Corrected Delta F/F
    
% Place the mask: only takes the cortex and trims imaging window, so the analysis only happens inside the brain and everything else is excluded (set to NaN, so they get ignored). 
brain_mask = applyCortexMask(analysis_dir, mouse_name, side, run_date, run_name, movie_blue); %(where output lives, name,L/R hemisphere, ..)

[movie_blue, movie_green] = removeBackground(movie_blue, movie_green, background_blue, background_green);% we want true fluorescence only, so this removes the background signals in the image by subtracting the background from its respective movie
[df_blue, df_green] = convertToDFF(movie_blue, movie_green, Frame_stimOn, brain_mask, Nimages, Nstim); % Conversion to df/f -- convert pixel brightness into percent change relative to baseline fluorescence
% alt: [df_blue_linear, df_green_linear] = convertToDFF_linear(movie_blue, movie_green, Frame_stimOn, brain_mask, Nimages, Nstim);

[df_blue, df_green] = applyLowPassFilter(df_blue, df_green, FrameRate_PCO, parameters(4), Ny, Nx, Nstim, Nimages); % Lowpass filter - removes noise
[df_movie,df_green_reshaped]  = regressAndSubtract(df_blue, df_green, Nstim, Nimages, Ny, Nx); % Hemodynamic correction - subtracts the blood flow-related changes in the light
% alt: df_movie = HemoCorrection(df_blue,df_green,movie_green(:,:,1:Frame_stimOn-1,:),5)

%} 

%}
