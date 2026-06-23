%{ 

Seizure_analysis1.m
Created on 06/082/26 by Azra Karatan, Wang Lab

The purpose of this file is to create a manual cortex mask of uploaded brain image tif file. 

[Green = hemodynamic, Blue = GCamp]

%}

%{ 
Display the first frame of the movie

info = imfinfo('0020.tif'); %open the tif
framecount = numel(info);  %count the frames, it is 2400
fprintf('framecount: %d\n', framecount);

frame1 = imread('0020.tif',1); %load the image
figure 	%open the image in a separate window in matlab
imagesc(frame1) %only open frame 1 out of 2400
axis image %the image is 500x500
colorbar
colormap gray 
%} 


%% Load the entire TIFF movie
filename = '0020.tif'; 
info = imfinfo(filename);
N = numel(info); %number of elements-returns the total # of elements contained

Ny = info(1).Height;
Nx = info(1).Width;

movie = zeros(Ny, Nx, N, 'uint16');

for k = 1:N
    movie(:,:,k) = imread(filename,k);
end

fprintf('Loaded movie: %d x %d x %d\n', Ny, Nx, N); %500x500x2400

%% Separate the channels
green  = movie(:,:,1:2:end); %counts the odds, starting from 1: #1200
blue = movie(:,:,2:2:end); %counts the evens, starting from 2: #1200

figure
subplot(1,2,1)
imagesc(green(:,:,1))
axis image
title('Green Channel (Hemodynamic Response)')

subplot(1,2,2)
imagesc(blue(:,:,1))
axis image
title('Blue Channel (GCamp Fluorescence)')

colormap gray

%% Manually make the cortex mask for frame 1

make_new_mask = true;

% First, check if mask file already exists
if isfile('brain_mask.mat')
    choice = questdlg( ...
        'brain_mask.mat already exists. Overwrite it?', ...
        'Overwrite Mask?', ...
        'Yes', 'No', 'No');

    if strcmp(choice, 'No')
        disp('Mask was not overwritten. Previous exists.')

        choice = questdlg( ...
            'Open existing mask?', ...
            'Existing Mask', ...
            'Yes', 'No', 'Yes');

        if strcmp(choice, 'Yes')
            
            %{ 
            this was for a check, not used now

            S = load('brain_mask.mat');

            if isfield(S,'brain_mask')
                brain_mask = S.brain_mask;

            elseif isfield(S,'mask_global')
                brain_mask = S.mask_global;
                save('brain_mask.mat','brain_mask')   % clean up old variable name
                disp('Old mask_global found and converted to brain_mask.')

            else
                error('brain_mask.mat contains neither brain_mask nor mask_global.')
            end
            %}

            disp('Loaded existing mask.')
            make_new_mask = false;
        end
    end
end

% Only if needed, make new mask and save it
if make_new_mask

    figure
    imagesc(green(:,:,1))
    axis image
    colormap gray
    title('Draw cortex boundary')

    h = drawfreehand; %creates a Freehand object and enables interactive drawing of a circular region-of-interest (ROI) on the current axes.

    % Save the mask once it's done
    uiwait(msgbox(sprintf([ ...
        'Adjust the mask if needed, then click OK to save.\n' ...
        '(To adjust: double-click to add waypoints and drag)'])));

    % Create and display mask :)
    brain_mask = createMask(h);
    save('brain_mask.mat','brain_mask')
    disp('Mask saved!')
end

% Load the mask and display, if it was created correctly
if isfile('brain_mask.mat')
    load('brain_mask.mat')

    %Show the brain mask
    %on its own
    figure
    imagesc(brain_mask)
    axis image
    colormap gray
    title('Brain Mask')
        
    %overlayed on the original brain image (green, blue)
    figure
    title('Brain Masks')
    subplot(1,2,1)
    imagesc(green(:,:,1))
    axis image
    subtitle('Green Channel (Hemodynamic Response)')
    hold on
    contour(brain_mask, [1 1], 'r', 'LineWidth', 2)

    subplot(1,2,2)
    imagesc(blue(:,:,1))
    axis image
    subtitle('Blue Channel (GCamp Fluorescence)')
    
    colormap gray
    hold on %hold on allows you to add a second line plot without deleting the existing line plot, and hold off deletes the first plot when you make a new one
    contour(brain_mask, [1 1], 'r', 'LineWidth', 2)

    %open the mat file itself for inspection/control --- the 000011110000's map ---
    %openvar('brain_mask')
else
    error('Error: brain_mask.mat was not found. (probably not saved correctly)')
end
