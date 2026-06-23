function [allen_boundaries,allen_transformed,mask_global,BehaviorFrames,trials_with_wipes,trials_wo_wipes,trials_with_groom,N_with_wipes,N_wo_wipes,N_with_groom,tform,df_map,df_map_with_wipes,df_map_wo_wipes,TC,TC_with_wipes,TC_wo_wipes] =...
    Analyze_one_WF_run(line_num,data_dir,behavior_dir,analysis_dir,mouse_name,run_date,run_name,Nstim,Nx,Ny,FrameRate_PCO,FrameRate_Basler,trial_duration,delay,stim_duration,stim_type,stim_param,DLC_resnet,side,parameters,save_videos)

% This function analyze single wide field runs

if nargin < 21 || isempty(save_videos)
    save_videos = false;
end

tic;
load('allen_atlas_reordered_resampled.mat');
barrel_idx=[6,12,20]; %A2,C2,E2
barrel_roi_atlas_idx = barrel_idx(:); % Knutsen / allen_BC_right row index; same order as dim 2 of TC_all_trials_barrels
barrel_roi_names = {'A2','C2','E2'};
barrel_roi_tc_desc = ['TC_all_trials_barrels is [Nimages x nBarrelROIs x Nstim]. Index b matches barrel_roi_names{b} / barrel_roi_atlas_idx(b). ', ...
    'Each value is mean df_movie (same ΔF/F as Allen) over the 5×5 patch; mask uses NaN outside ROI so nanmean is over patch pixels only.'];
% barrel_cortex_pixel_timecourse.mat: only Allen ROI S1bf_R (standard 56-ROI index 56 when names align).
barrel_cortex_pixel_roi_name = 'S1bf_R';           % case-insensitive match to cell_areanames{i}
barrel_cortex_pixel_roi_linear_index = 56;       % if name not found, use area_num(this index)
max_abs=parameters(1);
min_abs=parameters(2);
BehaviorFrames = [];
trials_with_wipes = [];
trials_wo_wipes = [];
trials_with_groom = [];
N_with_wipes = [];
N_wo_wipes = [];
N_with_groom = [];
N_with_whisk = [];
PCO_Frame_start_wipe = [];
Frame_start_wipe = [];
mean_df_movie = [];
mean_df_movie_with_wipes = [];
mean_df_movie_wo_wipes = [];
mean_df_green_with_wipes = [];
mean_df_green_wo_wipes = [];
TC_with_wipes_green = [];
TC_wo_wipes_green = [];
df_map = [];
df_map_with_wipes = [];
df_map_wo_wipes = [];
TC = [];
TC_with_wipes = [];
TC_wo_wipes = [];
TC_motion_trigged = [];
median_shift = 0;
max_shift = 0;
TC_all_trials = [];
TC_all_trials_with_wipes = [];
TC_all_trials_wo_wipes = [];
df_movie_with_wipes= [];
df_movie_wo_wipes= [];
TC_all_trials_barrels = [];
TC_all_trials_barrels_with_wipes = [];
TC_all_trials_barrels_wo_wipes = [];

background_blue = 1650;
background_green = 1600;


n_maps = 9;
save_data = 1;
close_figs = 1;  % set 0 while debugging so "close all" at end does not wipe figures

% Behavior files (DLC csv, labeled video, results_behavior.mat) are usually stored once per session
% under ...\Behavior\...\run1\. Von Fray (and thermal stim_*) call this function with nested
% run_name (e.g. run1\0.02); DLC is not duplicated under each gram folder, so look up session root.
rn_wf = char(string(run_name));
dlc_fn = [char(string(DLC_resnet)) '.csv'];
dlc_nested = fullfile(analysis_dir, 'Behavior', mouse_name, side, run_date, rn_wf, dlc_fn);
tok = regexp(rn_wf, '[^\\/]+', 'match');
if numel(tok) >= 2
    rn_beh_parent = tok{1};
    dlc_parent = fullfile(analysis_dir, 'Behavior', mouse_name, side, run_date, rn_beh_parent, dlc_fn);
else
    rn_beh_parent = rn_wf;
    dlc_parent = dlc_nested;
end
if isfile(dlc_nested)
    behavior_run_name = rn_wf;
    DLC_csv_name = dlc_nested;
elseif isfile(dlc_parent)
    behavior_run_name = rn_beh_parent;
    DLC_csv_name = dlc_parent;
else
    behavior_run_name = rn_wf;
    DLC_csv_name = dlc_nested;
end

mkdir([analysis_dir 'Wide_Field\' mouse_name '\' side '\' run_date '\' run_name]);

Nimages = FrameRate_PCO * trial_duration;

Frame_stimOn = floor(delay*FrameRate_PCO);
Frame_stimOff = floor((delay+stim_duration)*FrameRate_PCO);

time_vec = 0:1/FrameRate_PCO:(Nimages-1)/FrameRate_PCO;

% Run lists may use mixed case or trailing spaces; behavior branches must still match
stim_type_norm = lower(strtrim(char(string(stim_type))));

if strcmp(stim_param(end),'R')    
    stim_L = 0;    
elseif strcmp(stim_param(end),'L')
    stim_L = 1;
else
    stim_L = 0;     
end

% Load and analyze behavior

if strcmp(stim_type_norm,'piezo')
    
    behavior_analysis = 0;
    
elseif strcmp(stim_type_norm, 'vonfray') || strcmp(stim_type_norm, 'pole')
    
    BeResults_FilePath = fullfile(analysis_dir, 'Behavior', mouse_name, side, run_date, behavior_run_name, 'results_behavior.mat');
    labeled_mp4 = fullfile(analysis_dir, 'Behavior', mouse_name, side, run_date, behavior_run_name, [char(string(DLC_resnet)) '_labeled.mp4']);
    if ~isfile(DLC_csv_name)
        behavior_analysis = 0;
        warning('Analyze_one_WF_run:NoDLCCsv', ...
            'Behavior/wipes disabled: DLC csv not found (nested + session):\n  %s\n  %s', dlc_nested, dlc_parent);
    elseif ~isfile(BeResults_FilePath)
        behavior_analysis = 0;
        warning('Analyze_one_WF_run:NoBehaviorMat', ...
            'Behavior/wipes disabled: results_behavior.mat missing — run behavior first:\n  %s', BeResults_FilePath);
    elseif ~isfile(labeled_mp4)
        behavior_analysis = 0;
        warning('Analyze_one_WF_run:NoLabeledVideo', ...
            'Behavior/wipes disabled: labeled video missing:\n  %s', labeled_mp4);
    else
        behavior_analysis = 1;
        load(BeResults_FilePath);
        Frame_start_wipe = round(PCO_Frame_start_wipe / (FrameRate_Basler/FrameRate_PCO));
        v = VideoReader(labeled_mp4);
        mouse_portrait = readFrame(v);
        movie_Be = read(v);
        movie_Be = squeeze(single(movie_Be(37:536,221:720,1:3,1:(FrameRate_Basler*2/FrameRate_PCO):(FrameRate_Basler*trial_duration*Nstim))));
        movie_Be = reshape(movie_Be,[Ny,Nx,3,Nimages/2,Nstim]);
        if stim_L == 1
            mouse_portrait = flip(mouse_portrait,2);
            movie_Be = flip(movie_Be,2);
        end
    end
       
elseif strcmp(stim_type_norm,'thermal')
    
    BeResults_FilePath = fullfile(analysis_dir, 'Behavior', mouse_name, side, run_date, behavior_run_name, 'results_behavior.mat');
    labeled_mp4 = fullfile(analysis_dir, 'Behavior', mouse_name, side, run_date, behavior_run_name, [char(string(DLC_resnet)) '_labeled.mp4']);
    if ~isfile(DLC_csv_name)
        behavior_analysis = 0;
        warning('Analyze_one_WF_run:NoDLCCsv', ...
            'Behavior/wipes disabled: DLC csv not found:\n  %s\n  %s', dlc_nested, dlc_parent);
    elseif ~isfile(BeResults_FilePath)
        behavior_analysis = 0;
        warning('Analyze_one_WF_run:NoBehaviorMat', ...
            'Behavior/wipes disabled: results_behavior.mat missing:\n  %s', BeResults_FilePath);
    elseif ~isfile(labeled_mp4)
        behavior_analysis = 0;
        warning('Analyze_one_WF_run:NoLabeledVideo', ...
            'Behavior/wipes disabled: labeled video missing:\n  %s', labeled_mp4);
    else
        behavior_analysis = 2;
        load(BeResults_FilePath);
        Frame_start_wipe = round(PCO_Frame_start_wipe / (FrameRate_Basler/FrameRate_PCO));
        v = VideoReader(labeled_mp4);
        mouse_portrait = readFrame(v);
        movie_Be = read(v);
        if FrameRate_Basler*trial_duration*Nstim ~= v.NumFrames
            warning('video contains fewer frames than expected, dropping last trial!');
            Nstim = Nstim -1
        end
        movie_Be = squeeze(single(movie_Be(37:536,171:670,1:3,1:(FrameRate_Basler*2/FrameRate_PCO):(FrameRate_Basler*trial_duration*Nstim))));
        movie_Be = reshape(movie_Be,[Ny,Nx,3,Nimages/2,Nstim]);
        if stim_L == 1
            mouse_portrait = flip(mouse_portrait,2);
            movie_Be = flip(movie_Be,2);
        end
    end
           
else
    
    warning('Analyze_one_WF_run:NoBehaviorStimType', ...
        'stim_type "%s" not recognized for behavior (expected piezo, vonfray, pole, thermal).', stim_type_norm);
    behavior_analysis = 0;
    
end

% Load widefield images
pco_fileList = find_tif_files([data_dir mouse_name '\' side '\' run_date '\' run_name ]);
for j=1:length(pco_fileList)
    pco_patterns= regexp(pco_fileList{j}, ['#' '|' '.tif'], 'split');
    pco_num_list(j)=string(pco_patterns{2});
end


movie_WF = zeros(Ny,Nx,Nimages,Nstim); %4D array -> 500x500x160x10 (pixel x pixel x # images in stack x # trials)

if length(pco_num_list)~=Nstim
    warning("pco files do not match the expected number of trials")
end

% Trials whose PCO stack has the wrong frame count (too short, empty, or missing) are flagged
% and dropped below via the existing bad_trials handling. Stacks that are longer than Nimages
% are truncated (real data kept). NEVER pad with a repeated frame: it biases baseline/df_f and
% contaminates mean_df_movie, which can change the plot scale by an order of magnitude.
bad_trials_par = false(Nstim,1);

parfor i = 1:Nstim
    Stackname_pattern = [data_dir mouse_name '\' side '\' run_date '\' run_name '\pco*#' char(pco_num_list(i)) '.tif'];
    match_stackname_files = dir(Stackname_pattern);
    if isempty(match_stackname_files)
        warning('Analyze_one_WF_run:PCOStackMissing', ...
            'Trial %d: no TIFF matched %s — dropping trial.', i, Stackname_pattern);
        bad_trials_par(i) = true;
        continue
    end

    % Extract the file name of the first matching file
    Stackname = fullfile( match_stackname_files(1).folder,  match_stackname_files(1).name);


    [MyMovie,~] = LoadStack(Stackname,'uint16');
    nFr = size(MyMovie, 3);
    if nFr == Nimages
        movie_WF(:,:,:,i) = MyMovie;
    elseif nFr > Nimages
        warning('Analyze_one_WF_run:PCOFramesLong', ...
            'Trial %d has %d frames > Nimages=%d; truncating to Nimages.\n  %s', ...
            i, nFr, Nimages, Stackname);
        movie_WF(:,:,:,i) = MyMovie(:,:,1:Nimages);
    else
        warning('Analyze_one_WF_run:PCOFramesShort', ...
            ['Trial %d has %d frames < Nimages=%d — dropping trial (will NOT pad).\n  %s\n' ...
            'Check acquisition/export for this level or adjust trial_duration/FrameRate_PCO in the run list.'], ...
            i, nFr, Nimages, Stackname);
        bad_trials_par(i) = true;
    end
end

%the following removes the first 2 frames of every tiff
movie_WF(:,:,1,:) = movie_WF(:,:,3,:);      % because hardware messes up the first 2 frames :(
movie_WF(:,:,2,:) = movie_WF(:,:,4,:);

% find missed frames trials
bad_trials = bad_trials_par;
if any(bad_trials)
    fprintf(1, 'Analyze_one_WF_run: dropping %d/%d trials with missing or short PCO stacks.\n', ...
        sum(bad_trials), Nstim);
end

 movie_WF(:,:,:,bad_trials) = [];

if behavior_analysis > 0

    movie_Be(:,:,:,:,bad_trials) = [];
    
    trials_wo_wipes(bad_trials) = [];
    N_wo_wipes = sum(trials_wo_wipes);
    
    trials_with_wipes(bad_trials) = [];
    N_with_wipes = sum(trials_with_wipes);
    
    trials_with_groom(bad_trials) = [];
    N_with_groom = sum(trials_with_groom);
    
    Frame_start_wipe(bad_trials) = [];

    ind_bad = [];
    for i = 1:length(bad_trials)
        if bad_trials(i)
            ind_bad = [ind_bad [((i-1)*Nimages/2+1):(i*Nimages/2)]];
        end
    end
%     lefthand_x(ind_bad) = [];
%     lefthand_y(ind_bad) = [];
    righthand_x(ind_bad) = [];
    righthand_y(ind_bad) = [];
end

Nstim = Nstim - sum(bad_trials);

% spatial smooting with Gaussian kernel
movie_WF = spatialSmooth(movie_WF, Ny, Nx, Nimages, Nstim); %%xxxxxxxx

% separate blue and green movies
[movie_blue, movie_green] = separateColorChannels(movie_WF);

%halve sampling rates and # iamges since the channels are split
Frame_stimOn = round(Frame_stimOn/2);
Frame_stimOff = round(Frame_stimOff/2);
Frame_start_wipe = round(Frame_start_wipe/2);
FrameRate_PCO = FrameRate_PCO/2;
Nimages = Nimages/2;
time_vec = 0:1/FrameRate_PCO:(Nimages-1)/FrameRate_PCO;

% Place mask on cortex and trim imaging window to keep only cortex (the rest = NaN)
mask_global = applyCortexMask(analysis_dir, mouse_name, side, run_date, run_name, movie_blue);

% remove background F in the absence of GCAMP
[movie_blue, movie_green] = removeBackground(movie_blue, movie_green, background_blue, background_green);

% Convert to df/f
[df_blue, df_green] = convertToDFF(movie_blue, movie_green, Frame_stimOn, mask_global, Nimages, Nstim);
% [df_blue_linear, df_green_linear] = convertToDFF_linear(movie_blue, movie_green, Frame_stimOn, mask_global, Nimages, Nstim);


% Lowpass filter
[df_blue, df_green] = applyLowPassFilter(df_blue, df_green, FrameRate_PCO, parameters(4), Ny, Nx, Nstim, Nimages);

[df_movie,df_green_reshaped]  = regressAndSubtract(df_blue, df_green, Nstim, Nimages, Ny, Nx);
%df_movie = HemoCorrection(df_blue,df_green,movie_green(:,:,1:Frame_stimOn-1,:),5)


[tfPath_before_align, ~] = wf_widefield_run_paths(analysis_dir, mouse_name, side, run_date, run_name, 'Transform.mat');
[tform, allen_boundaries, allen_transformed, cell_areanames, area_num, n_areas, fig_barrel_map] = alignAtlasToPCO(movie_WF, mask_global, analysis_dir, mouse_name, side, run_date, run_name);
[all_allen_x_transformed_right, all_allen_y_transformed_right] = transformPointsForward(tform, allen_BC_right(:,1), allen_BC_right(:,2));

% Define the size of the square (must match df_movie / allen_transformed Ny x Nx)
square_size = 5;
half_size = floor(square_size / 2);
mask_barrels = zeros(Ny, Nx, numel(barrel_idx));
% Loop through each point and set the corresponding 5x5 square to 1
for i = 1:size(all_allen_x_transformed_right(barrel_idx), 1)
    x = all_allen_x_transformed_right(barrel_idx(i));
    y = all_allen_y_transformed_right(barrel_idx(i));
    
    % Determine the range for the square, ensuring it stays within matrix bounds
    x_start = max(1, x - half_size);
    x_end = min(Nx, x + half_size);
    y_start = max(1, y - half_size);
    y_end = min(Ny, y + half_size);
    
    % Set the corresponding elements to 1
    mask_barrels(y_start:y_end, x_start:x_end,i) = 1;
end
% Match Allen ROI weighting: outside ROI must be NaN so nanmean averages only
% the 5×5 patch, not the full FOV (zeros would dilute amplitude by n_roi/n_pixels).
for ib = 1:size(mask_barrels, 3)
    mslice = single(mask_barrels(:, :, ib));
    mslice(mslice == 0) = NaN;
    mask_barrels(:, :, ib) = mslice;
end
% New manual registration only: save barrel-map figure under this stim folder (Transform.mat is saved in alignAtlasToPCO).
if isempty(tfPath_before_align) && ~isempty(fig_barrel_map) && isgraphics(fig_barrel_map)
    savefig(fig_barrel_map,[analysis_dir 'Wide_Field\' mouse_name '\' side '\' run_date '\' run_name '\fig_barrel_map.fig']);
    saveas(fig_barrel_map,[analysis_dir 'Wide_Field\' mouse_name '\' side '\' run_date '\' run_name '\fig_barrel_map.pdf']);
end

colormat_time =[[0 0 0];jet(n_areas/2)]; 
colormat_space = [[0 0 0];parula(255)];

% Compute df map at different times (every 0.5s between stim start and 2s after)
% If behavior has ben analyzed, separate trials with and without wipes

df_movie(df_movie == 0 | ~isfinite(df_movie)) = 0;
% divide 2s post stimulus into n_maps-1 

n_integration_frames = floor(2*FrameRate_PCO/(n_maps-1));

mean_df_movie = squeeze(mean(df_movie,4));
mean_df_green= squeeze(mean(df_green_reshaped,4));

fig_df_maps = figure('Position',[0 550 1500 400]);
set(fig_df_maps,'PaperOrientation','landscape');
set(fig_df_maps,'PaperType','uslegal');

% Color for the atlas overlay (same in every panel, regardless of data scale)
atlas_overlay_rgb = [0 0 0]; % change to [1 1 1] for white, [0 0 0.3] for dark blue
boundary_mask = logical(allen_boundaries);

if behavior_analysis == 0

    for i = 1:n_maps
        df_map{i} = squeeze(mean(mean_df_movie(:,:,(Frame_stimOn+(i-2)*n_integration_frames):(Frame_stimOn+(i-1)*n_integration_frames)),3));
    end

    % One color range over the whole row (all snapshots) with ±10% padding
    all_vals = cellfun(@(m) m(~boundary_mask), df_map, 'UniformOutput', false);
    all_vals = vertcat(all_vals{:});
    [min_plot, max_plot] = df_map_clim(all_vals);

    for i = 1:n_maps
        subplot(1,n_maps,i);
        imagesc(df_map{i}, [min_plot, max_plot]);
        hold on;
        overlay_atlas_boundaries_solid(boundary_mask, atlas_overlay_rgb);
        hold off;
        text(100,-100,[num2str((i-2)*n_integration_frames/FrameRate_PCO) '-' num2str((i-1)*n_integration_frames/FrameRate_PCO) 's'],'Color','white');
        axis off;
        axis equal;
        axis tight;
        colormap(colormat_space);
        if i == n_maps-1
            colorbar('Position',[0.91 0.65 0.005 0.2]);
        end
    end

else

    idx_with = find(trials_with_wipes);
    idx_without = find(trials_wo_wipes);
    mean_df_movie_with_wipes = squeeze(mean(df_movie(:,:,:,idx_with),4));
    mean_df_movie_wo_wipes = squeeze(mean(df_movie(:,:,:,idx_without),4));
    mean_df_green_with_wipes = squeeze(mean(df_green_reshaped(:,:,:,idx_with),4));
    mean_df_green_wo_wipes = squeeze(mean(df_green_reshaped(:,:,:,idx_without),4));

    if ~isempty(idx_with)
        df_movie_with_wipes = df_movie(:,:,:,idx_with);
    end
    if ~isempty(idx_without)
        df_movie_wo_wipes = df_movie(:,:,:,idx_without);
    end


    for i = 1:n_maps
        df_map_with_wipes{i} = squeeze(mean(mean_df_movie_with_wipes(:,:,(Frame_stimOn+(i-2)*n_integration_frames):(Frame_stimOn+(i-1)*n_integration_frames)),3));
        df_map_wo_wipes{i} = squeeze(mean(mean_df_movie_wo_wipes(:,:,(Frame_stimOn+(i-2)*n_integration_frames):(Frame_stimOn+(i-1)*n_integration_frames)),3));
    end

    % Row-wise color ranges (different magnitudes for wipes vs no-wipes).
    if N_with_wipes > 0
        vw_vals = cellfun(@(m) m(~boundary_mask), df_map_with_wipes, 'UniformOutput', false);
        [min_plot_w, max_plot_w] = df_map_clim(vertcat(vw_vals{:}));
    end
    if N_wo_wipes > 0
        vo_vals = cellfun(@(m) m(~boundary_mask), df_map_wo_wipes, 'UniformOutput', false);
        [min_plot_o, max_plot_o] = df_map_clim(vertcat(vo_vals{:}));
    end

    for i = 1:n_maps

        if N_with_wipes>0
            subplot(2,n_maps,i);
            imagesc(df_map_with_wipes{i}, [min_plot_w, max_plot_w]);
            hold on;
            overlay_atlas_boundaries_solid(boundary_mask, atlas_overlay_rgb);
            hold off;
            text(100,20,[num2str((i-2)*n_integration_frames/FrameRate_PCO) ' - ' num2str((i-1)*n_integration_frames/FrameRate_PCO) 's'],'Color','white');
            axis off;
            axis equal;
            axis tight;
            if i == n_maps
                colorbar('Position',[0.91 0.65 0.005 0.2]);
            end
        end

        if N_wo_wipes>0
            subplot(2,n_maps,i+n_maps);
            imagesc(df_map_wo_wipes{i}, [min_plot_o, max_plot_o]);
            hold on;
            overlay_atlas_boundaries_solid(boundary_mask, atlas_overlay_rgb);
            hold off;
            text(100,20,[num2str((i-2)*n_integration_frames/FrameRate_PCO) ' - ' num2str((i-1)*n_integration_frames/FrameRate_PCO) 's'],'Color','white');
            axis off;
            axis equal;
            axis tight;
            if i == n_maps
                colorbar('Position',[0.91 0.2 0.005 0.2]);
            end
        end
    end
end


% Compute TC in each Allen atlas region (separately for wipes / no wipes if available )

fig_TCs_all = figure('Position',[500 200 400 800]);
set(fig_TCs_all,'PaperOrientation','portrait');
set(fig_TCs_all,'PaperType','uslegal');
    
line_offset = max(mean_df_movie(:)/4);

if behavior_analysis == 0
    
    LH = subplot(1,2,1);
    RH = subplot(1,2,2);
   
    for i = 1:n_areas
        
        mask_temp = single(allen_transformed == area_num(i));
        mask_temp(mask_temp == 0) = NaN;
        
        TC(:,i) = squeeze(nanmean(nanmean(mean_df_movie.*repmat(mask_temp,[1,1,Nimages]))));

        TC_all_trials(:,i,:)= squeeze(nanmean(nanmean(df_movie.*repmat(mask_temp,[1,1,Nimages]))));


        if i <= n_areas/2
            axes(LH);
            plot(time_vec,TC(:,i)+(n_areas/2-i)*line_offset,'color',colormat_time(i+1,:));
            hold on;
        else
            axes(RH);
            plot(time_vec,TC(:,i)+(n_areas-i)*line_offset,'color',colormat_time(i-n_areas/2+1,:));
            hold on;
        end
    end
    % Barrel ROIs (same as behavior_analysis>0 path): required for trial_data.TC_all_trials_barrels / analyze_*_barrels
    for i = 1:size(mask_barrels, 3)
        TC_all_trials_barrels(:, i, :) = squeeze(nanmean(nanmean(df_movie .* repmat(squeeze(mask_barrels(:, :, i)), [1, 1, Nimages, Nstim]))));
    end
    save([analysis_dir 'Wide_Field\' mouse_name '\' side '\' run_date '\' run_name '\Results.mat'],'TC');

    
    axes(LH);
    ylim([min(TC(:)) n_areas/2*line_offset + max(TC(:))]);
    xline(delay);
    xline(delay+stim_duration);
    yticks([0:line_offset:((n_areas/2-1)*line_offset)]);
    yticklabels(flipud(cell_areanames(1:n_areas/2)));
    title('TC LH');
    
    axes(RH);
    ylim([min(TC(:)) n_areas/2*line_offset + max(TC(:))]);
    xline(delay);
    xline(delay+stim_duration);
    yticks([0:line_offset:((n_areas/2-1)*line_offset)]);
    yticklabels(flipud(cell_areanames((n_areas/2+1):end)));
    title(['TC RH, ' num2str(min(TC(:))) ' to ' num2str(max(TC(:)))]);
    
else 
           
    LH_with_wipes = subplot(2,2,1);
    RH_with_wipes = subplot(2,2,2);
    LH_wo_wipes = subplot(2,2,3);
    RH_wo_wipes = subplot(2,2,4);
    
    for i = 1:n_areas
        
        mask_temp = single(allen_transformed == area_num(i));
        mask_temp(mask_temp == 0) = NaN;
        
        TC_with_wipes(:,i) = squeeze(nanmean(nanmean(mean_df_movie_with_wipes.*repmat(mask_temp,[1,1,Nimages]))));
        TC_wo_wipes(:,i) = squeeze(nanmean(nanmean(mean_df_movie_wo_wipes.*repmat(mask_temp,[1,1,Nimages]))));
        TC_with_wipes_green(:,i) = squeeze(nanmean(nanmean(mean_df_green_with_wipes.*repmat(mask_temp,[1,1,Nimages]))));
        TC_wo_wipes_green(:,i) = squeeze(nanmean(nanmean(mean_df_green_wo_wipes.*repmat(mask_temp,[1,1,Nimages]))));
        TC_all_trials(:,i,:)= squeeze(nanmean(nanmean(df_movie.*repmat(mask_temp,[1,1,Nimages]))));


        
        if i <= n_areas/2
            axes(LH_with_wipes);
            plot(time_vec,TC_with_wipes(:,i)+(n_areas/2-i)*line_offset,'color',colormat_time(i+1,:));
            hold on;
            axes(LH_wo_wipes);
            plot(time_vec,TC_wo_wipes(:,i)+(n_areas/2-i)*line_offset,'color',colormat_time(i+1,:));
            hold on;
        else
            axes(RH_with_wipes);
            plot(time_vec,TC_with_wipes(:,i)+(n_areas-i)*line_offset,'color',colormat_time(i-n_areas/2+1,:));
            hold on;
            axes(RH_wo_wipes);
            plot(time_vec,TC_wo_wipes(:,i)+(n_areas-i)*line_offset,'color',colormat_time(i-n_areas/2+1,:));
            hold on;
        end
    end
    
    
    for i = 1:size(mask_barrels,3)
        TC_all_trials_barrels(:,i,:)=squeeze(nanmean(nanmean(df_movie.*repmat(squeeze(mask_barrels(:,:,i)),[1,1,Nimages,Nstim]))));

    end



%     TC_all_trials_barrels_with_wipes=TC_all_trials_barrels(:,:,find(trials_with_wipes));
%     TC_all_trials_barrels_wo_wipes=TC_all_trials_barrels(:,:,find(trials_wo_wipes));
    
    ylimits = [min(min(TC_with_wipes(:),TC_wo_wipes(:))) n_areas/2*line_offset + max(max(TC_with_wipes(:),TC_wo_wipes(:)))];
    
    axes(LH_with_wipes);
    ylim(ylimits);
    LabelFontSizeMultiplier = 0.7;
    xline(delay);
    xline(delay+stim_duration);
    yticks([0:line_offset:((n_areas/2-1)*line_offset)]);
    yticklabels(flipud(cell_areanames(1:n_areas/2)));
    title(['TC LH, Nwipes = ' num2str(N_with_wipes)]);
    
    axes(RH_with_wipes);
    ylim(ylimits);
    LabelFontSizeMultiplier = 0.7;
    xline(delay);
    xline(delay+stim_duration);
    yticks([0:line_offset:((n_areas/2-1)*line_offset)]);
    yticklabels(flipud(cell_areanames((n_areas/2+1):end)));
    title(['TC RH, range ' num2str(ylimits(1)) ' to ' num2str(ylimits(2))]);
    
    axes(LH_wo_wipes);
    ylim(ylimits);
    LabelFontSizeMultiplier = 0.7;
    xline(delay);
    xline(delay+stim_duration);
    yticks([0:line_offset:((n_areas/2-1)*line_offset)]);
    yticklabels(flipud(cell_areanames(1:n_areas/2)));
    title(['TC LH,Nwowipes = ' num2str(N_wo_wipes)]);
    
    axes(RH_wo_wipes);
    ylim(ylimits);
    LabelFontSizeMultiplier = 0.7;
    xline(delay);
    xline(delay+stim_duration);
    yticks([0:line_offset:((n_areas/2-1)*line_offset)]);
    yticklabels(flipud(cell_areanames((n_areas/2+1):end)));
    
end

% plot trajectories

if behavior_analysis>0
    fig_trajectories = figure('Position',[0 550 1500 280]);
    
    set(fig_trajectories,'PaperOrientation','landscape');
    set(fig_trajectories,'PaperType','uslegal');
    
    if stim_L == 0
        x = round(righthand_x'); % in DLC, the stimulus is always on the right
        y = round(righthand_y');
    else
        x = flip(round(righthand_x'));  % movies are flipped before analysis in DLC, need to flip back
        y = round(righthand_y');
        
        
    end
    z = zeros(size(x));

    mask_BC = single(allen_transformed == area_num(28 + stim_L*28)); % BC
    mask_BC(mask_BC == 0) = NaN;
    TC_BC = squeeze(nanmean(nanmean(df_movie.*repmat(mask_BC,[1,1,Nimages,Nstim]))));
    TC_BC = reshape(TC_BC,[1,Nimages*Nstim]);
    
    mask_DZ = single(allen_transformed == area_num(22 + stim_L*28)); % DZ
    mask_DZ(mask_DZ == 0) = NaN;
    TC_DZ = squeeze(nanmean(nanmean(df_movie.*repmat(mask_DZ,[1,1,Nimages,Nstim]))));
    TC_DZ = reshape(TC_DZ,[1,Nimages*Nstim]);
    
    mask_FL = single(allen_transformed == area_num(25 + stim_L*28)); % FL
    mask_FL(mask_FL == 0) = NaN;
    TC_FL = squeeze(nanmean(nanmean(df_movie.*repmat(mask_FL,[1,1,Nimages,Nstim]))));
    TC_FL = reshape(TC_FL,[1,Nimages*Nstim]);
    
    % Adaptive normalization limits for trajectory overlays per area
    min_df = min(TC_BC(:));
    max_df = max(TC_BC(:));
    range_df = max(max_df - min_df, eps);
    min_df_BC = min_df - 0.1 * range_df;
    max_df_BC = max_df + 0.1 * range_df;
    
    subplot(1,3,1);
    imagesc(mouse_portrait);
    hold on;
    col = (TC_BC - min_df_BC)/(max_df_BC - min_df_BC);
    surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.1);
    title('BC');
    
    subplot(1,3,2);
    imagesc(mouse_portrait);
    hold on;
    min_df = min(TC_DZ(:));
    max_df = max(TC_DZ(:));
    range_df = max(max_df - min_df, eps);
    min_df_DZ = min_df - 0.1 * range_df;
    max_df_DZ = max_df + 0.1 * range_df;
    col = (TC_DZ - min_df_DZ)/(max_df_DZ - min_df_DZ);
    surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.1);
    title('DZ');
    
    subplot(1,3,3);
    imagesc(mouse_portrait);
    hold on;
    min_df = min(TC_FL(:));
    max_df = max(TC_FL(:));
    range_df = max(max_df - min_df, eps);
    min_df_FL = min_df - 0.1 * range_df;
    max_df_FL = max_df + 0.1 * range_df;
    col = (TC_FL - min_df_FL)/(max_df_FL - min_df_FL);
    surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.1);
    title('FL');
    saveas(fig_trajectories,[analysis_dir 'Wide_Field\' mouse_name '\' side '\' run_date '\' run_name '\fig_trajectories' stim_type '_' stim_param '_' date '.pdf']);
    savefig(fig_trajectories,[analysis_dir 'Wide_Field\' mouse_name '\' side '\' run_date '\' run_name '\fig_trajectories' stim_type '_' stim_param '_' date '.fig']);
end


% Save (trial_data + Results before figure export so a plotting error cannot drop outputs)
if save_data

    out_dir = fullfile(analysis_dir, 'Wide_Field', mouse_name, side, run_date, run_name);
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    trial_mat = fullfile(out_dir, 'trial_data.mat');
    results_mat = fullfile(out_dir, 'Results.mat');

    save(trial_mat,...
        'mouse_name','run_date','run_name','stim_type','stim_param','Nstim','Nimages','Nx','Ny','trial_duration','delay','stim_duration','FrameRate_PCO','FrameRate_Basler',...
        'BehaviorFrames','trials_wo_wipes','N_wo_wipes','trials_with_wipes','N_with_wipes','trials_with_groom','N_with_groom',...
        'allen_transformed','allen_boundaries','mask_global','TC_all_trials','TC_all_trials_barrels',...
        'barrel_roi_atlas_idx','barrel_roi_names','barrel_roi_tc_desc',...
        'df_movie');  %to change
    fprintf(1, 'Saved %s\n', trial_mat);

    % Trial-by-trial df/f for S1bf_R only: (allen_transformed == area_num(i)), same as per-area TC loop.
    wf_run_dir = out_dir;
    mask_barrel_cortex_registered = false(Ny, Nx);
    TC_barrel_cortex_pixels = [];
    pixel_linear_idx = [];
    barrel_cortex_pixel_desc = '';
    mask_ok = false;
    barrel_cortex_pixel_atlas_index_used = [];
    barrel_cortex_pixel_atlas_area_num_value = [];
    barrel_cortex_atlas_area_names = {};
    mask_barrel_cortex_area_px_atlas_union = 0;
    mask_barrel_cortex_area_px_post_clip = 0;
    barrel_cortex_atlas_resized_to_df = false;
    barrel_cortex_pixel_roi_resolution = '';
    try
        Ta = allen_transformed;
        if ~isempty(Ta) && ~isequal(size(Ta), [Ny, Nx])
            Ta = imresize(double(Ta), [Ny, Nx], 'nearest');
            barrel_cortex_atlas_resized_to_df = true;
        end
        if isempty(Ta) || ~isequal(size(Ta), [Ny, Nx])
            barrel_cortex_pixel_desc = sprintf( ...
                'allen_transformed missing or size %s incompatible with df_movie [%d %d]', ...
                mat2str(size(allen_transformed)), Ny, Nx);
        else
            roi_nm = strtrim(char(string(barrel_cortex_pixel_roi_name)));
            ai = [];
            for k = 1:n_areas
                if strcmpi(strtrim(char(string(cell_areanames{k}))), roi_nm)
                    ai = k;
                    break
                end
            end
            if ~isempty(ai)
                barrel_cortex_pixel_roi_resolution = 'name_match';
            elseif isnumeric(barrel_cortex_pixel_roi_linear_index) && ...
                    isscalar(barrel_cortex_pixel_roi_linear_index) && ...
                    barrel_cortex_pixel_roi_linear_index >= 1 && ...
                    barrel_cortex_pixel_roi_linear_index <= n_areas
                ai = double(barrel_cortex_pixel_roi_linear_index);
                barrel_cortex_pixel_roi_resolution = sprintf('linear_index_%d_name_not_found', ai);
            end
            if isempty(ai)
                barrel_cortex_pixel_desc = sprintf( ...
                    'could not resolve ROI ''%s'' (linear fallback %g invalid for n_areas=%d)', ...
                    roi_nm, double(barrel_cortex_pixel_roi_linear_index), n_areas);
            else
                anm = char(string(cell_areanames{ai}));
                av = area_num(ai);
                mask_barrel_cortex_registered = (Ta == av);
                barrel_cortex_pixel_atlas_index_used = ai;
                barrel_cortex_pixel_atlas_area_num_value = av;
                barrel_cortex_atlas_area_names = {anm};
                mask_barrel_cortex_area_px_atlas_union = sum(mask_barrel_cortex_registered(:));
                if mask_barrel_cortex_area_px_atlas_union == 0
                    barrel_cortex_pixel_desc = sprintf( ...
                        'no pixels with allen_transformed==area_num(%d) (label %g, name ''%s'')', ...
                        ai, double(av), anm);
                else
                    mask_ok = true;
                    barrel_cortex_pixel_desc = sprintf( ...
                        'S1bf_R mask: cell_areanames{%d}=''%s'', area_num=%g, resolution=%s', ...
                        ai, anm, double(av), barrel_cortex_pixel_roi_resolution);
                end
            end
        end
    catch ME_bc
        barrel_cortex_pixel_desc = sprintf('barrel_mask_failed: %s', ME_bc.message);
    end
    if mask_ok && any(mask_barrel_cortex_registered(:))
        if ~isempty(mask_global) && isequal(size(mask_global), [Ny, Nx])
            cort = isfinite(mask_global) & (double(mask_global) ~= 0);
            mask_barrel_cortex_registered = mask_barrel_cortex_registered & cort;
        end
        if ~isempty(allen_boundaries) && isequal(size(allen_boundaries), [Ny, Nx])
            mask_barrel_cortex_registered = mask_barrel_cortex_registered & ~logical(allen_boundaries);
        end
        mask_barrel_cortex_area_px_post_clip = sum(mask_barrel_cortex_registered(:));
        if ~any(mask_barrel_cortex_registered(:))
            barrel_cortex_pixel_desc = [barrel_cortex_pixel_desc ' (empty after cortex/boundary clip)'];
        end
    end
    if mask_ok && any(mask_barrel_cortex_registered(:))
        pixel_linear_idx = find(mask_barrel_cortex_registered(:));
        if ~isempty(df_movie) && isequal(size(df_movie), [Ny, Nx, Nimages, Nstim])
            vm = reshape(df_movie, Ny * Nx, Nimages, Nstim);
            pix_mat = vm(pixel_linear_idx, :, :);
            TC_barrel_cortex_pixels = permute(pix_mat, [2, 1, 3]);
        else
            if isempty(df_movie)
                barrel_cortex_pixel_desc = [barrel_cortex_pixel_desc ' | df_movie_empty'];
            else
                barrel_cortex_pixel_desc = [barrel_cortex_pixel_desc sprintf(' | df_movie_size_%s_not_[Ny,Nx,Nimages,Nstim]', mat2str(size(df_movie)))];
            end
        end
        save(fullfile(wf_run_dir, 'barrel_cortex_pixel_timecourse.mat'), ...
            'mask_barrel_cortex_registered', 'pixel_linear_idx', 'TC_barrel_cortex_pixels', ...
            'time_vec', 'FrameRate_PCO', 'Nimages', 'Nx', 'Ny', 'Nstim', 'Frame_stimOn', 'delay', ...
            'barrel_cortex_pixel_desc', ...
            'barrel_cortex_pixel_roi_name', 'barrel_cortex_pixel_roi_linear_index', ...
            'barrel_cortex_pixel_atlas_index_used', 'barrel_cortex_pixel_atlas_area_num_value', ...
            'barrel_cortex_pixel_roi_resolution', 'barrel_cortex_atlas_area_names', ...
            'barrel_cortex_atlas_resized_to_df', ...
            'mask_barrel_cortex_area_px_atlas_union', 'mask_barrel_cortex_area_px_post_clip', ...
            'barrel_idx', 'barrel_roi_atlas_idx', 'barrel_roi_names', ...
            'mouse_name', 'run_date', 'run_name', 'stim_type', 'stim_param', '-v7.3');
        fprintf(1, 'Saved %s\n', fullfile(wf_run_dir, 'barrel_cortex_pixel_timecourse.mat'));
    end

    save(results_mat,...
        'mouse_name','run_date','run_name','stim_type','stim_param','Nstim','Nimages','Nx','Ny','trial_duration','delay','stim_duration','Frame_stimOn','FrameRate_PCO','FrameRate_Basler',...
        'BehaviorFrames','trials_wo_wipes','N_wo_wipes','trials_with_wipes','N_with_wipes','trials_with_groom','N_with_groom',...
        'allen_transformed','allen_boundaries','mask_global',...
        'mean_df_movie','mean_df_green','mean_df_movie_with_wipes','mean_df_movie_wo_wipes',...
        'mean_df_green_with_wipes','mean_df_green_wo_wipes',...
        'df_map_with_wipes','df_map_wo_wipes','TC_with_wipes','TC_wo_wipes','TC_with_wipes_green','TC_wo_wipes_green');
    fprintf(1, 'Saved %s\n', results_mat);

    try
        saveas(fig_TCs_all, fullfile(out_dir, ['fig_TCs_' stim_type '_' stim_param '_' date '.pdf']));
        saveas(fig_df_maps, fullfile(out_dir, ['DF_maps_' stim_type '_' stim_param '_' date '.pdf']));
    catch ME_fig
        warning('Analyze_one_WF_run:FigureExportFailed', '%s', ME_fig.message);
    end
    % DF maps: Illustrator-friendly SVG (vector axes/text + one embedded image per imagesc panel)
    try
        if isgraphics(fig_df_maps) && isvalid(fig_df_maps)
            df_svg = fullfile(out_dir, ['DF_maps_' stim_type '_' stim_param '_' date '.svg']);
            wf_debug_figure_vector_export(fig_df_maps, df_svg, 'svg', 1200, 'Mode', 'vector');
        end
    catch ME_svg
        warning('Analyze_one_WF_run:DFMapSVGExport', '%s', ME_svg.message);
    end

    if save_videos
        fprintf(1, 'Writing MP4 movies to: %s\n', out_dir);
        if behavior_analysis == 0
        
        % Adaptive normalization for mean_df_movie
        temp_min = min(mean_df_movie(:));
        temp_max = max(mean_df_movie(:));
        temp_range = max(temp_max - temp_min, eps);
        low = temp_min - 0.1 * temp_range;
        high = temp_max + 0.1 * temp_range;
        mean_df_movie_rescaled = (mean_df_movie - low)/(high - low);
        mean_df_movie_rescaled(mean_df_movie_rescaled>1) = 1;
        mean_df_movie_rescaled(mean_df_movie_rescaled<0) = 0;
        mean_movie_rescaled = ~allen_boundaries.*ceil(mean_df_movie_rescaled * (length(colormat_space) - 1));
        mean_movie_rescaled(isnan(mean_movie_rescaled)) = 0;
           
        mean_movie_rescaled_RGB = ind2rgb(reshape(mean_movie_rescaled,[Ny*Nx*Nimages,1]),colormat_space);
        mean_movie_rescaled_R = reshape(mean_movie_rescaled_RGB(:,1),[Ny,Nx,Nimages]);
        mean_movie_rescaled_G = reshape(mean_movie_rescaled_RGB(:,2),[Ny,Nx,Nimages]);
        mean_movie_rescaled_B = reshape(mean_movie_rescaled_RGB(:,3),[Ny,Nx,Nimages]);
        mean_movie_rescaled_RGB = cat(3,permute(mean_movie_rescaled_R,[1,2,4,3]),permute(mean_movie_rescaled_G,[1,2,4,3]),permute(mean_movie_rescaled_B,[1,2,4,3]));
        
        v = VideoWriter(fullfile(out_dir, 'mean_df_movie'), 'MPEG-4');
        v.FrameRate = round(FrameRate_PCO/2);
        try
            open(v);
            writeVideo(v, mean_movie_rescaled_RGB);
            close(v);
            fprintf(1, 'Wrote %s.mp4\n', fullfile(out_dir, 'mean_df_movie'));
        catch ME_vid
            try, close(v); catch, end
            warning('Analyze_one_WF_run:VideoWriteFailed', '%s', ME_vid.message);
        end
           
        else
              
        % Adaptive normalization for mean_df_movie_with_wipes
        temp_min = min(mean_df_movie_with_wipes(:));
        temp_max = max(mean_df_movie_with_wipes(:));
        temp_range = max(temp_max - temp_min, eps);
        low = temp_min - 0.1 * temp_range;
        high = temp_max + 0.1 * temp_range;
        mean_df_movie_rescaled = (mean_df_movie_with_wipes - low)/(high - low);
        mean_df_movie_rescaled(mean_df_movie_rescaled>1) = 1;
        mean_df_movie_rescaled(mean_df_movie_rescaled<0) = 0;       
        mean_movie_rescaled = ~allen_boundaries.*ceil(mean_df_movie_rescaled * (length(colormat_space) - 1));
        mean_movie_rescaled(isnan(mean_movie_rescaled)) = 0;
           
        mean_movie_rescaled_RGB = ind2rgb(reshape(mean_movie_rescaled,[Ny*Nx*Nimages,1]),colormat_space);
        mean_movie_rescaled_R = reshape(mean_movie_rescaled_RGB(:,1),[Ny,Nx,Nimages]);
        mean_movie_rescaled_G = reshape(mean_movie_rescaled_RGB(:,2),[Ny,Nx,Nimages]);
        mean_movie_rescaled_B = reshape(mean_movie_rescaled_RGB(:,3),[Ny,Nx,Nimages]);
        mean_movie_rescaled_RGB = cat(3,permute(mean_movie_rescaled_R,[1,2,4,3]),permute(mean_movie_rescaled_G,[1,2,4,3]),permute(mean_movie_rescaled_B,[1,2,4,3]));
        
        v = VideoWriter(fullfile(out_dir, 'mean_df_movie_with_wipes'), 'MPEG-4');
        v.FrameRate = round(FrameRate_PCO/2);
        try
            open(v);
            writeVideo(v, mean_movie_rescaled_RGB);
            close(v);
            fprintf(1, 'Wrote %s.mp4\n', fullfile(out_dir, 'mean_df_movie_with_wipes'));
        catch ME_vid
            try, close(v); catch, end
            warning('Analyze_one_WF_run:VideoWriteFailed', '%s', ME_vid.message);
        end
        
        % Adaptive normalization for mean_df_movie_wo_wipes
        temp_min = min(mean_df_movie_wo_wipes(:));
        temp_max = max(mean_df_movie_wo_wipes(:));
        temp_range = max(temp_max - temp_min, eps);
        low = temp_min - 0.1 * temp_range;
        high = temp_max + 0.1 * temp_range;
        mean_df_movie_rescaled = (mean_df_movie_wo_wipes - low)/(high - low);
        mean_df_movie_rescaled(mean_df_movie_rescaled>1) = 1;
        mean_df_movie_rescaled(mean_df_movie_rescaled<0) = 0;       
        mean_movie_rescaled = ~allen_boundaries.*ceil(mean_df_movie_rescaled * (length(colormat_space) - 1));
        mean_movie_rescaled(isnan(mean_movie_rescaled)) = 0;
          
        mean_movie_rescaled_RGB = ind2rgb(reshape(mean_movie_rescaled,[Ny*Nx*Nimages,1]),colormat_space);
        mean_movie_rescaled_R = reshape(mean_movie_rescaled_RGB(:,1),[Ny,Nx,Nimages]);
        mean_movie_rescaled_G = reshape(mean_movie_rescaled_RGB(:,2),[Ny,Nx,Nimages]);
        mean_movie_rescaled_B = reshape(mean_movie_rescaled_RGB(:,3),[Ny,Nx,Nimages]);
        mean_movie_rescaled_RGB = cat(3,permute(mean_movie_rescaled_R,[1,2,4,3]),permute(mean_movie_rescaled_G,[1,2,4,3]),permute(mean_movie_rescaled_B,[1,2,4,3]));
                
        v = VideoWriter(fullfile(out_dir, 'mean_df_movie_wo_wipes'), 'MPEG-4');
        v.FrameRate = round(FrameRate_PCO/2);
        try
            open(v);
            writeVideo(v, mean_movie_rescaled_RGB);
            close(v);
            fprintf(1, 'Wrote %s.mp4\n', fullfile(out_dir, 'mean_df_movie_wo_wipes'));
        catch ME_vid
            try, close(v); catch, end
            warning('Analyze_one_WF_run:VideoWriteFailed', '%s', ME_vid.message);
        end
        
        % Export full df/f movies for trials with wipes and without wipes
        if any(trials_with_wipes)
            df_with = df_movie(:,:,:,logical(trials_with_wipes));
            temp_min = min(df_with(:));
            temp_max = max(df_with(:));
            temp_range = max(temp_max - temp_min, eps);
            low = temp_min - 0.1 * temp_range;
            high = temp_max + 0.1 * temp_range;
            df_with_norm = (df_with - low) / (high - low);
            df_with_norm(df_with_norm > 1) = 1;
            df_with_norm(df_with_norm < 0) = 0;
            df_with_norm(isnan(df_with_norm)) = 0;
            
            df_with_uint8 = uint8(df_with_norm * (length(colormat_space) - 1));
            df_with_RGB = ind2rgb(reshape(df_with_uint8,[Ny*Nx*Nimages*size(df_with,4),1]),colormat_space);
            df_with_R = reshape(df_with_RGB(:,1),[Ny,Nx,Nimages*size(df_with,4)]);
            df_with_G = reshape(df_with_RGB(:,2),[Ny,Nx,Nimages*size(df_with,4)]);
            df_with_B = reshape(df_with_RGB(:,3),[Ny,Nx,Nimages*size(df_with,4)]);
            df_with_RGB = cat(3,permute(df_with_R,[1,2,4,3]),permute(df_with_G,[1,2,4,3]),permute(df_with_B,[1,2,4,3]));
            
            v = VideoWriter(fullfile(out_dir, 'df_movie_with_wipes'), 'MPEG-4');
            v.FrameRate = round(FrameRate_PCO/2);
            try
                open(v);
                writeVideo(v, df_with_RGB);
                close(v);
                fprintf(1, 'Wrote %s.mp4\n', fullfile(out_dir, 'df_movie_with_wipes'));
            catch ME_vid
                try, close(v); catch, end
                warning('Analyze_one_WF_run:VideoWriteFailed', '%s', ME_vid.message);
            end
        end
        
        if any(trials_wo_wipes)
            df_without = df_movie(:,:,:,logical(trials_wo_wipes));
            temp_min = min(df_without(:));
            temp_max = max(df_without(:));
            temp_range = max(temp_max - temp_min, eps);
            low = temp_min - 0.1 * temp_range;
            high = temp_max + 0.1 * temp_range;
            df_without_norm = (df_without - low) / (high - low);
            df_without_norm(df_without_norm > 1) = 1;
            df_without_norm(df_without_norm < 0) = 0;
            df_without_norm(isnan(df_without_norm)) = 0;
            
            df_without_uint8 = uint8(df_without_norm * (length(colormat_space) - 1));
            df_without_RGB = ind2rgb(reshape(df_without_uint8,[Ny*Nx*Nimages*size(df_without,4),1]),colormat_space);
            df_without_R = reshape(df_without_RGB(:,1),[Ny,Nx,Nimages*size(df_without,4)]);
            df_without_G = reshape(df_without_RGB(:,2),[Ny,Nx,Nimages*size(df_without,4)]);
            df_without_B = reshape(df_without_RGB(:,3),[Ny,Nx,Nimages*size(df_without,4)]);
            df_without_RGB = cat(3,permute(df_without_R,[1,2,4,3]),permute(df_without_G,[1,2,4,3]),permute(df_without_B,[1,2,4,3]));
            
            v = VideoWriter(fullfile(out_dir, 'df_movie_without_wipes'), 'MPEG-4');
            v.FrameRate = round(FrameRate_PCO/2);
            try
                open(v);
                writeVideo(v, df_without_RGB);
                close(v);
                fprintf(1, 'Wrote %s.mp4\n', fullfile(out_dir, 'df_movie_without_wipes'));
            catch ME_vid
                try, close(v); catch, end
                warning('Analyze_one_WF_run:VideoWriteFailed', '%s', ME_vid.message);
            end
        end
        
        
%         movie_temp = ~allen_boundaries .* (df_movie - min(df_movie(:))) / (max(df_movie(:))-min(df_movie(:)));
        % Adaptive normalization for df_movie overlays
        temp_min = min(df_movie(:));
        temp_max = max(df_movie(:));
        temp_range = max(temp_max - temp_min, eps);
        low = temp_min - 0.1 * temp_range;
        high = temp_max + 0.1 * temp_range;
        movie_temp = ~allen_boundaries .* (df_movie - low) / (high - low);
       
        movie_temp(movie_temp>1) = 1;
        movie_temp(movie_temp<0) = 0;
        movie_temp(isnan(movie_temp)) = 0;
                 
        movie_temp = ceil(movie_temp * (length(colormat_space) - 1));
        movie_temp = uint8(movie_temp);
       
        movie_temp_RGB = ind2rgb(reshape(movie_temp,[Ny*Nx*Nimages*Nstim,1]),colormat_space);
        
        movie_temp_R = reshape(movie_temp_RGB(:,1),[Ny,Nx,Nimages*Nstim]);
        movie_temp_G = reshape(movie_temp_RGB(:,2),[Ny,Nx,Nimages*Nstim]);
        movie_temp_B = reshape(movie_temp_RGB(:,3),[Ny,Nx,Nimages*Nstim]);
                
        movie_temp_RGB = cat(3,permute(movie_temp_R,[1,2,4,3]),permute(movie_temp_G,[1,2,4,3]),permute(movie_temp_B,[1,2,4,3]));
        
        movie_Be = (movie_Be - min(movie_Be(:)))/(max(movie_Be(:)) - min(movie_Be(:)));
        
        movie_Be(20:40,20:40,1,Frame_stimOn:Frame_stimOff,:) = 1;
        movie_Be(20:40,20:40,2,Frame_stimOn:Frame_stimOff,:) = 0;
        movie_Be(20:40,20:40,3,Frame_stimOn:Frame_stimOff,:) = 0;
        
        movie_Be = reshape(movie_Be,[Ny,Nx,3,Nimages*Nstim]);
                
        movie_temp = [movie_temp_RGB,movie_Be];
        v = VideoWriter(fullfile(out_dir, 'Movie_WF_and_Be'), 'MPEG-4');
        v.FrameRate = round(FrameRate_PCO/2);
        try
            open(v);
            writeVideo(v, movie_temp);
            close(v);
            fprintf(1, 'Wrote %s.mp4\n', fullfile(out_dir, 'Movie_WF_and_Be'));
        catch ME_vid
            try, close(v); catch, end
            warning('Analyze_one_WF_run:VideoWriteFailed', '%s', ME_vid.message);
        end
       
        end        
    end
end

if close_figs
    close all;
end

toc

end

function [vmin, vmax] = df_map_clim(vals)
%DF_MAP_CLIM  Robust CLim for a row of df maps.
%   Uses the 1st–99th percentile of finite values (outlier-tolerant) so a
%   handful of edge / regression-artefact pixels cannot crush the colormap.
%   Symmetrizes around zero when the data straddle zero, which matches how
%   we interpret df/f (positive & negative responses comparable).
    vals = double(vals(:));
    vals = vals(isfinite(vals));
    if isempty(vals)
        vmin = -eps;
        vmax = eps;
        return;
    end
    lo = prctile(vals, 1);
    hi = prctile(vals, 99);
    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
        lo = min(vals);
        hi = max(vals);
    end
    if lo < 0 && hi > 0
        a = max(abs(lo), abs(hi));
        vmin = -a;
        vmax =  a;
    else
        vmin = lo;
        vmax = hi;
    end
    if vmin == vmax
        vmin = vmin - eps;
        vmax = vmax + eps;
    end
end

function overlay_atlas_boundaries_solid(boundary_mask, rgb)
%OVERLAY_ATLAS_BOUNDARIES_SOLID  Draw atlas boundaries as a solid color, constant across scales.
%   boundary_mask: logical 2-D mask where true marks boundary pixels.
%   rgb: 1x3 color in [0,1], e.g. [0 0 0] black, [1 1 1] white.
    if nargin < 2 || isempty(rgb)
        rgb = [0 0 0];
    end
    if ~any(boundary_mask(:))
        return;
    end
    [H, W] = size(boundary_mask);
    overlay_img = repmat(reshape(rgb, [1 1 3]), [H, W, 1]);
    h = image(overlay_img);
    set(h, 'AlphaData', double(boundary_mask));
end
