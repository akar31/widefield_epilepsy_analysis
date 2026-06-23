function mask_global = applyCortexMask(analysis_dir, mouse_name, side, run_date, run_name, movie_blue)
% Masks.mat: same resolution as Transform.mat (nested stim folder, else session run1\...)
[mask_file_path, ~] = wf_widefield_run_paths(analysis_dir, mouse_name, side, run_date, run_name, 'brain_mask.mat');

if ~isempty(mask_file_path)
    load(mask_file_path, 'mask_global');
else
%     % Create a new mask if the file doesn't exist
%     [mask_global] = Make_cortex_mask(squeeze(movie_blue(:,:,1,1)));
%     mask_global(mask_global == 0) = NaN;
%     
%     % Save the newly created mask
%     save(mask_file_path, 'mask_global');
        % Open log file in append mode
        logFileID2 = fopen('missing_masks_lists.txt', 'a');  % 'a' for appending data to the file
        if logFileID2 == -1
            error('Failed to open log file.');
        end
        
        % Format the error message
        errorMsg2 = sprintf('no mask found dor: %s %s ', mouse_name,run_date);
        
        % Print error message to command window
        fprintf(errorMsg2);
        
        % Write error message to log file and close it
        fprintf(logFileID2, errorMsg2);
        fclose(logFileID2);
        
        % Continue with the next iteration of the loop
      
end

end
