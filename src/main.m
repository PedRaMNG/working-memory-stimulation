cd(fileparts(mfilename('fullpath')));
try
    tic;
    %% Initialization
    % close all; 
    clearvars; clc;     
    rng(42);
    fprintf('Initializing...');
    [params, paramss] = model_parameters(1);
    paramss_fields = fieldnames(paramss);
    [params, paramss] = show_parameter_gui(params, paramss, paramss_fields);
    params = model_parameters(1, params, paramss);
    clear paramss;
	model  = init_model();
    
    %% Simulation
    fprintf('\nInitialization completed.\n--- Starting Simulation ---\n');
    disp(['Iteration: 1/', num2str(params.n), ...
      ' (', num2str((1 / (params.n - 1)) * 100, '%.1f'), '% complete)']);
    
    % sendTelegramMessage('Simulation has started');
    [model] = simulate_model(model, params);

    % Final timing and notification
    totalTime = toc;
    [hours, minutes, seconds] = convertTime(totalTime);
    timeStr = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
    
    % Send notifications
    % sendTelegramMessage(sprintf('Simulation completed!\nTotal runtime: %s', timeStr));
    
%% Compute memory performance
    if params.simPattern ~= 3
        [memory_performance] = ...
            compute_memory_performance(model.images, model.V_line_E, model.T_Iapp);
        
        fprintf('Mean memory performance: %0.4f\n', memory_performance.mean_performance);
        fmt = repmat(' %0.4f', 1, numel(memory_performance.learned_pattern_similarities));
        fprintf(['Memory performance per image: ', fmt, '\n'], ...
            memory_performance.learned_pattern_similarities);
        % For experiment 1 (Enhance) and 3 (augmentation), you have to ...
        % ... use compute_memory_performance2 wtih some manual changes instructed in that file.
    end

%% Plot or Save data
	% run SaveData.m
    % run main_video_fast.m
    run main_video_rgb.m
    % run plot_settings.m

    fprintf('\nTotal simulation time: %s\n', timeStr);

catch ME
    if (strcmp(ME.identifier,'MATLAB:nomem'))
        error(['Out of memory. ' ...
            'Please, increase the amount of available memory or ' ...
            'change the simulation settings to not save unnecessary data. ' ...
            '\nThe required amount of RAM is 64 GB for experiment 1 ' ...
            'and 32 GB for experiments 2 and 3. ' ...
            '\n\nMake sure to clear data in the workspace before resimulation.'], 0);
    else
        rethrow(ME);
    end
end


function [hours, minutes, seconds] = convertTime(timeInSeconds)
    hours = floor(timeInSeconds / 3600);
    minutes = floor((timeInSeconds - hours * 3600) / 60);
    seconds = floor(timeInSeconds - hours * 3600 - minutes * 60);
end

function adjusted_matrix = adjust_initial_values(input_matrix)
    % Determine the dominant initial value
    initial_values = input_matrix(:, 1);
    unique_values = unique(initial_values);
    value_counts = histc(initial_values, unique_values);
    [max_count, idx] = max(value_counts);
    dominant_value = unique_values(idx);
    
    % Create an adjusted matrix
    adjusted_matrix = input_matrix;
    
    % Adjust each row to have the dominant initial value
    for i = 1:size(input_matrix, 1)
        adjustment = dominant_value - input_matrix(i, 1)-min(input_matrix(:, 1));
        % adjustment = dominant_value - input_matrix(i, 1);
        adjusted_matrix(i, :) = input_matrix(i, :) + adjustment;
    end
end
