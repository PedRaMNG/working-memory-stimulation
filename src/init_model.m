function model = init_model()

    model = struct;
    params = model_parameters();
    model.T = 0:params.step:params.t_end;
    model.T = single(model.T);

    %% Zone
    model.Pre_EE  = zeros(1, params.quantity_connections_EE, 'int8');
    model.Post_EE = zeros(1, params.quantity_connections_EE, 'int8');
    model.Pre_InE  = zeros(1, params.quantity_connections_InE, 'int8');
    model.Post_InE = zeros(1, params.quantity_connections_InE, 'int8');
    model.Pre_ExI  = zeros(1, params.quantity_connections_ExI, 'int8');
    model.Post_ExI = zeros(1, params.quantity_connections_ExI, 'int8');

    %% Neurons
    model.V_line_E = zeros(params.quantity_neurons_E, params.n, 'double');
    model.V_line_I = zeros(params.quantity_neurons_I, params.n, 'double');
    model.V_line_E(:, 1) = params.V_0;
    model.V_line_I(:, 1) = params.V_0;
    model.G = zeros(params.quantity_neurons_E, params.n, 'double'); %can be reduced
    model.U_line_E = zeros(params.quantity_neurons_E, 1, 'double');
    model.U_line_I = zeros(params.quantity_neurons_I, 1, 'double');
    model.U_line_E(:, 1) = params.U_0;
    model.U_line_I(:, 1) = params.U_0;
    model.Isyn_line_EE = zeros(params.quantity_neurons_E, 1,'double');
    model.Isyn_line_InE = zeros(params.quantity_neurons_E, 1,'double');
    model.Isyn_line_ExI = zeros(params.quantity_neurons_I, 1,'double');

    %% Neuron activity
    model.neuron_astrozone_activity_Glu = zeros(params.mastro, params.nastro, params.n, 'double'); %can be reduced
    model.neuron_astrozone_spikes = zeros(params.mastro, params.nastro, params.n, 'int8');

    %% Astrocytes
    model.Ineuro = zeros(params.mastro, params.nastro, params.n, 'int8');
    model.W_WM = zeros(params.mastro, params.nastro, params.n, 'logical'); 
    model.neuron_astrozone_input_timings = ones(params.mastro, params.nastro, params.n, 'double');
    model.I_app_astro = zeros(params.mastro, params.nastro, 1, 'logical');

    model.Ca = zeros(params.mastro, params.nastro, params.n, 'double'); %can be reduced
    model.Ca_expand = zeros(params.mneuro_E, params.nneuro_E, 1, 'double'); %can be reduced
    model.H = zeros(params.mastro, params.nastro, params.n, 'double');
    model.IP3 = zeros(params.mastro, params.nastro, params.n, 'double');
    model.Ca(:, :, 1) = params.ca_0;
    model.H(:, :, 1) = params.h_0;
    model.IP3(:, :, 1) = params.ip3_0;
    model.Gli_global = zeros(params.mastro, params.nastro, params.n, 'double');
    model.Gli_global_expand = zeros(params.mneuro_E, params.nneuro_E, 1, 'double');
    
    %% Iapp for video 
    model.Iapp_v_full = zeros(params.mneuro_E, params.nneuro_E, params.n, 'uint8');

    %% Prepare model
    [model.images, model.Damage_E, model.Damage_I, model.Damage_astro] = load_images(params.simPattern, params.mneuro_I);
    model.I_poisson_noise = make_poisson_noise();
    [model.Iapp, model.T_Iapp, model.T_Iapp_met, model.T_record_met] = make_experiment(model.images);
    
    %% creating connections
    [model.Pre_EE, model.Post_EE] = create_connections(1);
    [model.Pre_InE, model.Post_InE] = create_connections(2);
    [model.Pre_ExI, model.Post_ExI] = create_connections(3);
    
    [Post2_EE, indi] = sort(model.Post_EE);
    % model.Pre2_EE = model.Pre_EE(indi); % for testing, can be removed

    model.numofPost_EE = zeros(params.quantity_neurons_E, 1, 'double');
    for i = 1:params.quantity_connections_EE
        model.numofPost_EE(Post2_EE(i)) = model.numofPost_EE(Post2_EE(i)) + 1;
    end

    [Post2_ExI, ~] = sort(model.Post_ExI);
    numofPost_ExI = zeros(params.quantity_neurons_I, 1, 'double');
    for i = 1:params.quantity_connections_ExI
        numofPost_ExI(Post2_ExI(i)) = numofPost_ExI(Post2_ExI(i)) + 1;
    end

    [Post2_InE, ~] = sort(model.Post_InE);
    numofPost_InE = zeros(params.quantity_neurons_E, 1, 'double');
    for i = 1:params.quantity_connections_InE
        numofPost_InE(Post2_InE(i)) = numofPost_InE(Post2_InE(i)) + 1;
    end
    
    %% Network Damage mode
    [model.ST_S_EE, model.ST_N_EE] ...
        = impair(params.quantity_connections_EE, params, model.Damage_E, model.numofPost_EE, model.Post_EE);
    [model.ST_S_ExI] ...
        = impair(params.quantity_connections_ExI, params, model.Damage_I, numofPost_ExI, model.Post_ExI);
    [model.ST_S_InE] ...
        = impair(params.quantity_connections_InE, params, model.Damage_E, numofPost_InE, model.Post_InE);

    %% Stimulator
    
    [model.electrode_sensitivity_neuron_I_rec, ...
     model.electrode_sensitivity_neuron_E_rec] = electrode(params, 0); %rec

    [model.electrode_sensitivity_neuron_I_stim, ...
     model.electrode_sensitivity_neuron_E_stim] = electrode(params, 1); %stim
    % quick plot code
    % % model.tumVec = 1 - double(model.Damage_astro)/255; %temp
    % % imshow(model.tumVec);
    % % model.tumVec(model.tumVec > 0.004) = 1; %temp
    % % figure(); imshow(model.tumVec);

    if params.stimulation_mode == 0
        % Manual
        if params.stimulation_Augmentation
            
            model.manual_electrode_neuro_I_stim = sum(model.electrode_sensitivity_neuron_I_stim(:,params.electrode_neuron_stimulate), 2);
            model.manual_electrode_neuro_I_stim_timed = zeros(params.mneuro_I*params.nneuro_I, params.n, 'double');

            model.manual_electrode_neuro_E_stim_timed = zeros(params.mneuro_E*params.nneuro_E, params.n, 'double');

            for stim_events = 1:numel(params.stimulation_time_neuro_E_Augmentation)
                for i = params.stim_start_neuro_E_G(stim_events):params.stim_end_neuro_E_G(stim_events)
                    
                    model.manual_electrode_neuro_E_stim = ...
                    sum(model.electrode_sensitivity_neuron_E_stim(:, params.electrode_neuron_stimulate_G{stim_events}), 2);

                    model.manual_electrode_neuro_E_stim_timed(:,i) = ...
                    model.manual_electrode_neuro_E_stim .* params.stim_neuro_amp_E;

                end
            end
            % quick plot video
            % manual_electrode_neuro_E_stim_timed_sqr = reshape(model.manual_electrode_neuro_E_stim_timed, ...
            %                             params.mneuro_E,params.nneuro_E,params.n);
            % show_video(manual_electrode_neuro_E_stim_timed_sqr);
        else
            model.manual_electrode_neuro_E_stim = sum(model.electrode_sensitivity_neuron_E_stim(:,params.electrode_neuron_stimulate), 2);
            % reshaped_electrode_E = reshape(model.manual_electrode_neuro_E_stim, [params.mneuro_E params.nneuro_E]);
            % figure(); surf(reshaped_electrode_E'); view(90,90); title('Stimulating Electrodes Excitatory');
            model.manual_electrode_neuro_I_stim = sum(model.electrode_sensitivity_neuron_I_stim(:,params.electrode_neuron_stimulate), 2);
            % reshaped_electrode_I = reshape(model.manual_electrode_neuro_I_stim, [params.mneuro_I params.nneuro_I]);
            % figure(); surf(reshaped_electrode_I'); view(90,90); title('Stimulating Electrodes Inhibitory');
            % imwrite(model.manual_electrode_neuro_E_stim, 'SavedData/matrix_image.jpg');
            model.manual_electrode_neuro_E_stim_timed = zeros(params.mneuro_E*params.nneuro_E, params.n, 'double');
            for i = params.stim_start_neuro_E:params.stim_end_neuro_E
                model.manual_electrode_neuro_E_stim_timed(:,i) = ...
                model.manual_electrode_neuro_E_stim .* params.stim_neuro_amp_E;
            end

            model.manual_electrode_neuro_I_stim_timed = zeros(params.mneuro_I*params.nneuro_I, params.n, 'double');
            for stim_events = 1:numel(params.stim_start_neuro_I)
                for i = params.stim_start_neuro_I(stim_events):params.stim_end_neuro_I(stim_events)
                    model.manual_electrode_neuro_I_stim_timed(:,i) = ...
                    model.manual_electrode_neuro_I_stim .* params.stim_neuro_amp_I;
                end
            end
        end
    else
        % Auto
        model.electrode_neuron_E_rec = zeros(params.n_electrode_record, params.n);
        electrode_neuron_E_rec = zeros(params.n_electrode_record,1);
        model.electrode_neuron_I_rec = zeros(params.n_electrode_record, params.n);
        electrode_neuron_I_rec = zeros(params.n_electrode_record, 1);
        model.electrode_neuron_total_rec = zeros(params.n_electrode_record, params.n);
        electrode_neuron_total_rec = zeros(params.n_electrode_record, 1);
        model.auto_electrode_neuro_E_stim_timed = zeros(params.mneuro_E*params.nneuro_E, params.n, 'double');
        model.auto_electrode_neuro_I_stim_timed = zeros(params.mneuro_I*params.nneuro_I, params.n, 'double');
        model.auto_electrode_astro_stim_timed = zeros(params.mastro, params.nastro, params.n, 'double');

        model.E_rec_passed = zeros(params.electrode_grid_record(1), params.electrode_grid_record(2), params.n, 'logical'); 
        model.Total_E_rec_passed = zeros(params.electrode_grid_stim(1),params.electrode_grid_stim(2), 'int8');
        for n = 1:params.n_electrode_record
            electrode_neuron_E_rec(n,1) = sum(model.electrode_sensitivity_neuron_E_rec(:,n)*params.V_0);
            electrode_neuron_I_rec(n,1) = sum(model.electrode_sensitivity_neuron_I_rec(:,n)*params.V_0);
            electrode_neuron_total_rec(n,1) = electrode_neuron_E_rec(n,1) + electrode_neuron_I_rec(n,1);
        end
        model.adjusted_matrix = adjust_initial_values(electrode_neuron_total_rec);
    end 

    if params.stimulation_mode && params.stimulation_neuro_I 
        % Finding image target for disruption
        model.image_target_dist = zeros(params.n_electrode_record, numel(params.disruption_pattern));
        image_number = 1;
        for imG = params.disruption_pattern
            inputs = double(cell2mat(model.images(imG))) / 255;
            inputs = reshape(inputs, params.mneuro_E*params.nneuro_E,1);
            inputs = 1-inputs;
            for e = 1:params.n_electrode_record 
                electrode_activity = sum(inputs .* ...
                                        model.electrode_sensitivity_neuron_E_rec(:,e));  % Calculate activity
                
                if electrode_activity > 20
                    model.image_target_dist(e,image_number) = 1;  % Mark the electrode as active (1)
                else
                    model.image_target_dist(e,image_number) = 0;
                end
            
            end
            image_number = image_number + 1;
        end
        model.image_target_dist_sqr = reshape(model.image_target_dist,...
                                             [params.electrode_grid_record(1), ...
                                             params.electrode_grid_record(2), ...
                                             numel(params.disruption_pattern)]);

        model.flag_dist = zeros(1,params.n);
        model.dist_timeSteps = [];
        
        model.integral_electrode_rec_dist = zeros(params.n_electrode_record, params.n);
        model.integral_electrode_rec3 = zeros(params.n_electrode_record, params.n);
        
        model.integral_electrode_rec_dist = zeros(params.n_electrode_record, params.n/params.dis_sampling);
        model.integral_electrode_rec3 = zeros(params.n_electrode_record, params.n/params.dis_sampling);
    end

    % Electrode record calibration
    function adjusted_matrix = adjust_initial_values(input_matrix)
        % Determine the dominant initial value
        initial_values = input_matrix(:, 1);
        unique_values = unique(initial_values);
        value_counts = histc(initial_values, unique_values);
        [max_count, idx] = max(value_counts);
        dominant_value = unique_values(idx);
        
        % Create an adjusted matrix
        adjusted_matrix = zeros(size(input_matrix,1),1);
        
        % Adjust each row to have the dominant initial value
        for j = 1:size(input_matrix, 1)
            adjusted_matrix(j, 1) = dominant_value - input_matrix(j, 1);
        end
        % offseting everything to zero
        adjusted_matrix(:,1) = adjusted_matrix(:,1) ...
                             - input_matrix(1, 1) ...
                             - adjusted_matrix(1, 1);
    end

end



