%% Video with a single color map for all variables

% Record
optVideo1.electrode_neuron_E_rec = zeros(params.n_electrode_record, params.n);
optVideo1.electrode_neuron_I_rec = zeros(params.n_electrode_record, params.n);
optVideo1.electrode_neuron_total_rec = zeros(params.n_electrode_record, params.n);

for n = 1 : params.n_electrode_record %opt.electrode_record
    for i = 1:params.n
        optVideo1.electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
        optVideo1.electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
        optVideo1.electrode_neuron_total_rec(n,i) = optVideo1.electrode_neuron_E_rec(n,i) + optVideo1.electrode_neuron_I_rec(n,i);
    end
    disp(['Electrode: ', num2str(n), ' / ', num2str(params.n_electrode_record)]);
end
%normalizing the values
optVideo1.electrode_neuron_total_rec = adjust_initial_values(optVideo1.electrode_neuron_total_rec);

% Stim
if params.stimulation_mode == 0 % Manual
    if params.stimulation_neuro_I == 1
        optVideo1.electrode_neuro_stim_timed = model.manual_electrode_neuro_I_stim_timed;
    else
        optVideo1.electrode_neuro_stim_timed = model.manual_electrode_neuro_E_stim_timed;
    end
else % Auto
    if params.stimulation_neuro_I == 1
        optVideo1.electrode_neuro_stim_timed = model.auto_electrode_neuro_I_stim_timed;
    else
        optVideo1.electrode_neuro_stim_timed = model.auto_electrode_neuro_E_stim_timed;
    end
end

[model.video1] = make_video(model.Ca_expand, ...
                            model.V_line_E, ...
                            model.V_line_I, ...
                            model.Iapp_v_full, ...
                            model.T_record_met(1:params.n), ...
                            optVideo1.electrode_neuron_total_rec, ...
                            optVideo1.electrode_neuro_stim_timed, ...
                            params);
optVideo1.colrs = [255, 255, 255;
            110, 109, 98;
            29, 37, 8; 
            0.6*255,0.1*255,0.2*255]/255; %164, 57, 70     
optVideo1.locs = [0,0.33,0.66,1];
optVideo1.limits = [0, 255]; % 350

optVideo1.cmap = makeGradient(256, optVideo1.colrs, optVideo1.locs, 1);
optVideo1.fps = 100;
optVideo1.full_screen = false;
optVideo1.figure_position = [100, 100, 1050, 550];
show_video(model.video1, optVideo1);

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