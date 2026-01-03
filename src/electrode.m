function [electrode_sensitivity_neuron_I, ...
          electrode_sensitivity_neuron_E] = electrode(params, is_stim_or_rec)

    electrode_coordinates_E = ...
    calculate_electrode_coordinates(params, params.mneuro_E, params.nneuro_E, is_stim_or_rec);

    [electrode_sensitivity_neuron_E, ...
     electrode_sensitivity_neuron_I] = ...
        calculate_electrode_sensitivity_neuron(params, electrode_coordinates_E, is_stim_or_rec);

    electrode_sensitivity_neuron_E = reshape(electrode_sensitivity_neuron_E, [], size(electrode_sensitivity_neuron_E, 3));
    electrode_sensitivity_neuron_I = reshape(electrode_sensitivity_neuron_I, [], size(electrode_sensitivity_neuron_I, 3));
end

function electrode_positions = calculate_electrode_coordinates(params, mneuro, nneuro, is_stim_or_rec)
    if is_stim_or_rec == 1 %stimulation
        electrode_grid = params.electrode_grid_stim;
    else % record
        electrode_grid = params.electrode_grid_record;
    end
    electrode_positions = zeros(prod(electrode_grid), 2);
    
    % without offset 
    % x_sub_net = mneuro/electrode_grid(1);
    % y_sub_net = nneuro/electrode_grid(2);
    % for i = 1:electrode_grid(1) %x
    %     for j = 1:electrode_grid(2) %y
    %         electrode_positions((i-1)*electrode_grid(2) + j, 1) = (i - 0.5)*x_sub_net;
    %         electrode_positions((i-1)*electrode_grid(2) + j, 2) = (j - 0.5)*y_sub_net;
    %     end
    % end  

    % with offset
    % Define the offset values (you can adjust these as needed)
    x_offset = 2.5;  % example offset in x direction
    y_offset = 2.5;  % example offset in y direction

    % Calculate the spacing between electrodes
    x_sub_net = (mneuro - 2*x_offset) / electrode_grid(1);
    y_sub_net = (nneuro - 2*y_offset) / electrode_grid(2);

    % Compute the electrode positions with the specified offsets
    for i = 1:electrode_grid(1) % x
        for j = 1:electrode_grid(2) % y
            electrode_positions((i-1)*electrode_grid(2) + j, 1) = x_offset + (i - 0.5) * x_sub_net;
            electrode_positions((i-1)*electrode_grid(2) + j, 2) = y_offset + (j - 0.5) * y_sub_net;
        end
    end


end

function [electrode_sensitivity_neuron_E, electrode_sensitivity_neuron_I] =...
           calculate_electrode_sensitivity_neuron(params, electrode_positions, is_stim_or_rec)
    if is_stim_or_rec == 1 %stimulation
        electrode_lambda = params.electrode_lambda_stim;
    else % record
        electrode_lambda = params.electrode_lambda_record;
    end
    mneuro_E = params.mneuro_E;
    nneuro_E = params.nneuro_E;
    mneuro_I = params.mneuro_I;
    nneuro_I = params.nneuro_I;
    electrode_sensitivity_neuron_E = zeros(mneuro_E, nneuro_E, size(electrode_positions, 1));
    electrode_sensitivity_neuron_I = zeros(mneuro_I, nneuro_I, size(electrode_positions, 1));
    % Calculate the scaling factors for the layers
    scale_x_E = max(mneuro_E, mneuro_I) / mneuro_E;
    scale_y_E = max(nneuro_E, nneuro_I) / nneuro_E;
    scale_x_I = max(mneuro_E, mneuro_I) / mneuro_I;
    scale_y_I = max(nneuro_E, nneuro_I) / nneuro_I;
    % Calculate the sensitivity for the excitatory network
    for electrode = 1:size(electrode_positions, 1)
        for i = 1:mneuro_E
            for j = 1:nneuro_E
                d = sqrt((i * scale_x_E - electrode_positions(electrode, 1))^2 + (j * scale_y_E - electrode_positions(electrode, 2))^2);
                electrode_sensitivity_neuron_E(i, j, electrode) = exp(-0.5 * (d / electrode_lambda)^2);
            end
        end
    end
    % Calculate the sensitivity for the inhibitory network
    for electrode = 1:size(electrode_positions, 1)
        for i = 1:mneuro_I
            for j = 1:nneuro_I
                d = sqrt((i * scale_x_I - electrode_positions(electrode, 1))^2 + (j * scale_y_I - electrode_positions(electrode, 2))^2);
                electrode_sensitivity_neuron_I(i, j, electrode) = exp(-0.5 * (d / electrode_lambda)^2);
            end
        end
    end
end