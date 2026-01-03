function [neuron_astrozone_activity_Glu, I_app_astro] ...
    = get_neuron_astrozone_activity(G, ...
                t_Iapp_met,...
                array_Iapp,...
                I_stim)
    
    params = model_parameters();
    I_app_astro = zeros(params.mastro, params.nastro);
    Total_Glu = reshape(G, params.mneuro_E, params.nneuro_E);
    neuron_astrozone_activity_Glu = zeros(params.mastro, params.nastro);

    if I_stim == 0
        I_stim = zeros(params.mneuro_E, params.nneuro_E, 'uint8');
    else
        I_stim = reshape(I_stim, params.mneuro_E, params.nneuro_E);
    end


    if t_Iapp_met == 0
        Iapp = zeros(params.mneuro_E, params.nneuro_E, 'uint8');
    else
        Iapp = array_Iapp(:, :, t_Iapp_met); % for the timeline of applied input
    end

    for j = 1:params.mastro
        for k = 1:params.nastro
            % Extracting the corresponding 2x2 region
            region = Total_Glu(j:j+1, k:k+1);
            % Calculating the sum of the 2x2 region
            region_sum = sum(region(:));
            % Storing the results
            neuron_astrozone_activity_Glu(j, k) = region_sum;
            
            I_app_astro(j,k) = (sum((Iapp(j:j+1, k:k+1) > 0), 'all') >= 2) || ...
                               (sum((I_stim(j:j+1, k:k+1) > 0.2), 'all') >= 2);

        end
    end

end




