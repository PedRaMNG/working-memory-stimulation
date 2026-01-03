function [Ca, h, IP3, I_WM, array_I_neuro, Gli_global, ...
    neuron_astrozone_input_timings] = ...
    step_astrocytes(neuron_astrozone_activity_Glu, ...
                    array_I_neuro, i, Ca, h, IP3, I_WM, Gli_global,...
                    neuron_astrozone_input_timings, I_app_astro)
	
    params          = model_parameters();
    diffusion_Ca    = zeros(params.mastro, params.mastro, 'double');
    diffusion_IP3   = zeros(params.mastro, params.mastro, 'double');

	for j = 1 : params.mastro
        for k = 1 : params.nastro
            if (neuron_astrozone_input_timings(j,k,i)) && (neuron_astrozone_activity_Glu(j, k) >= params.Glu_memorize)
                
                neuron_astrozone_input_timings(j,k, i:i+params.input_window) = 0;
                shift = fix(params.t_neuro_astro / params.step) - 1;
                array_I_neuro(j, k, i : i + shift) = params.amplitude_neuro_astro;

            end
            

            % Computes diffusion of calcium and IP3 in astrocyte network
            if (j == 1) && (k == 1)                             % Corner top left
                diffusion_Ca  = Ca(j + 1, k)  + Ca(j,k + 1)  - 2 * Ca(j,k);
                diffusion_IP3 = IP3(j + 1, k) + IP3(j,k + 1) - 2 * IP3(j,k);
            elseif (j == params.mastro) && (k == params.nastro) % Corner bottom right
                diffusion_Ca  = Ca(j - 1, k)  + Ca(j,k - 1)  - 2 * Ca(j,k);
                diffusion_IP3 = IP3(j - 1, k) + IP3(j,k - 1) - 2 * IP3(j,k);
            elseif (j == 1) && (k == params.nastro)             % Corner top right
                diffusion_Ca  = Ca(j + 1, k)  + Ca(j,k - 1)  - 2 * Ca(j,k);
                diffusion_IP3 = IP3(j + 1, k) + IP3(j,k - 1) - 2 * IP3(j,k);
            elseif (j == params.mastro) && (k == 1)             % Corner bottom left
                diffusion_Ca  = Ca(j - 1, k)  + Ca(j,k + 1)  - 2 * Ca(j,k);
                diffusion_IP3 = IP3(j - 1, k) + IP3(j,k + 1) - 2 * IP3(j,k);
            elseif j == 1                                       % First top row
                diffusion_Ca  = Ca(j + 1, k)  + Ca(j, k - 1)  + Ca(j,k + 1)  - 3 * Ca(j,k);
                diffusion_IP3 = IP3(j + 1, k) + IP3(j, k - 1) + IP3(j,k + 1) - 3 * IP3(j,k);
            elseif j == params.mastro                           % Last bottom row
                diffusion_Ca  = Ca(j - 1, k)  + Ca(j, k - 1)  + Ca(j,k + 1)  - 3 * Ca(j,k);
                diffusion_IP3 = IP3(j - 1, k) + IP3(j, k - 1) + IP3(j,k + 1) - 3 * IP3(j,k);
            elseif k == 1                                       % First left column
                diffusion_Ca  = Ca(j - 1, k)  + Ca(j + 1, k)  + Ca(j,k + 1)  - 3 * Ca(j,k);
                diffusion_IP3 = IP3(j - 1, k) + IP3(j + 1, k) + IP3(j,k + 1) - 3 * IP3(j,k);
            elseif k == params.nastro                           % Last right column
                diffusion_Ca  = Ca(j - 1, k)  + Ca(j + 1, k)  + Ca(j,k - 1)  - 3 * Ca(j,k);
                diffusion_IP3 = IP3(j - 1, k) + IP3(j + 1, k) + IP3(j,k - 1) - 3 * IP3(j,k);
            elseif (j > 1) && (j < params.mastro) && (k > 1) && (k < params.nastro) % remaining columns and ...
                                                                                    % rows in the middle
                diffusion_Ca  = Ca(j - 1, k)  + Ca(j + 1, k)  + Ca(j, k - 1)  + Ca(j,k + 1)  - 4 * Ca(j,k);
                diffusion_IP3 = IP3(j - 1, k) + IP3(j + 1, k) + IP3(j, k - 1) + IP3(j,k + 1) - 4 * IP3(j,k);
            end

            %% Astrocyte model
            X			= [Ca(j, k) h(j, k) IP3(j, k)];
            I_neuro		= array_I_neuro(j, k, i);
            % Solving astrocyte equasions using Rung-Kutta method in "runge_astro" function:
            
            w1		= runge_astro(0, X,                      I_neuro, diffusion_Ca, diffusion_IP3, Gli_global(j,k));
            w2		= runge_astro(0, X + params.u2   .* w1', I_neuro, diffusion_Ca, diffusion_IP3, Gli_global(j,k));
            w3		= runge_astro(0, X + params.u2   .* w2', I_neuro, diffusion_Ca, diffusion_IP3, Gli_global(j,k));
            w4		= runge_astro(0, X + params.step .* w3', I_neuro, diffusion_Ca, diffusion_IP3, Gli_global(j,k));

            X			= X + params.u6 .* (w1' + 2 .* w2' + 2 .* w3' + w4');
            Ca(j, k)	= X(1);
            h(j, k)		= X(2);
            IP3(j, k)	= X(3);
            			
            %% Astrocyte detection of synaptic events
            bnh = rem(i, params.shift_window_astro_watch);	% astrocytes' timestep
            if (Ca(j, k) > params.ca_threshold_global) ...
                && (bnh == 0) ...
                && (neuron_astrozone_activity_Glu(j, k) >= params.Glu_recall_global) ...
                && (I_app_astro(j,k)) 
                
                I_WM(j, k, i : i + params.impact_astro) = 1; %This timing is constantly being updated
            end
            
            
            % Whole-cell gliotransmitter for WM 
            if (Ca(j, k) > params.ca_threshold_global) && (I_WM(j,k,i)) && (neuron_astrozone_activity_Glu(j, k) >= params.Glu_recall_global)
                % gliotransmitter release from astrocytes
                Gli_global(j,k) = Gli_global(j,k) + params.step*(params.r_Gli_global - (Gli_global(j,k)/params.tau_Gli)); 
            else
                % gliotransmitter decay rate
                Gli_global(j,k) = Gli_global(j,k) + params.step*(-1*(Gli_global(j,k)/(params.tau_Gli*0.4)));      
            end
            Gli_global(j,k) = min(Gli_global(j,k), 1);

        end
	end
end