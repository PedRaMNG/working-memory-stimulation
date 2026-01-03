function [auto_electrode_neuro_I_stim_timed, ...
          electrode_neuron_total_rec,...
          integral_electrode_rec_dist, ...
          integral_electrode_rec3] =...
                stimulator_dist(params, model, i, ...
                                electrode_neuron_total_rec, ...
                                auto_electrode_neuro_I_stim_timed, ...
                                integral_electrode_rec_dist, ...
                                integral_electrode_rec3)

    for n = 1 : params.n_electrode_record 
        model.electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
        model.electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
        electrode_neuron_total_rec(n,i) = model.electrode_neuron_E_rec(n,i) + model.electrode_neuron_I_rec(n,i) + model.adjusted_matrix(n, 1);
    end
    % electrode_neuron_total_rec = adjust_initial_values(electrode_neuron_total_rec);

    if rem(i, params.dis_sampling) == 0 ...
        && i > 0 && model.flag_dist(1,i) == 0 ...
        && (params.disruptBothTrainTest == 1 || i < params.test_start_time/params.step)

        rec_sample = i/params.dis_sampling;
        % In-line Integral function -------------------
        integr = zeros(length(1 : params.n_electrode_record), 1);
        for idx = 1:params.n_electrode_record
            n = idx; 
            for j = (i-params.dis_sampling+1):(i-1) % (i-80):(i-1)
                subarea = 0.5 * params.step * (electrode_neuron_total_rec(n, j) + electrode_neuron_total_rec(n, j + 1));
                integr(idx) = integr(idx) + subarea;
            end
        end
        integral_electrode_rec_dist(:, rec_sample) = integr;
        % end of Integral --------------
        
        integral_electrode_rec2 = integral_electrode_rec_dist(:,rec_sample);
        integral_electrode_rec3(:,rec_sample) = (integral_electrode_rec2 > 5); %threshold based on the recoding interval


        % checking the similaries
        for img = 1:numel(params.disruption_pattern)
            comparison = (model.image_target_dist(:,img) == integral_electrode_rec3(:,rec_sample)); % Element-wise comparison
            similarity = sum(comparison(:)) / numel(model.image_target_dist(:,img)) * 100;
            
            % % if matched, activate stimulating electrodes for x second
            if similarity >= 95 
                
                model.dist_timeSteps = [model.dist_timeSteps, i]; % save the time
                
                % Setting stimulating electrode based on the image target
                for j = 1:params.electrode_grid_stim(1)
                    for k = 1:params.electrode_grid_stim(2)
                        % Extract the current 2x2 block
                        block = model.image_target_dist_sqr((2*j-1):(2*j), (2*k-1):(2*k), img);
                        % Check if there is any active recording electrode 
                        has_one = any(block(:) == 1);
                        % Activates related stimulating electrode
                        if has_one
                            stimulating_electrodes(j, k) = 1;
                        else
                            stimulating_electrodes(j, k) = 0;
                        end
                    end
                end

                % stimulating_electrodes = stimulating_electrodes';
                indices = find(stimulating_electrodes == 1);

                Auto_electrode_neuro_I_stim = ...
                sum(model.electrode_sensitivity_neuron_I_stim(:,indices), 2);

                for stim_I_dur = i:i+params.stim_neuro_I_auto_dur
                    auto_electrode_neuro_I_stim_timed(:,stim_I_dur) = ...
                        Auto_electrode_neuro_I_stim .* params.stim_neuro_amp_I;
                    model.flag_dist(1,i:i+params.stim_neuro_I_auto_dur/params.step) = 1; 
                end

            end
        end
    end
end

% plots
% integral_electrode_rec_sqr = reshape(model.integral_electrode_rec_dist, params.electrode_grid_record(1),...
%                               params.electrode_grid_record(2),params.n/params.dis_sampling);
% integral_electrode_rec_sqr = permute(integral_electrode_rec_sqr, [2, 1, 3]);
% show_video(integral_electrode_rec_sqr);

% integral_electrode_rec3_sqr = reshape(model.integral_electrode_rec3, params.electrode_grid_record(1),...
%                               params.electrode_grid_record(2),params.n/params.dis_sampling);
% integral_electrode_rec3_sqr = permute(integral_electrode_rec3_sqr, [2, 1, 3]);
% show_video(integral_electrode_rec3_sqr);

% auto_electrode_neuro_I_stim_timed_sqr = reshape(model.auto_electrode_neuro_I_stim_timed, ...
%                                         params.mneuro_I,params.nneuro_I,params.n);
% auto_electrode_neuro_I_stim_timed_sqr = permute(auto_electrode_neuro_I_stim_timed_sqr, [2, 1, 3]);
% show_video(auto_electrode_neuro_I_stim_timed_sqr);
