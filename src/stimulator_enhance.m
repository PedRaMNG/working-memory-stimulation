function [auto_electrode_neuro_E_stim_timed] = stimulator_enhance(params, model, time_step)
   electrode_neuron_E_rec = zeros(params.n_electrode_record, time_step);
   electrode_neuron_I_rec = zeros(params.n_electrode_record, time_step);
   electrode_neuron_total_rec = zeros(params.n_electrode_record, time_step);
   auto_electrode_neuro_E_stim_timed = zeros(params.mneuro_E * params.nneuro_E, params.n);
   
   for n = 1 : params.n_electrode_record 
      for i = 1:time_step
      electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
      electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
      electrode_neuron_total_rec(n,i) = electrode_neuron_E_rec(n,i) + electrode_neuron_I_rec(n,i);
      end
   end

   electrode_neuron_total_rec = adjust_initial_values(electrode_neuron_total_rec);
   
   %%
   train_patterns = size (params.learn_order,2);
   integral_interval = cell(1, train_patterns);
   for i = 1:train_patterns 
      integral_interval{i} = [model.T_Iapp(i,1), model.T_Iapp(i,2)];
      integral_electrode_rec(:,i) = calculate_integral(params, electrode_neuron_total_rec, integral_interval{i}, [1 : params.n_electrode_record]);
   end

   %% creating first mask from training phase
   for set_idx = 1:train_patterns
      integral_electrode_rec2 = reshape(integral_electrode_rec(:, set_idx), [params.electrode_grid_record(1) params.electrode_grid_record(2)]);
      integral_electrode_rec2 = integral_electrode_rec2';
      max_value = max(integral_electrode_rec2(:)); 
      integral_electrode_rec2 = (integral_electrode_rec2 > (max_value*5)/10);
      integral_electrode_input_mask(:,:,set_idx) = integral_electrode_rec2;
   end

   %% creating second mask from the testing phase
   shift = size(params.learn_order,2);
   train_patterns = size (params.learn_order + shift,2);
   integral_interval = cell(1, train_patterns);
   for i = 1:train_patterns 
      integral_interval{i} = [model.T_Iapp(i+shift,1)+1500, model.T_Iapp(i+shift,2)+1000];
      integral_electrode_rec(:,i) = calculate_integral(params, electrode_neuron_total_rec, integral_interval{i}, [1 : params.n_electrode_record]);
   end

   % Calculate the mean of the max values
   max_values = zeros(1, train_patterns);
   for set_idx = 1:train_patterns
       integral_electrode_rec2 = reshape(integral_electrode_rec(:, set_idx), [params.electrode_grid_record(1), params.electrode_grid_record(2)]);
       integral_electrode_rec2 = integral_electrode_rec2';
       integral_electrode_rec2 = integral_electrode_rec2 .* integral_electrode_input_mask(:, :, set_idx);
       max_values(set_idx) = max(integral_electrode_rec2(:)); 
   end
   mean_max_value = mean(max_values);

   for set_idx = 1:train_patterns
      integral_electrode_rec2 = reshape(integral_electrode_rec(:, set_idx), [params.electrode_grid_record(1) params.electrode_grid_record(2)]);
      integral_electrode_rec2 = integral_electrode_rec2';
      integral_electrode_rec2 = integral_electrode_rec2.*integral_electrode_input_mask(:,:,set_idx);
%       max_value = max(integral_electrode_rec2(:)); 

      integral_electrode_rec2 = ...
                     (integral_electrode_rec2 > (mean_max_value*0.5)/10) ... %0.5
                     & (integral_electrode_rec2 < (mean_max_value*7.5)/10);  %7 7.5
   
      
      integral_electrode_rec3(:,:,set_idx) = integral_electrode_rec2;
   end

   % Sum all
   sum_all_integrals = sum(integral_electrode_rec3, 3);

   % stim resutls

   for j = 1:params.electrode_grid_stim(1)
      for k = 1:params.electrode_grid_stim(2)
          % Extract the current 2x2 block
          block = sum_all_integrals((2*j-1):(2*j), (2*k-1):(2*k));
          
          % Check if there is at least one 4 in the block
          has_two_or_more = any(block(:) >= 2);
          
          % Check if there are at least two 3s in the block
          num_ones = sum(block(:) == 1);
          has_two_or_more_ones = num_ones >= 2;
          
          % Set Total_E_rec_passed based on the conditions
          if has_two_or_more || has_two_or_more_ones
              Total_E_rec_passed(j, k) = 1;
          else
              Total_E_rec_passed(j, k) = 0;
          end
      end
  end

   Total_E_rec_passed = Total_E_rec_passed';
   indices = find(Total_E_rec_passed == 1);
   electrode_sensitivity_neuron_E_stim_sqr = reshape(model.electrode_sensitivity_neuron_E_stim, params.mneuro_E, params.nneuro_E, []);
   imageVar_stim = sum(electrode_sensitivity_neuron_E_stim_sqr(:,:,indices), 3);

   imageVar_stim2 = reshape(imageVar_stim, [], size(imageVar_stim, 3));
    for j = params.stim_start_neuro_E:params.stim_end_neuro_E
          auto_electrode_neuro_E_stim_timed(:,j) = ...
              imageVar_stim2 .* params.stim_neuro_amp_E;
    end
end   


%% Supplimentary Functions
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


function integral_electrode_rec = calculate_integral(params, electrode_neuron_total_rec, range, dimension)
   sampletime = params.step;

   % Initialize integral accumulator
   integr = zeros(length(dimension), 1);
   
   % Compute integral using trapezoidal rule
   for idx = 1:length(dimension)
       n = dimension(idx);
       for i = range(1):(range(2) - 1)
           subarea = 0.5 * sampletime * (electrode_neuron_total_rec(n, i) + electrode_neuron_total_rec(n, i + 1));
           integr(idx) = integr(idx) + subarea;
       end
   end

   % Output the computed integral
   integral_electrode_rec = integr;
end
