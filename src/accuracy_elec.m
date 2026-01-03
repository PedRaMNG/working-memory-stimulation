function [similarity, tline, target] = accuracy_elec(params, model, patTers)
   %% Initialization 
   % this code is optimized for the 71x71 network.
   % Step 1
   electrode_neuron_E_rec = zeros(params.n_electrode_record, params.n);
   electrode_neuron_I_rec = zeros(params.n_electrode_record, params.n);
   electrode_neuron_total_rec = zeros(params.n_electrode_record, params.n);
   for n = 1 : params.n_electrode_record 
         for i = 1:params.n
         electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
         electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
         electrode_neuron_total_rec(n,i) = electrode_neuron_E_rec(n,i) + electrode_neuron_I_rec(n,i);
         end
         disp(['Electrode: ', num2str(n), ' / ', num2str(params.n_electrode_record)]);
   end
   electrode_neuron_total_rec = adjust_initial_values(electrode_neuron_total_rec);

   % step 2
   similarity = zeros(params.n, length(patTers));
   target = zeros(params.n_electrode_record, numel(patTers));
   image_number = 1;
   for imG = patTers
      inputs = double(cell2mat(model.images(imG))) / 255;
      inputs = reshape(inputs, params.mneuro_E*params.nneuro_E,1);
      inputs = 1-inputs;
      for e = 1:params.n_electrode_record 
            electrode_activity = sum(inputs .* model.electrode_sensitivity_neuron_E_rec(:,e));  % Calculate activity
            if electrode_activity > 20
               target(e,image_number) = 1;  % Mark the electrode as active (1)
            else
               target(e,image_number) = 0;
            end
      end
      image_number = image_number + 1;
   end

   %% Similarity calculaiton
   disp('Processing time steps ...');
   window = 80;      % Set the window size
   shiftWindow = 20; % Set the shift step
   cout = 1;
   for i = (window + 1):shiftWindow:params.n
      % Integral  -------------------
      integr = zeros(length(1 : params.n_electrode_record), 1);
      for idx = 1:params.n_electrode_record
         n = idx; 
         for j = (i - window):(i - 1)
               subarea = 0.5 * params.step * (electrode_neuron_total_rec(n, j) + electrode_neuron_total_rec(n, j + 1));
               integr(idx) = integr(idx) + subarea;
         end
      end
      % Store the integral results for each window slide
      integral_electrode_rec_dist(:, cout) = integr;
      % Process the calculated integral
      integral_electrode_rec2 = integral_electrode_rec_dist(:, cout);
      integral_electrode_rec3(:, cout) = (integral_electrode_rec2 > 5);
      % checking the similarities
      for img = 1:numel(patTers)
         comparison = (target(:, img) == integral_electrode_rec3(:, cout)); % Element-wise comparison
         similarity(cout, img) = sum(comparison(:)) / numel(target(:, img)) * 100;
      end
      cout = cout + 1;  % Increment counter for storing results
   end

   tline = 1:size(similarity,1);
   %
   internal_plots1 = 0;
   internal_plots2 = 0;
   if internal_plots1
      %% plot similarity
      figure();
      vectime = [0, params.t_end*1000]
      hold on;
      for k = 1:length(patTers)
            siGnal = similarity(:,k);
            plot(tline*params.step*10, siGnal, ...
               'Linewidth', 0.5, ...
               'DisplayName', num2str(k));
      end
      ax = gca; ax.YGrid = 'on';
      ylabel('Accuracy');
      ylim([55, 110]);
      xlim(vectime);
      ax = gca; ax.YGrid = 'on';
      lgd = legend;
      title(lgd, 'image');
      xlabel('Time (s)');
      disp('Finished!');
   end
   if internal_plots2
      %% plot Target
      target_sqr = reshape(target,...
                  params.electrode_grid_record(1), ...
                  params.electrode_grid_record(2), ...
                  numel(patTers));

      target_sqr = permute(target_sqr, [2, 1, 3]);
      num_patterns = numel(patTers);
      % Determine the subplot grid size
      subplot_rows = ceil(sqrt(num_patterns));
      subplot_cols = ceil(num_patterns / subplot_rows);
      figure('Name', 'Electrode Grid Visualization', 'Position', [100, 100, 1200, 800]);
      for ii = 1:num_patterns
         subplot(subplot_rows, subplot_cols, ii);
         imagesc(squeeze(target_sqr(:,:,ii)));
         colormap(gca, 'gray');  % Use grayscale colormap
         axis equal tight;
         title(['Pattern ', num2str(patTers(ii))]);
         xlabel('Column');
         ylabel('Row');
         % Add colorbar
         c = colorbar;
         c.Label.String = 'Activation';
         % Optionally, you can add grid lines to better visualize individual electrodes
         hold on;
         for j = 0.5:1:params.electrode_grid_record(2)+0.5
            line([j j], [0.5 params.electrode_grid_record(1)+0.5], 'Color', 'r', 'LineStyle', ':');
         end
         for j = 0.5:1:params.electrode_grid_record(1)+0.5
            line([0.5 params.electrode_grid_record(2)+0.5], [j j], 'Color', 'r', 'LineStyle', ':');
         end
         hold off;
      end
      % Adjust the layout
      sgtitle('Electrode Grid Activation Patterns');
      set(gcf, 'Color', 'w');  % Set figure background to white
   end

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