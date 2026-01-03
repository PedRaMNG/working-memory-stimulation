function [params, paramss] = show_parameter_gui(params, paramss, param_names)
   % Calculate the required dimensions based on the parameters
   param_count = length(param_names);
   row_height = 30;
   padding = 20;
   button_height = 30;
   min_width = 300;

   % Calculate maximum lengths for parameter names and values
   max_name_length = max(cellfun(@length, param_names));

   % Calculate max value length, handling both numeric and string values
   value_lengths = zeros(size(param_names));
   for i = 1:length(param_names)
         param_value = paramss.(param_names{i});
         if isnumeric(param_value)
            value_lengths(i) = length(num2str(param_value));
         elseif ischar(param_value) || isstring(param_value)
            value_lengths(i) = length(char(param_value));
         end
   end
   max_value_length = max(value_lengths);

   % Calculate widths
   name_width = max(100, max_name_length * 8);  % Approx. 8 pixels per character
   value_width = max(150, max_value_length * 8);
   fig_width = max(min_width, name_width + value_width + 3 * padding);

   fig_height = param_count * row_height + 2 * padding + button_height;

   % Create a figure for the GUI
   fig = figure('Name', 'Parameter Editor', ...
                  'Position', [300, 300, fig_width, fig_height], ...
                  'MenuBar', 'none', ...
                  'ToolBar', 'none', ...
                  'NumberTitle', 'off');

   % Create UI controls for each parameter
   edit_fields = cell(1, param_count);
   for i = 1:param_count
         param_name = param_names{i};
         y_pos = fig_height - (i * row_height) - padding + row_height/2;
         
         % Create label
         uicontrol('Style', 'text', ...
                  'Position', [padding, y_pos, name_width, 20], ...
                  'String', param_name, ...
                  'HorizontalAlignment', 'left');
         
         % Get current parameter value
         param_value = paramss.(param_name);
         
         % Convert value to string for display
         if isnumeric(param_value)
            display_value = num2str(param_value);
         else
            display_value = char(param_value);
         end
         
         % Create edit field
         edit_fields{i} = uicontrol('Style', 'edit', ...
                                    'Position', [name_width + 2*padding, y_pos, value_width, 20], ...
                                    'String', display_value, ...
                                    'Tag', param_name, ...
                                    'UserData', class(param_value));  % Store original data type
   end

   % Add a "Continue" button
   uicontrol('Style', 'pushbutton', ...
               'Position', [fig_width/2 - 50, padding, 100, button_height], ...
               'String', 'Continue', ...
               'Callback', @(~,~) uiresume(fig));

   % Wait for the user to close the GUI
   uiwait(fig);

   % Update params with new values
   for i = 1:param_count
         param_name = param_names{i};
         if ishandle(edit_fields{i})
            new_value_str = get(edit_fields{i}, 'String');
            original_type = get(edit_fields{i}, 'UserData');
            
            % Convert string back to original data type
            try
               switch original_type
                     case 'double'
                        new_value = str2double(new_value_str);
                        if isnan(new_value)
                           error('Invalid numeric value');
                        end
                     case {'char', 'string'}
                        new_value = string(new_value_str);
                     otherwise
                        new_value = new_value_str;
               end
               
               params.(param_name) = new_value;  % Update params
               paramss = rmfield(paramss, param_name);  % Remove from paramss
            catch
               warning(['Invalid input for ' param_name '. Keeping original value.']);
            end
         else
            warning(['Handle for ' param_name ' is not valid. Keeping original value.']);
         end
   end

   % Close the figure if it's still open
   if ishandle(fig)
         close(fig);
   end
end