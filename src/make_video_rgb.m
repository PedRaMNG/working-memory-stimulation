function [video] = make_video_rgb (Iapp, V_line_E, V_line_I, Record, Glu, Ca, Gli, stimulation, t_record_met, params, section_properties)
    % make_video function that applies different colormaps and ranges to sections

    V_E = reshape(V_line_E, params.mneuro_E, params.nneuro_E, []);
    V_I = reshape(V_line_I, params.mneuro_I, params.nneuro_I, []);
    if params.stimulation_neuro_I == 1
        stimulation = reshape(stimulation, params.mneuro_I, params.nneuro_I, []);
    else
        stimulation = reshape(stimulation, params.mneuro_E, params.nneuro_E, []);
    end
    Record = reshape(Record, params.electrode_grid_record(1), params.electrode_grid_record(2), []);


    Ca_v = Ca(:,:,t_record_met > -1);
    V_v = V_E(:,:,t_record_met > -1);
    Iapp_v = Iapp(:,:,t_record_met > -1);
    Glu_v = Glu(:,:,t_record_met > -1);
    Gli_v = Gli(:,:,t_record_met > -1);


    % Resize V_I and Record and Glu and Gli if needed. In this case, resize V_I and Record to Iapp_v size.
    target_size = [size(Iapp_v, 1), size(Iapp_v, 2)];
    V_I_v = V_I(:,:,t_record_met > -1);
    V_I_v_resized = zeros(target_size(1), target_size(2), size(V_I_v, 3), 'like', V_I_v);
    for t = 1:size(V_I_v, 3)
        V_I_v_resized(:,:,t) = imresize(V_I_v(:,:,t), target_size);
    end
    V_I_v = V_I_v_resized;

    Record_v = Record(:,:,t_record_met > -1);
    Record_v_resized = zeros(target_size(1), target_size(2), size(Record_v, 3), 'like', Record_v);
    for t = 1:size(Record_v, 3)
       Record_v_resized(:,:,t) = imresize(Record_v(:,:,t), target_size);
    end
    Record_v = Record_v_resized;
    Record_v = permute(Record_v, [2,1,3]);


    Glu_v_resized = zeros(target_size(1), target_size(2), size(Glu_v, 3), 'like', Glu_v);
    for t = 1:size(Glu_v, 3)
        Glu_v_resized(:,:,t) = imresize(Glu_v(:,:,t), target_size);
    end
    Glu_v = Glu_v_resized;


    Gli_v_resized = zeros(target_size(1), target_size(2), size(Gli_v, 3), 'like', Gli_v);
    for t = 1:size(Gli_v, 3)
        Gli_v_resized(:,:,t) = imresize(Gli_v(:,:,t), target_size);
    end
    Gli_v = Gli_v_resized;


    % Resize stimulation matrix
    if params.stimulation_neuro_I == 1
        stimulation_v = stimulation(:,:,t_record_met > -1);
        stimulation_v_resized = zeros(target_size(1), target_size(2), size(stimulation_v, 3), 'like', stimulation_v);
        for t = 1:size(stimulation_v, 3)
            stimulation_v_resized(:,:,t) = imresize(stimulation_v(:,:,t), target_size);
        end
        stimulation_v = stimulation_v_resized;
    else
        stimulation_v = stimulation(:,:,t_record_met > -1);
    end


    sections_data = {Iapp_v, V_v, V_I_v, Record_v, Glu_v, Ca_v, Gli_v, stimulation_v}; % reorder sections_data
    processed_sections_rgb = cell(size(sections_data)); % Cell array to store RGB frames

    for i = 1:length(sections_data)
        section = sections_data{i};
        props = section_properties{i};
        cmap = props.cmap;
        limits = double(props.limits);  % Get limits from props
        num_frames = size(section, 3);
        section_rgb_frames = zeros(size(section,1), size(section,2), 3, num_frames, 'uint8'); % Initialize

        for t = 1:num_frames
            frame_2d = section(:,:,t); % Extract 2D frame

            % Check if limits are [0, 0] (or effectively zero)
            if all(abs(limits) < 1e-6)  % Use a small tolerance for floating-point comparison
                % Default normalization (dynamic range)
                section_min = min(frame_2d(:));
                section_range = max(frame_2d(:)) - section_min;
                if section_range == 0
                    normalized_frame = zeros(size(frame_2d));
                else
                    normalized_frame = (frame_2d - section_min) / section_range;
                end
            else
                % Apply normalization using provided limits
                data_min = limits(1);
                data_max = limits(2);
                normalized_frame = (frame_2d - data_min) / (data_max - data_min);
                normalized_frame = max(0, min(1, normalized_frame)); % Clamp to [0, 1]
            end

            % Convert to RGB
            indexed_frame = round(normalized_frame * (size(cmap, 1) - 1)) + 1;
            section_rgb_frames(:,:,:,t) = uint8(ind2rgb(uint8(indexed_frame), cmap) * 255);
        end
        processed_sections_rgb{i} = section_rgb_frames; % Store
        fprintf('Dataset %d / %d\r', i, length(sections_data));
    end

    % Concatenate, now for 2x4 grid
    Iapp_v_rgb = processed_sections_rgb{1};
    V_v_rgb = processed_sections_rgb{2};
    V_I_v_rgb = processed_sections_rgb{3};
    Record_v_rgb = processed_sections_rgb{4};
    Glu_v_rgb = processed_sections_rgb{5};
    Ca_v_rgb = processed_sections_rgb{6};
    Gli_v_rgb = processed_sections_rgb{7};
    stimulation_v_rgb = processed_sections_rgb{8};


    video_frames_rgb = zeros(size(Iapp_v_rgb,1)*2, size(Iapp_v_rgb,2)*4, 3, size(Iapp_v_rgb,4), 'uint8'); % adjust width to *4
    for t = 1:size(Iapp_v_rgb, 4)
        top_row = horzcat(Iapp_v_rgb(:,:,:,t), V_v_rgb(:,:,:,t),  V_I_v_rgb(:,:,:,t), Record_v_rgb(:,:,:,t)); % 4 elements in top row
        bottom_row = horzcat(Glu_v_rgb(:,:,:,t), Ca_v_rgb(:,:,:,t), Gli_v_rgb(:,:,:,t), stimulation_v_rgb(:,:,:,t)); % 4 elements in bottom row
        video_frames_rgb(:,:,:,t) = vertcat(top_row, bottom_row);
    end
    video = video_frames_rgb;
end