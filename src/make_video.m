function [video] = make_video (Ca, V_line_E, V_line_I, Iapp, t_record_met, Record, stimulation_E, params)

    V_E = reshape(V_line_E, params.mneuro_E, params.nneuro_E, []);
    V_I = reshape(V_line_I, params.mneuro_I, params.nneuro_I, []);
    if params.stimulation_neuro_I == 1
        stimulation_E = reshape(stimulation_E, params.mneuro_I, params.nneuro_I, []); 
    else
        stimulation_E = reshape(stimulation_E, params.mneuro_E, params.nneuro_E, []); 
    end
    Record = reshape(Record, params.electrode_grid_record(1), params.electrode_grid_record(2), []);
    
    
    
    % Iapp = Iapp ./ 10;
    % V_E = V_E ./ 100;
    Ca_v = Ca(:,:,t_record_met > -1);
    Ca1 = min(Ca_v,[],'all');
    Ca_v = (Ca_v - Ca1);
    Ca1 = max(Ca_v,[],'all');
    Ca_v = Ca_v ./ Ca1 .* 255;
    Ca_v = uint8(Ca_v);
    
    V_v = V_E(:,:,t_record_met > -1);
    V_v1 = min(V_v,[],'all');
    V_v = (V_v - V_v1);
    V_v1 = max(V_v,[],'all');
    V_v = V_v ./ V_v1 .* 255;
    V_v = uint8(V_v);
    
    Iapp_v = Iapp(:,:,t_record_met > -1);
    Iapp_v = single(Iapp_v);
    Iapp_v1 = min(Iapp_v,[],'all');
    Iapp_v = (Iapp_v - Iapp_v1);
    Iapp_v1 = max(Iapp_v,[],'all');
    Iapp_v = Iapp_v ./ Iapp_v1 .* 255;
    Iapp_v = uint8(Iapp_v);
    
    % Resize V_I to match the size of the other variables
    target_size = [size(Iapp_v, 1), size(Iapp_v, 2)];
    V_I_v = V_I(:,:,t_record_met > -1);
    V_I_v_resized = zeros(target_size(1), target_size(2), size(V_I_v, 3), 'like', V_I_v);
    for t = 1:size(V_I_v, 3)
        V_I_v_resized(:,:,t) = imresize(V_I_v(:,:,t), target_size);
    end
    V_I_v = V_I_v_resized;
    
    % Normalize V_I
    V_I_v1 = min(V_I_v,[],'all');
    V_I_v = (V_I_v - V_I_v1);
    V_I_v1 = max(V_I_v,[],'all');
    V_I_v = V_I_v ./ V_I_v1 .* 255;
    V_I_v = uint8(V_I_v);
    
    % Resize Record to match the size of the other variables
    target_size = [size(Iapp_v, 1), size(Iapp_v, 2)];
    Record_v = Record(:,:,t_record_met > -1);
    Record_v_resized = zeros(target_size(1), target_size(2), size(Record_v, 3), 'like', Record_v);
    for t = 1:size(Record_v, 3)
       Record_v_resized(:,:,t) = imresize(Record_v(:,:,t), target_size);
    end
    Record_v = Record_v_resized;
    
    % Process Record Data (Record)
    Record_v1 = min(Record_v,[],'all');
    Record_v = (Record_v - Record_v1);
    Record_v1 = max(Record_v,[],'all');
    Record_v = Record_v ./ Record_v1 .* 255;
    Record_v = uint8(Record_v);
    Record_v = permute(Record_v, [2,1,3]);
    
    % Resize stimulation_E matrix
    if params.stimulation_neuro_I == 1
        stimulation_E_v = stimulation_E(:,:,t_record_met > -1);
        stimulation_E_v_resized = zeros(target_size(1), target_size(2), size(stimulation_E_v, 3), 'like', stimulation_E_v);
        for t = 1:size(stimulation_E_v, 3)
            stimulation_E_v_resized(:,:,t) = imresize(stimulation_E_v(:,:,t), target_size);
        end
        stimulation_E_v = stimulation_E_v_resized;

        stimulation_E_v1 = min(stimulation_E_v,[],'all');
        stimulation_E_v = (stimulation_E_v - stimulation_E_v1);
        stimulation_E_v1 = max(stimulation_E_v,[],'all');
        stimulation_E_v = stimulation_E_v ./ stimulation_E_v1 .* 255;
        stimulation_E_v = uint8(stimulation_E_v);

    else
        % Process stimulation_E Data (stimulation_E)
        stimulation_E_v = stimulation_E(:,:,t_record_met > -1);
        stimulation_E_v1 = min(stimulation_E_v,[],'all');
        stimulation_E_v = (stimulation_E_v - stimulation_E_v1);
        stimulation_E_v1 = max(stimulation_E_v,[],'all');
        stimulation_E_v = stimulation_E_v ./ stimulation_E_v1 .* 255;
        stimulation_E_v = uint8(stimulation_E_v);
    end
    
    
    % Combine all processed data into video frames with two rows
    top_row = horzcat(Iapp_v, V_v, Ca_v);
    bottom_row = horzcat(Record_v, V_I_v, stimulation_E_v);
    video = vertcat(top_row, bottom_row);
    end