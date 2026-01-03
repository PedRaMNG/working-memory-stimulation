%% Video RGB with custom color map for each variable, this function needs a lot of rams
%% for a faster video generation you can run main_video_fast.m
%% Record and stimulation data
% Record
optVideo2.electrode_neuron_E_rec = zeros(params.n_electrode_record, params.n);
optVideo2.electrode_neuron_I_rec = zeros(params.n_electrode_record, params.n);
optVideo2.electrode_neuron_total_rec = zeros(params.n_electrode_record, params.n);

for n = 1 : params.n_electrode_record %opt.electrode_record
    for i = 1:params.n
        optVideo2.electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
        optVideo2.electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
        optVideo2.electrode_neuron_total_rec(n,i) = optVideo2.electrode_neuron_E_rec(n,i) + optVideo2.electrode_neuron_I_rec(n,i);
    end
    fprintf('Electrode: %d / %d\r', n, params.n_electrode_record);
end
%normalizing the values
optVideo2.electrode_neuron_total_rec = adjust_initial_values(optVideo2.electrode_neuron_total_rec);

% Stim
if params.stimulation_mode == 0 % Manual
    if params.stimulation_neuro_I == 1
        optVideo2.electrode_neuro_stim_timed = model.manual_electrode_neuro_I_stim_timed;
    else
        optVideo2.electrode_neuro_stim_timed = model.manual_electrode_neuro_E_stim_timed;
    end
else % Auto
    if params.stimulation_neuro_I == 1
        optVideo2.electrode_neuro_stim_timed = model.auto_electrode_neuro_I_stim_timed;
    else
        optVideo2.electrode_neuro_stim_timed = model.auto_electrode_neuro_E_stim_timed;
    end
end

%%
% Define colormaps and ranges for each section
optVideo2.colrsIn = [255, 255, 255; 237,31,30]/255;
optVideo2.locsIn = [0, 1];
optVideo2.limitsIn = [0, 1];
optVideo2.cmapIn = makeGradient(256, optVideo2.colrsIn, optVideo2.locsIn, 0);
%colrs1 = colrs1(end:-1:1, :); % hot colormap for Iapp


% Voltage E and I
optVideo2.colrsV  = [1 0 0; 1 1 1; 0 0 0];
optVideo2.locsV	 = [-80, -72, 30];
optVideo2.limitsV = [-80, +30];
optVideo2.cmapV   = makeGradient(256, optVideo2.colrsV, optVideo2.locsV, 0);


% Calcium
% optVideo2.locsCa = [0, 0.05,0.11, max(max(model.Ca(:)))];

optVideo2.max_ca = max(max(model.Ca(:)));optVideo2.min_ca = min(min(model.Ca(:)));
optVideo2.locsCa = [optVideo2.min_ca, optVideo2.max_ca*0.3,optVideo2.max_ca*0.7, optVideo2.max_ca];

optVideo2.limitsCa	= [params.ca_0, 0.9];
optVideo2.colrsCa = [255, 255, 255; 29, 67, 80; 237,31,30; 104, 0, 104]/255;
optVideo2.cmapCa   = makeGradient(256, optVideo2.colrsCa, optVideo2.locsCa, 0);


% Record
% optVideo2.cmapRec = parula(256);
% optVideo2.limitsRec = [0, 0];
% optVideo2.limitsRec = [0, max(max(optVideo2.electrode_neuron_total_rec))];
optVideo2.colrsRec = [255, 255, 255; 29, 67, 80; 237,31,30]/255;
optVideo2.locsRec = [0, max(max(optVideo2.electrode_neuron_total_rec))*0.25, max(max(optVideo2.electrode_neuron_total_rec))*0.5];
optVideo2.cmapRec = makeGradient(256, optVideo2.colrsRec, optVideo2.locsRec, 0);
optVideo2.limitsRec = [0, max(max(optVideo2.electrode_neuron_total_rec))*0.5];


% Stim
optVideo2.colrsStim = [255, 255, 255; 29, 67, 80; 237,31,30]/255;
optVideo2.locsStim = [0, 0.5, 1];
optVideo2.cmapStim = makeGradient(256, optVideo2.colrsStim, optVideo2.locsStim, 0);
optVideo2.limitsStim = [0, max(max(optVideo2.electrode_neuro_stim_timed))];


%Gli
optVideo2.colrsGli = [255, 255, 255; 29, 37, 8; 0.6*255,0.1*255,0.2*255]/255;
optVideo2.locsGli = [0, 0.9, 1];
optVideo2.cmapGli = makeGradient(256, optVideo2.colrsGli, optVideo2.locsGli, 0);
optVideo2.limitsGli = [0, 1];

%Neuron-Astrozone Glutamate (Glu)
optVideo2.colrsAZ = [255, 255, 255; 29, 37, 8; 29, 67, 80; 29, 67, 80; 237,31,30; 104, 0, 104]/255;
if max(max(model.neuron_astrozone_activity_Glu(:))) > params.Glu_memorize*1.5
    optVideo2.max_AZ = max(max(model.neuron_astrozone_activity_Glu(:)));
else
    optVideo2.max_AZ = params.Glu_memorize*2;
end
optVideo2.locsAZ = [0,params.Glu_recall_global*0.9, params.Glu_recall_global, params.Glu_memorize*1.1, params.Glu_memorize*1.5, optVideo2.max_AZ];
optVideo2.cmapAZ = makeGradient(256, optVideo2.colrsAZ, optVideo2.locsAZ, 0); % rename cmap to cmapAZ to avoid confusion

optVideo2.limitsAZ = [0, optVideo2.max_AZ]; % rename limits to limitsAZ


video_properties.section_properties = {
    struct('cmap', optVideo2.cmapIn, 'limits', optVideo2.limitsIn),... % Iapp_v
    struct('cmap', optVideo2.cmapV, 'limits', optVideo2.limitsV),...   % V_v
    struct('cmap', optVideo2.cmapV, 'limits', optVideo2.limitsV),...   % V_I_v (using same colormap as V_v)
    struct('cmap', optVideo2.cmapRec, 'limits', optVideo2.limitsRec),...% Record_v
    struct('cmap', optVideo2.cmapAZ, 'limits', optVideo2.limitsAZ),... % Glu
    struct('cmap', optVideo2.cmapCa, 'limits', optVideo2.limitsCa),...  % Ca_v
    struct('cmap', optVideo2.cmapGli, 'limits', optVideo2.limitsGli),... % Gli
    struct('cmap', optVideo2.cmapStim, 'limits', optVideo2.limitsStim) % stimulation_v
};

disp('Compiling RGB video...');
[model.video2] = make_video_rgb(model.Iapp_v_full,... % Iapp_v
                                model.V_line_E,...     % V_v
                                model.V_line_I,...     % V_I_v
                                optVideo2.electrode_neuron_total_rec,... % Record_v
                                model.neuron_astrozone_activity_Glu,... % Glu
                                model.Ca_expand,...    % Ca_v
                                model.Gli_global,...   % Gli
                                optVideo2.electrode_neuro_stim_timed,... % stimulation_v
                                model.T_record_met(1:params.n),...
                                params, video_properties.section_properties);


%%

video_properties.video_fps = 100; % Set desired FPS
video_properties.full_screen = false;
video_properties.figure_position = [100, 100, 1050, 550]; % [left, bottom, width, height]

show_video_rgb(model.video2, video_properties);



%%
% video_properties.colrsAZ = [255, 255, 255; 29, 37, 8; 29, 67, 80; 29, 67, 80; 237,31,30; 104, 0, 104]/255; 
% if max(max(model.neuron_astrozone_activity_Glu(:))) > params.Glu_memorize*1.5
%     video_properties.max_AZ = max(max(model.neuron_astrozone_activity_Glu(:)));
% else
%     video_properties.max_AZ = params.Glu_memorize*2;
% end
% video_properties.locsAZ = [0,params.Glu_recall_global*0.9, params.Glu_recall_global, params.Glu_memorize*1.1, params.Glu_memorize*1.5, video_properties.max_AZ];
% video_properties2.cmap = makeGradient(256, video_properties.colrsAZ, video_properties.locsAZ, 0);
% 
% video_properties2.limits = [0, video_properties.max_AZ]; 
% video_properties2.fps = 100;
% video_properties2.full_screen = false;
% video_properties2.figure_position = [100, 100, 790/3, 550/2]; % [left, bottom, width, height]
% 
% show_video(model.neuron_astrozone_activity_Glu, video_properties2);
%%
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
        adjusted_matrix(i, :) = input_matrix(i, :) + adjustment;
    end
end