close all;
% clc;

%%  PLOT CONFIGURATION
% Main configuration file for plotting neural network simulation results
% This file controls which plots are generated and their visual settings
%
% Usage: Set values to 0 to disable plotting, or 1-5 to enable with pre-defined window sizes
% Each plot type is organized by the type of data being visualized:
%   - Single neuron/astrocyte plots (individual cell activity)
%   - Network activity maps (spatial-temporal patterns)
%   - Electrode recordings and analysis

%% IMPORTANT:
% TIP 1
% Check the parameters you are trying to plot in init_model.m, they need to
% be saved for params.n.
% TIP 2
% if you faced index error, use this technique
    % len = min(length(Time), length(siGnal(tline)));
    % plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...

%%  SINGLE NEURON/ASTROCYTE PLOTS 

opt.plotNeuronVoltage           = 0;%n		%1% Single Neuron Voltage, Glu, Freq, Input
opt.plotCalIP3                  = 0;%a  	%7% Calcium and IP3
opt.plotIappIsumIsyn            = 0;%n  	%9% Iapp, Isyn_EE (will not plot since it is not saving it on default)

%%  NETWORK ACTIVITY MAPS 

opt.plotIappMAP                 = 0;		%5% Applied current raster plot
opt.plotSpikeRaster             = 0;		%3% Voltage raster plot
opt.plotGluMAP                  = 0;		%4% Glutamate raster plot
    opt.plotGluMAPNeuro         = 0;  % Additional option for neuron-specific glutamate mapping
opt.plotCalMAP                  = 0;		%6% Calcium raster plot
opt.plotG_WM                    = 0;		%12% G_Working Memory
opt.plotIneuroAstro             = 0;		%14% I_neuro-astro (J_omega)
opt.Gli_global                  = 0;        %17% Gli whole-cell

%%  ELECTRODE RECORDINGS & ANALYSIS 

opt.plotElectrodeRecord         = 0;      %2% Electrode Recording
opt.plotElectrodeRecord_integral= 1;      %20% Electrode Recording Integral Analysis - Should be checked
opt.ElectrodesVSdamage          = 0;      %18% Electrodes vs Damage Overlay
opt.plotElectrodesVSinputs      = 0;      %19% Electrodes vs Input Overlay
opt.plot_Damage_VS_Inputs       = 0;      %21% Damage vs Input Comparison
opt.plotAccuracy_elec           = 0;      %22% Correlation based on Electrodes (optimized for the 71x71 network)

%%  3D VISUALIZATIONS IN 2D

opt.plotGluAstroZone3D          = 0;
opt.plotCalcium3D               = 0; 
opt.plotIneuronAstro_train_3D   = 0; 
opt.Gli_global3D                = 0;

%% ===================== some manual plots =====================
% figure();surf(model.neuron_astrozone_activity_Glu(:,:,2240)');
% figure();surf(model.neuron_astrozone_activity_Glu(:,:,5250)');
% figure();imshow(model.neuron_astrozone_activity_Glu(:,:,2500)>3.8);
% figure();surf(model.neuron_astrozone_activity_Glu(:,:,8460)');
% figure();imshow(model.neuron_astrozone_activity_Glu(:,:,8460)>3.5);

% v_test = model.V_line_E(:,2000);
% for e = 1:params.n_electrode_record 
%         electrode_activity_input(:,e) = sum(v_test .* ...
%                                  model.electrode_sensitivity_neuron_E_rec(:,e));  % Calculate activity    
% end
% electrode_activity_input2 = reshape(electrode_activity_input, 8, 8);
% figure();surf(electrode_activity_input2);

%%  GLOBAL PLOT SETTINGS 
% These settings apply to multiple plot types

% Time settings
opt.Time			   = 0;%[2.25 2.6];	% Time window for plots (0 = default: full simulation)

% Visual settings
opt.Opacity            = 0.3; 	% Transparency for overlay images (0-1)

% Integral analysis settings (for electrode signal integration)
opt.integralTime       = [0.7 0.85];
opt.integralTime2      = [0.85 1.038];
opt.integralThreshold1 = 23.8;
opt.integralThreshold2 = 10;

% Electrode configuration
opt.electrode_record = [5,7,10,14,20,40];	% Recording electrode indices
% opt.electrode_record = [1:64];
% opt.electrode_record = [55,56,63,64];

opt.electrode_stim   = [6,7, 9,10,11];		% Stimulating electrode indices
% opt.electrode_stim = params.electrode_neuron_stimulate;
% opt.electrode_stim = [1:16];
% opt.electrode_stim = [16]; 

%%  CELL SELECTION FOR SINGLE CELL PLOTS 
% Specify which individual neurons and astrocytes to plot in single-cell analyses

% Astrocyte selection for single-cell plots [row, column] format
opt.numofAstro	= [1,1; 2,2; 3,3; 10,10; 20,20; 30,30];	

% Neuron selection for single-cell plots [row, column] format
opt.numofNeuro	= [1,1; 2,2; 3,3; 4,4; 5,5; 6,6; 7,7; 8,8; 9,9; 10,10];

%%  PLOT-SPECIFIC VISUAL SETTINGS 
% Color schemes, limits, and other visual parameters for each plot type

% Plot 1 settings - Single Neuron Voltage Analysis
opt.voltageFrequencyBandwidth   = 500;		% Bandwidth samples (Frequency window)
opt.voltageFrequencyShiftShift  = 10;		% Bandwidth Shift Each Sample (Shift frequency window)
opt.voltageFrequencySmoother    = 0.01;		% Smoother of frequency plot

% Plot 3 settings - Voltage Raster Plot
opt.colrs3	= [1 0 0; 1 1 1; 0 0 0];	% Color scheme: [red, white, black]
opt.locs3	= [-80, -72, 30];		% Color transition points
opt.limits3	= [-80, +30];		% Color axis limits

% Plot 4 settings - Glutamate Maps (Neurons and Astrocyte Zones)
% opt.colrs4	= [1 1 1; 
%                 0.9727 0.8242 0.5; 
%                 0.5547 0.6602 0.8555; 
%                 0.1914 0.1797 0.4023];
% opt.locs4	= [0, 0.68, 0.72, 2.6];
opt.limits4	= 0; %[-80, +30];

% opt.colrs4 = [1 1 1; 
%         196/255 195/255 176/255; %
%         ([85 4 0] + 120)/255; % light red
%         54/255 19/255 0/255]; % black red 
% opt.colrs4 = [255, 255, 255;
%              29, 67, 80;
%              243, 144, 100;
%              179, 65, 59 ]/255; %243, 144, 79
opt.colrs4 = [255, 255, 255; % with sexy red colors
             29, 37, 8;
             29, 67, 80;
             237,31,30;
             237,31,36 ]/255; %243, 144, 79

opt.locs4 = [0, 0.68, 0.72,1.73 ,2.9];	% Color transition points

% opt.colrs44 = [255, 255, 255;
%              29, 67, 80; %gray ish
%              127, 139, 75; % green
%              127, 139, 75; % green
%              0, 220, 220;% light blue
%              205, 124, 120; %light red
%              54, 19, 0]/255; %dark red
% opt.locs44 = [0, ...
%             params.Glu_recall_global*0.9,...            
%             params.Glu_recall_global,...
%             params.Glu_memorize, ...
%             params.Glu_memorize*1.1, ...
%             params.Glu_memorize*1.3, ...
%             max(max(model.neuron_astrozone_activity_Glu(:)))];

% Active color scheme for Plot 4 - Astrocyte Zone Glutamate Map
opt.colrs44 = [255, 255, 255; % with sexy red colors
             29, 37, 8; 
             29, 67, 80;
             29, 67, 80; 
             237,31,30; 
             104, 0, 104]/255; 

% Dynamic color limits based on actual data values
if max(max(model.neuron_astrozone_activity_Glu(:))) > params.Glu_memorize*1.5
    max_AZ = max(max(model.neuron_astrozone_activity_Glu(:)));
else
    max_AZ = params.Glu_memorize*2;
end

opt.locs44 = [0, ...
            params.Glu_recall_global*0.9,...            
            params.Glu_recall_global,...
            params.Glu_memorize*1.1, ...
            params.Glu_memorize*1.5, ...
            max_AZ];

% Plot 5 settings - Applied Current Map
opt.colrs5	= [1 1 1; 0.5 0.5 0.5; 0 0 0];	% Color scheme: [white, gray, black]
opt.locs5	= [0, 7, 10];			        % Color transition points
opt.limits5	= 0; %[5, 10];		            % Color axis limits (0 = auto)

% Plot 6 settings - Calcium Map
% opt.colrs6	= [1 1 1; 0.9727 0.8242 0.5; 0.5547 0.6602 0.8555; 0.1914 0.1797 0.4023];
% opt.locs6	= [0, 0.1, 0.3, 0.9];
opt.limits6	= [params.ca_0, 0.9];

% opt.colrs6 = [1 1 1; 
%         196/255 195/255 176/255; %
%         [85/255 4/255 0/255]*1.3; % light red
%         54/255 19/255 0/255]; % black red 
% opt.colrs6 = [255, 255, 255;
%              29, 67, 80;
%              243, 144, 100;
%              179, 65, 59 ]/255;

opt.colrs6 = [255, 255, 255;
             29, 67, 80;
             237,31,30; %85, 4, 0;
             104, 0, 104]/255;

opt.locs6 = [0, 0.05,0.11, max(max(model.Ca(:)))];

% Plot 12 settings - Working Memory Current Map (I_WM)
opt.colrs12	 = [1 1 1; 0.6,0.1,0.2];	% Color scheme: [white, red]
opt.locs12	 = [0, 1];
opt.limits12 = 0; %[5, 10];

% Plot 14 settings - Neuron-to-Astrocyte Current Map (I_neuro_astro)
opt.colrs14	 = [1 1 1; 0.6,0.1,0.2];	% Color scheme: [white, red]
opt.locs14	 = [0, 1];
opt.limits14 = 0; %[5, 10];

% Plot 17 settings - Gliotransmitter Global Map
opt.Gli_global_clr = [255, 255, 255; 
                      29, 37, 8; 
                      0.6*255,0.1*255,0.2*255]/255;
opt.Gli_global_loc = [0, 0.9, 1];

% Plot 18/21 settings - Damage vs Input Comparison
opt.DamVSinpu_clr = [255, 255, 255; 29, 67, 80; 237,31,30]/255;
opt.DamVSinpu_loc = [0, 255/2, 255];

%%  GENERAL MAP COLOR SETTINGS 
% Default color schemes for various map visualizations

opt.MAPcolor1 = [0.5, 0.8, 0.9]; %[Red, Green, Blue] - Light blue-green
opt.MAPcolor2 = [1, 1, 1];       %[Red, Green, Blue] - White

%% Call plot function
plot_functions(model, params, opt);

%%
%     Iastro_neuro = reshape(model.I_WM(1,1,:),1,60000);
%     Post = model.Post_EE;
%     Pre = model.Pre_EE;
%     T_Iapp_met = (model.T_Iapp_met>=1);
%     Iapp = reshape(model.Iapp_v_full(1,1,:),1,60000);
%     V = model.V_line_E(1,:);
%     I_Sum = model.Isumm(1,:);
%     U_E = model.U_line_E(1,:);
%     Glu = model.G(1,:);
%     Isyn_EE = model.Isyn_line_EE(1,:);
    
%% --------------------


% timess = T_Iapp(:,1) + 500;
% for i = timess
% 	pp = 1 - double(model.video(:,:,i))/255;
% 	aa(:,:,1) = pp * 101;
% 	aa(:,:,2) = pp * 167;
% 	aa(:,:,3) = pp * 151;
% 	image(aa/255);
%     axis equal;
%     xlim([-10, 250]);
% 	input('');
% end


%% ST map

%     connection = zeros(params.quantity_neurons_E, params.quantity_neurons_E);
%     Ppre = reshape(model.Pre_EE, params.N_connections_EE, params.quantity_neurons_E);
%     Ppost = reshape(model.Post_EE, params.N_connections_EE, params.quantity_neurons_E);
%     tumVec = (double(model.Damage_E(:))/255)*0 + 1;
% 
%     for i = 1:params.quantity_neurons_E
%         connection(Ppost(:, i), i) = tumVec(Ppost(:, i));
%     end
%     
%     %connection = 1 - connection;
%     
%     colrs = [255, 255, 255;         % red gradient ------
%          141, 168, 218;
%          0, 0, 0]/255;
%     locs = [0,                          ...
%             0.6*max(connection(:)),     ...
%             max(connection(:))];
%     pro.limits  = [min(connection(:)), max(connection(:))];
%     pro.cmap    = makeGradient(256, colrs, locs, 0);
%     
%     figure();
%         plot(0,0); hold on;
%         imagesc(connection);
%         colormap(pro.cmap);
%         set(gca,'XTick',[], 'ZTick', [], 'YTick', []);
%         xlim([0, params.quantity_neurons_E]);
%         ylim([0, params.quantity_neurons_E]);
        

%% Histogram of ST

%     figure();
%     histogram(model.ST_S_EE);



%% plot input patterns

%     vv = double(model.Iapp(:,:,8)) / 8;
%     img(:,:,1) = (1 - vv * 0.4);
%     img(:,:,2) = (1 - vv * 0.9);
%     img(:,:,3) = (1 - vv * 0.9);
%     
%     figure();
%     % plot(0,0); hold on;
%     image(img);
%     
%     axis equal;
%     xlim([0.5 79.5]);
%     ylim([0.5 79.5]);


%% Plot map of signals with custom gradiant

%     HH	= zeros(params.mastro*params.nastro, params.n);
%     k	= 1;
%     for j = 1:params.nastro
%         for i = 1:params.mastro
%             HH(k,:) = reshape(model.Ca(i,j,:), [1, params.n]);
%             k = k + 1;
%         end
%     end

%% countor
% zmin = floor(min(differenceS(:))); 
% zmax = ceil(max(differenceS(:)));
% zinc = (zmax - zmin) / 8; % 5, 8
% zlevs = zmin:zinc:zmax;
% figure();contour(differenceS,zlevs);


%% --------------------------------------

% num_with    = memory_performance_num_with.freq_images;
% num_witho   = memory_performance_num_witho.freq_images;
% num_diff    = num_with - num_witho;
% 
% strip_with  = memory_performance_strip_with.freq_images;
% strip_witho = memory_performance_strip_witho.freq_images;
% strip_diff  = strip_with - strip_witho;
% 
% vert_with   = memory_performance_vert_with.freq_images;
% vert_witho  = memory_performance_vert_witho.freq_images;
% vert_diff   = vert_with - vert_witho;


%% Frequency
% shots = sum(vert_with, 3);
% colrs = [255, 255, 255;
%              110, 109, 98;
%              29, 67, 80;
%              243, 144, 79]/255; % 204, 17, 19
% locs = [0,      ...
%             0.30,   ...
%             0.7,   ...
%             1];

% colrs = [49, 46, 103;
%          141, 168, 218;
%          255, 255, 255;
%          250, 70, 70;
%          54,  19, 0]/255;
% 
%     locs = [min(shots(:)),      ...
%             min(shots(:))*0.2,  ...
%             0,                  ...
%             max(shots(:))*0.2,	...
%             max(shots(:))];
% pro.limits  = [min(shots(:)), max(shots(:))];
% pro.cmap    = makeGradient(256, colrs, locs, 0);
    
% for ii = 1:size(shots, 3)
%     figure();
%         x = 1:79;
%         y = 1:79;
%         [X,Y] = meshgrid(x,y);
%         surf(X,Y, shots(end:-1:1,:,ii));
%         set(gca,'XTick',[], 'ZTick', [], 'YTick', []);
%         shading interp;
%         colormap(pro.cmap);
%         view(0,90);
%         colorbar;
%         caxis(pro.limits);
% %         axis equal;
%         axis([0, size(shots, 1)+1, 0, size(shots, 2)+1]);
% end

%%
% All_stim_electrodes_neuro = sum(model.electrode_sensitivity_neuron_E_stim, 3);
% All_rec_electrodes_neuro = sum(model.electrode_sensitivity_neuron_E_rec, 3);
% figure(); surf(All_stim_electrodes_neuro'); view(90,90); title('All Stimulating Electrodes');
% figure(); surf(All_rec_electrodes_neuro'); view(90,90); title('All Stimulating Electrodes');
% imwrite(All_stim_electrodes_neuro, 'SavedData/all_stim_electrodes_43x43.jpg');
% imwrite(All_rec_electrodes_neuro, 'SavedData/all_rec_electrodes_43x43.jpg');