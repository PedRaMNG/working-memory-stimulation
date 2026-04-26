function ress = plot_functions(model, params, opt)

    %% =============================== MAIN PLOTTING FUNCTION ===============================
    % This function generates various visualizations of neural network simulation results
    % based on the configuration options specified in the opt structure
    %
    % Inputs:
    %   model - Structure containing simulation data (voltages, currents, etc.)
    %   params - Structure containing simulation parameters
    %   opt - Structure containing plot configuration options
    %
    % Output:
    %   ress - Return value (currently always 1)

    %% ===================== EXTRACT PLOT CONFIGURATION OPTIONS =====================
    % Extract all plot configuration flags from the opt structure
    % These determine which plots will be generated

    % Single neuron/astrocyte plot flags
    plotNeuronVoltage       = opt.plotNeuronVoltage;      %1% Single Neuron Voltage and Glutamate
    plotCalIP3              = opt.plotCalIP3;             %7% Calcium and IP3
    plotIappIsumIsyn        = opt.plotIappIsumIsyn;       %9% Iapp, Isyn_EE, Isum

    % Network activity map plot flags
    plotSpikeRaster         = opt.plotSpikeRaster;        %3% Voltage Map
    plotGluMAP              = opt.plotGluMAP;             %4% Glutamate Map
    plotIappMAP             = opt.plotIappMAP;            %5% I applied Map
    plotCalMAP              = opt.plotCalMAP;             %6% Calcium Map and mean of that
    plotG_WM                = opt.plotG_WM;               %12% I working memory
    plotIneuroAstro         = opt.plotIneuroAstro;        %14% I neuro astro
    plotGli_global          = opt.Gli_global;             %17% Gli global

    % Electrode recording and analysis plot flags
    plotElectrodeRecord     = opt.plotElectrodeRecord;    %2% 
    plotElectrodeRecord_integral = opt.plotElectrodeRecord_integral; %20%
    plotElectrodesVSdamage  = opt.ElectrodesVSdamage;     %18% Electrodes vs damage overlay
    plotElectrodesVSinputs  = opt.plotElectrodesVSinputs; %19% Electrodes vs inputs overlay
    plot_Damage_VS_Inputs   = opt.plot_Damage_VS_Inputs;  %21% Damage vs inputs comparison
    plotAccuracy_elec       = opt.plotAccuracy_elec;      %22% correlation with electrodes

    % 3D visualization plot flags
    plotGluAstroZone3D      = opt.plotGluAstroZone3D;    %23% 3D Glutamate astrozone
    plotCalcium3D           = opt.plotCalcium3D;         %24% 3D Calcium
    plotIneuronAstro_train_3D = opt.plotIneuronAstro_train_3D; %25% 3D Training visualization
    plotGli_global3D        = opt.Gli_global3D;          %26% 3D Gli global

    %% ===================== TIME AND DATA SETUP =====================
    % Initialize time vectors and extract common plotting parameters

    vectime = opt.Time;                   % Gets time interval from configuration

    % Time variables for plotting
    Time = model.T(2:end);               % Time vector (excluding initial point)
    tline = 1:params.n;                 % Time indices for full simulation

    % Handle time window specification
    if vectime == 0
        vectime = [0, params.t_end];     % Default: full simulation time
    end

    % Convert time to sample indices
    vectimen = floor(vectime / params.step);

    % Extract single-cell plotting parameters
    % Part 1 vars - Single Neuron voltage and Glutamate analysis
    voltageFrequencyBandwidth = opt.voltageFrequencyBandwidth;
    voltageFrequencyBandShift = opt.voltageFrequencyShiftShift;
    voltageFrequencySmoother = opt.voltageFrequencySmoother;

    % Part 7, 8 vars: Astrocyte selection [j, k]
    numofAstro = opt.numofAstro;

    % Part 1 & 9 vars: Neuron selection [j, k]
    numofNeuro = opt.numofNeuro;

    % Map color settings
    MAPcolor1 = opt.MAPcolor1;
    MAPcolor2 = opt.MAPcolor2;

    % Visual settings
    backOpacity = opt.Opacity;

    %% ===================== PLOTTING STYLE SETUP =====================
    % Define colors, line styles, and create background images for consistent styling

    % Color palette for multiple line plots (R-G-B format)
    clrs = [0, 114, 189;                 % Blue
                    217, 83, 25;                 % Orange
                    237, 177, 32;                 % Yellow
                    126, 47, 142;                 % Purple
                    119, 172, 48;                 % Green
                    77, 190, 238;                 % Cyan
                    162, 20, 47] / 255;           % Red
    noc = size(clrs, 1);                  % Number of colors

    % Line style options for differentiating multiple traces
    linestyle = {'-', ...                 % Solid line
                    '--', ...                 % Dashed line
                    '-.', ...                 % Dash-dotted line
                    ':'}; ...                 % Dotted line
    nols = numel(linestyle);              % Number of line styles

    % Create background images for damage overlay (for excitatory and inhibitory neurons)
    dam = damageBackground(model.Damage_E, params.n);
    damage = makeImage(dam, MAPcolor2);
    dam_I = damageBackground(model.Damage_I, params.n);
    damage_I = makeImage(dam_I, MAPcolor2);
    %% ===================== SINGLE CELL ANALYSIS PLOTS =====================
    % These plots show detailed activity from individual neurons or astrocytes

    %% Part 1 - Single Neuron Analysis: Voltage, Glutamate, and Frequency
    % Shows detailed electrical activity and neurotransmitter levels for selected neurons
    if plotNeuronVoltage ~= 0
        disp('Started - plotNeuronVoltage');
        ssi = setsizee(plotNeuronVoltage);
        figure(); ssi();
        
        % Subplot 1: Membrane Voltage
        subplot(311);
        hold on;
        plotSecAreas(model, params.step, [-200, 200]);
        miny = 0; maxy = 0;

        % Plot voltage traces for each selected neuron
        for k = 1:size(numofNeuro, 1)
            neurInd = sub2ind([params.mneuro_E, params.nneuro_E], numofNeuro(k, 1), numofNeuro(k, 2));
            siGnal = model.V_line_E(neurInd, :);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 0.5, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofNeuro(k, 1)), ', ', num2str(numofNeuro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        ylabel('Voltage (mV)');
        ylim([1.1 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
        
        % Subplot 2: Glutamate Concentration
        subplot(312);
        hold on;
        plotSecAreas(model, params.step, [-100, 100]);
        miny = 0; maxy = 0;
        neuron_astrozone_activity_Glu_sqr = reshape(model.neuron_astrozone_activity_Glu, (params.mneuro_E-1)*(params.nneuro_E-1), params.n);
        
        % Plot glutamate traces for each selected neuron
        for k = 1:size(numofNeuro, 1)
            neurInd = sub2ind([params.mneuro_E, params.nneuro_E], numofNeuro(k, 1), numofNeuro(k, 2));
            siGnal = neuron_astrozone_activity_Glu_sqr(neurInd, :);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal(tline)));
            len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 1, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofNeuro(k, 1)), ', ', num2str(numofNeuro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        ylabel('Glutamate');
        ylim([0.8 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
        
        % Subplot 3: Firing Frequency Analysis
        subplot(313);
        hold on;
        plotSecAreas(model, params.step, [-100, 100]);
        miny = 0; maxy = 0;

        % Calculate and plot firing frequency for each selected neuron
        for k = 1:size(numofNeuro, 1)
            neurInd = sub2ind([params.mneuro_E, params.nneuro_E], numofNeuro(k, 1), numofNeuro(k, 2));
            [~, siGnal] = onlinefreq(model.V_line_E(neurInd, :)', ...
                voltageFrequencyBandwidth, ...
                voltageFrequencyBandShift, ...
                params.step, ...
                voltageFrequencySmoother, ...
                params.neuron_fired_thr); siGnal = siGnal';

            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 1, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofNeuro(k, 1)), ', ', num2str(numofNeuro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        lgd = legend;
        title(lgd, 'Neurons(m,n)');
        ylabel('Fire frequency (Hz)');
        ylim([0.8 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
        xlabel('Time (s)');
        disp('Finished!');
    end

    %% ===================== ELECTRODE RECORDING ANALYSIS =====================
    % These plots show electrode recordings and related analyses

    %% Part 2 - Electrode Recording Analysis
    if plotElectrodeRecord ~= 0 
        disp('Started - plotElectrodeRecord');
        ssi = setsizee(plotElectrodeRecord);
        %-----
        if params.stimulation_mode %automatic stimulation
            % electrode_neuron_E_rec = model.electrode_neuron_E_rec;
            % electrode_neuron_I_rec = model.electrode_neuron_I_rec;
            % electrode_neuron_total_rec ...
            %     = electrode_neuron_E_rec ...
            %     + electrode_neuron_I_rec ...
            %     + model.adjusted_matrix(:, 1);
            electrode_neuron_total_rec = model.electrode_neuron_total_rec;
        else %manual stimulation
            electrode_neuron_E_rec = zeros(params.n_electrode_record, params.n);
            electrode_neuron_I_rec = zeros(params.n_electrode_record, params.n);
            electrode_neuron_total_rec = zeros(params.n_electrode_record, params.n);
            for n = 1 : params.n_electrode_record %opt.electrode_record
                for i = 1:params.n
                    electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
                    electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
                    electrode_neuron_total_rec(n,i) = electrode_neuron_E_rec(n,i) + electrode_neuron_I_rec(n,i);
                    % electrode_neuron_total_rec(n,i) = electrode_neuron_E_rec(n,i) + electrode_neuron_I_rec(n,i) + model.adjusted_matrix(n, 1);
                end
            end
            %normalizing the values
            electrode_neuron_total_rec = adjust_initial_values(electrode_neuron_total_rec);

        end
        

        %-----
        xMat = repmat(tline', 1, numel(opt.electrode_record)); 
        y = 0:0.001:0.001*(numel(opt.electrode_record)-1);
        yMat = repmat(y, numel(tline'), 1); 
        % zMat = electrode_neuron_E_rec(opt.electrode_record,:)'; 
        zMat = electrode_neuron_total_rec(opt.electrode_record,:)'; 
        
        %---
        figure(); ssi();
        subplot(231);
        electrode_sensitivity_neuron_E_rec_sqr = reshape(model.electrode_sensitivity_neuron_E_rec, params.mneuro_E, params.nneuro_E, []);
        imageVar1 = sum(electrode_sensitivity_neuron_E_rec_sqr(:,:,opt.electrode_record), 3);
        imageVar01 = 1 - imageVar1;
        colormap(hot); % Apply the colormap to the current figure
        imageVar01_colorized = ind2rgb(im2uint8(imageVar01), hot); % Colorize imageVar01 with the colormap
        imshow(imageVar01_colorized);
        hold on;
        % Overlay the second image (damage) with transparency
        damage = double(model.Damage_E) / 255;
        h = imshow(damage);
        set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
        % colormap(gray);
        title('Overlay with Transparency');
        hold off;
        %---
        subplot(232);
        % colormap(parula); 
        imageVar2 = sum(electrode_sensitivity_neuron_E_rec_sqr(:,:,opt.electrode_record), 3);
        surf(imageVar2');
        colorbar; % colormap(Turbo)
        shading interp;
        title('Position of electrodes');
        view(40,40);
        %---
        subplot(233);
        hold on;
        for k = 1:numel(opt.electrode_record)
            plot3(xMat(:,k), yMat(:,k), zMat(:,k),...
            'Linewidth', 0.5, ...
            'color', clrs(rem(k - 1, noc) + 1, :), ...
            'DisplayName', ['[' num2str(opt.electrode_record(k)) ']']);
        end
        grid; lgd = legend;
        title(lgd, "Electrode's index");
        xlabel('x'); ylabel('y'); zlabel('z');
        title('Voltage of selected electrodes');
        view(40,40); 

        %----
        subplot(2,3,[4 6]);
        hold on;
        plotSecLines(model, params.step, [-1600, 200]);
        % plot(Time, electrode_neuron_E_rec(1,tline), 'color', [0.4940 0.1840 0.5560]);
        miny = 0; maxy = -20000;
        for k = 1:size(opt.electrode_record, 2)
            % neurInd = sub2ind([params.mneuro_E, params.nneuro_E], numofNeuro(k, 1), numofNeuro(k, 2));
            % siGnal = electrode_neuron_E_rec(opt.electrode_record(k),tline);
            siGnal = electrode_neuron_total_rec(opt.electrode_record(k),tline);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 0.5, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[' num2str(opt.electrode_record(k)) ']']);
        end
        lgd = legend;
        title(lgd, "Electrode's index");
        title('Voltage of selected electrodes');
        ylabel('Mean Voltage');
        xlabel('Time (s)');
        ylim([miny, maxy]);
        xlim(vectime);
        disp('Finished!');
    end

    %% ===================== NETWORK ACTIVITY MAPS =====================
    % These plots show spatial-temporal activity patterns across the entire network

    %% Part 3 - Voltage Raster Plot
    % Shows membrane potential of all neurons over time as a raster plot
    if plotSpikeRaster ~= 0
        plotSpikeRasterWiththreshold = 0;
        if plotSpikeRasterWiththreshold
            disp('Started - plotSpikeRaster');
            disp('Excitatory plot');
            % Excitatory plot
            ssi = setsizee(plotSpikeRaster);
            BBC = model.V_line_E;
            BBC = BBC > -40;
            BBC = BBC(:, tline);
            
            figure(); ssi();
            hold on;
            imagesc(BBC);
            colormap([1 1 1; 0 0 0]);  % Just black and white
            set(gca, 'Color', 'white');
            colorbar;
            time_steps = double(model.T_Iapp);
            time_steps = time_steps(:);
            scaled_time_steps = time_steps*params.step;
            for i = 1:length(time_steps)
                xline(time_steps(i), '--r','LineWidth', 1);
            end
            im = image(damage);
            im.AlphaData = backOpacity;
            xlabel('Time (s)');
            ylabel('Excitatory Neurons');
            title('Voltage map - Excitatory network');
            xlim(vectimen);
            ylim([0.5, params.quantity_neurons_E + 0.5]);
            
            % Inhibitory plot
            disp('Inhibitory plot');
            ssi = setsizee(plotSpikeRaster);
            BBC = model.V_line_I;
            BBC = BBC > -40;
            figure(); ssi(); 
            hold on;
            imagesc(BBC);
            colormap([1 1 1; 0 0 0]);  % Just black and white
            set(gca, 'Color', 'white');
            colorbar;
            time_steps = double(model.T_Iapp);
            time_steps = time_steps(:);
            scaled_time_steps = time_steps*params.step;
            for i = 1:length(time_steps)
                xline(time_steps(i), '--r','LineWidth', 1);
            end
            im = image(damage_I);
            im.AlphaData = backOpacity;
            plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
            xlabel('Time (s)');
            ylabel('Inhibitory Neurons');
            title('Voltage map - Inhibitory network');
            xlim(vectimen);
            ylim([0.5, params.quantity_neurons_I + 0.5]);

            disp('Finished!');
        else
            disp('Started - plotSpikeRaster');
            disp('Excitatory plot');
            % Excitatory plot
            ssi = setsizee(plotSpikeRaster);
            BBC = model.V_line_E;
            BBC = BBC(:, tline);
            figure(); ssi();
            % plot(0, 0); 
            hold on;
            colrs = opt.colrs3;
            locs = opt.locs3;
            cmap = makeGradient(256, colrs, locs, 0);
            imagesc(BBC);
            colormap(cmap);
            if opt.limits3
                caxis(opt.limits3);
            end
            colorbar;
            plotSecLines(model, 1, [0, size(BBC, 1) + 1]);

            im = image(damage);
            im.AlphaData = backOpacity;
            xlabel('Time (s)');
            ylabel('Excitatory Neurons');
            title('Voltage map - Excitatory network');
            xlim(vectimen);
            ylim([0.5, params.quantity_neurons_E + 0.5]);

            % Inhibitory plot
            disp('Inhibitory plot');
            ssi = setsizee(plotSpikeRaster);
            BBC = model.V_line_I;
            figure(); ssi();
            % plot(0, 0); 
            hold on;
            
            colrs = opt.colrs3;
            locs = opt.locs3;
            cmap = makeGradient(256, colrs, locs, 0);
            
            imagesc(BBC);

            colormap(cmap);
            colorbar;
            if opt.limits3
                caxis(opt.limits3);
            end

            im = image(damage_I);
            im.AlphaData = backOpacity;
            plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
            xlabel('Time (s)');
            ylabel('Inhibitory Neurons');
            title('Voltage map - Inhibitory network');
            xlim(vectimen);
            ylim([0.5, params.quantity_neurons_I + 0.5]);

            disp('Finished!');
        end
        
    end

    %% Part 4 - Glutamate Maps
    % Shows glutamate levels in neurons and astrocyte zones over time
    if plotGluMAP ~= 0
        if opt.plotGluMAPNeuro
            % from neurons
            disp('Started - plotGluMAP');
            ssi = setsizee(plotGluMAP);
            BBC = model.G;
            BBC = BBC(:, tline);
            figure(); ssi();
            plot(0, 0); hold on;
            colrs = opt.colrs4;
            locs = opt.locs4;
            cmap = makeGradient(256, colrs, locs, 0);
            imagesc(BBC);
            colormap(cmap);
            colorbar;

            % if opt.limits4
            %     caxis(opt.limits4);
            % end
            % im = image(damage);
            % im.AlphaData = backOpacity;

            plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
            xlabel('Time (s)');
            ylabel('Neurons');
            title('Glutamate Maps (Neurons)');
            xlim(vectimen);
            ylim([0.5, params.quantity_neurons_E + 0.5]);
            disp('Finished!');
        end
        %-------------------------
        disp('Started - plotGluAZMAP');
        % from neurons
        ssi = setsizee(plotGluMAP);
        BBC = model.neuron_astrozone_activity_Glu;
        BBC = reshape(BBC, params.mastro * params.nastro, size(BBC,3));
        BBC = BBC(:, tline);
        figure(); ssi();
        plot(0, 0); hold on;
        colrs = opt.colrs44;
        locs = opt.locs44;
        cmap = makeGradient(256, colrs, locs, 0);
        imagesc(BBC);
        colormap(cmap);
        colorbar;
        plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
        xlabel('Time (s)');
        ylabel('Astrocytes');
        title('Glu astrozone (astrocytes)');
        xlim(vectimen);
        ylim([0.5, params.mastro * params.nastro + 0.5]);
        disp('Finished!');
    end

    %% Part 5 - Applied Current Map
    % Shows stimulation pattern across all neurons over time
    if plotIappMAP ~= 0
        disp('Started - plotIappMAP');
        ssi = setsizee(plotIappMAP);
        HH = zeros(params.quantity_neurons_E, params.n);
        k = 1;
        for j = 1:params.nneuro_E
            for i = 1:params.mneuro_E
                HH(k, :) = reshape(model.Iapp_v_full(i, j, :), [1, params.n]);
                k = k + 1;
            end
        end
        BBC = HH(:, tline);
        figure(); ssi();
        plot(0, 0); hold on;
        colrs = opt.colrs5;
        locs = opt.locs5;
        cmap = makeGradient(10, colrs, locs, 0);
        imagesc(BBC);
        colormap(cmap);
        colorbar;
        if opt.limits5
            caxis(opt.limits5);
        end
        im = image(damage);
        im.AlphaData = backOpacity;
        plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
        xlabel('Time (s)');
        ylabel('Neurons');
        title('I Applied Map');
        xlim(vectimen);
        ylim([0.5, params.quantity_neurons_E + 0.5]);
        disp('Finished!');
    end

    %% Part 6 - Calcium Map
    % Shows calcium concentration across all astrocytes over time
    if plotCalMAP ~= 0
        disp('Started - plotCalMAP');
        ssi = setsizee(plotCalMAP);
        HH = zeros(params.mastro * params.nastro, params.n);
        k = 1;
        for j = 1:params.nastro
            for i = 1:params.mastro
                HH(k, :) = reshape(model.Ca(i, j, :), [1, params.n]);
                k = k + 1;
            end
        end
        BBC = HH(:, tline);
        figure(); ssi();
        % subplot(211);
        colrs = opt.colrs6;
        locs = opt.locs6;
        cmap = makeGradient(256, colrs, locs, 0);
        % cmap = jet(256);
        imagesc(BBC);
        colormap(cmap);
        colorbar;

        if sum(opt.limits6)
            caxis(opt.limits6);
        end

        hold on;
        plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
        % plotSecLines(model, 1, [0, size(BBC, 1)]);
        ylabel('Neurons');
        title('Calcium');
        xlim(vectimen);
        % subplot(212);
        % plot(Time, sum(BBC) / size(BBC, 1));
        % plot(Time(1:23999), sum(BBC) / size(BBC, 1));
        % ylabel('Mean of Ca^{2+} (\mu{}M)');
        % xlabel('Time (s)');
        % xlim(vectime);
        disp('Finished!');
    end

    %% Part 7 - Single Astrocyte Analysis: Calcium and IP3 Dynamics
    % Shows calcium signaling and IP3 concentration for selected astrocytes
    if plotCalIP3 ~= 0
        disp('Started - plotCalIP3');
        ssi = setsizee(plotCalIP3);

        figure(); ssi();
        subplot(211);
        hold on;
        plotSecAreas(model, params.step, [0, 10]);
        miny = 0; maxy = 0;

        for k = 1:size(numofAstro, 1)
            siGnal = reshape(model.Ca(numofAstro(k, 1), numofAstro(k, 2), :), [1, params.n]);
            % siGnal = reshape(model.Gli_global(numofAstro(k, 1), numofAstro(k, 2), :), [1, params.n]);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 1, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofAstro(k, 1)), ', ', num2str(numofAstro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        lgd = legend;
        title(lgd, 'Astrocytes(j,k)');
        ylabel('Calcium');
        ylim([0.8 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
        subplot(212);
        hold on;
        plotSecAreas(model, params.step, [0, 30]);
        miny = 0; maxy = 0;

        for k = 1:size(numofAstro, 1)
            siGnal = reshape(model.IP3(numofAstro(k, 1), numofAstro(k, 2), :), [1, params.n]);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 1, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofAstro(k, 1)), ', ', num2str(numofAstro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        ylabel('IP3');
        ylim([0.8 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
        xlabel('Time (s)');

        disp('Finished!');
    end

    %% Part 9 - Single Neuron Current Analysis
    % Shows various current components for selected neurons
    if plotIappIsumIsyn ~= 0
        disp('Started - plotIappIsumIsyn');
        ssi = setsizee(plotIappIsumIsyn);

        figure(); ssi();
        subplot(211);
        hold on;
        plotSecAreas(model, params.step, [-1000, 1000]);
        miny = 0; maxy = 0;

        for k = 1:size(numofNeuro, 1)
            siGnal = reshape(model.Iapp_v_full(numofNeuro(k, 1), numofNeuro(k, 2), :), [1, params.n]);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal));
            % len = min(length(Time), length(siGnal(tline)));
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 1, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofNeuro(k, 1)), ', ', num2str(numofNeuro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        lgd = legend;
        title(lgd, 'Neurons(m,n)');
        ylabel('I applied');
        ylim([0.8 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
        subplot(212);
        hold on;
        plotSecAreas(model, params.step, [-1000, 1000]);
        miny = 0; maxy = 0;

        for k = 1:size(numofNeuro, 1)
            neurInd = sub2ind([params.mneuro_E, params.nneuro_E], numofNeuro(k, 1), numofNeuro(k, 2));
            siGnal = model.Isyn_line_EE(neurInd, :);
            miny = min([siGnal, miny]); maxy = max([siGnal, maxy]);
            len = min(length(Time), length(siGnal));
            % len = min(length(Time), length(siGnal(tline))); % you have to save "model.Isyn_line_EE" for params.n in the init_model
            plot(Time(1:len), siGnal(tline(1:len)), linestyle{rem(k - 1, nols) + 1}, ...
                'Linewidth', 1, ...
                'color', clrs(rem(k - 1, noc) + 1, :), ...
                'DisplayName', ['[', num2str(numofNeuro(k, 1)), ', ', num2str(numofNeuro(k, 2)), ']']);
        end

        ax = gca; ax.YGrid = 'on';
        ylabel('I syn');
        % ylim([1.1 * miny, 1.1 * maxy + eps]);
        xlim(vectime);
       
        disp('Finished!');
    end

    %% Part 12 - Working Memory Current Map (I_WM)
    % Shows working memory current across astrocyte network
    if plotG_WM ~= 0
        disp('Started - plotG_WM');
        ssi = setsizee(plotG_WM);

        HH = zeros(params.mastro * params.nastro, params.n);
        k = 1;

        for j = 1:params.nastro

            for i = 1:params.mastro
                % HH(k, :) = reshape(params.Eta_WM * model.Gli_global(i, j, 1:params.n), [1, params.n]);
                HH(k, :) = reshape(model.I_WM(i, j, 1:params.n), [1, params.n]);
                k = k + 1;
            end

        end

        BBC = HH(:, tline);
        figure(); ssi();
        plot(0, 0); hold on;
        colrs = opt.colrs12;
        locs = opt.locs12;
        cmap = makeGradient(2, colrs, locs, 0);
        imagesc(BBC);
        colormap(cmap);
        colorbar;

        if opt.limits12
            caxis(opt.limits12);
        end

        plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
        xlabel('Time (s)');
        ylabel('Astrocytes');
        title('I WM');
        xlim(vectimen);
        ylim([0.5, params.mastro * params.nastro + 0.5]);
        disp('Finished!');
    end

    %% Part 14 - Neuron-to-Astrocyte Current Map (I_neuro_astro)
    % Shows coupling current from neurons to astrocytes
    if plotIneuroAstro ~= 0
        disp('Started - plotIneuroAstro');
        ssi = setsizee(plotIneuroAstro);

        % HH = zeros(params.mastro * params.nastro, params.n);
        % k = 1;
        % for j = 1:params.nastro
        %     for i = 1:params.mastro
        %         HH(k, :) = reshape(model.Ineuro(i, j, :), [1, params.n]);
        %         k = k + 1;
        %     end
        % end
        HH = reshape(model.Ineuro(:,:,1:params.n), [params.mastro * params.nastro, params.n]);
        BBC = HH(:, tline);
        figure(); ssi();
        plot(0, 0); hold on;
        colrs = opt.colrs14;
        locs = opt.locs14;
        cmap = makeGradient(2, colrs, locs, 0);
        % cmap = jet(256);
        imagesc(BBC);
        colormap(cmap);
        colorbar;
        if opt.limits14
            caxis(opt.limits14);
        end
        plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
        xlabel('Time (s)');
        ylabel('Astrocytes');
        title('I Neuro Astro');
        xlim(vectimen);
        ylim([0.5, params.mastro * params.nastro + 0.5]);
        disp('Finished!');
    end

    %% Part 17 - Gliotransmitter Global Map
    % Shows whole-cell gliotransmitter activity across astrocyte network
    if plotGli_global ~= 0
        disp('Started - plotGli_global');
        ssi = setsizee(plotGli_global);

        HH = zeros(params.mastro * params.nastro, params.n);
        k = 1;
        for j = 1:params.nastro
            for i = 1:params.mastro
                HH(k, :) = reshape(model.Gli_global(i, j, 1:params.n), [1, params.n]);
                k = k + 1;
            end
        end
        BBC = HH(:, tline);
        figure(); ssi();
        plot(0, 0); hold on;
        colrs = opt.Gli_global_clr;
        locs = opt.Gli_global_loc;
        cmap = makeGradient(256, colrs, locs, 0);
        imagesc(BBC);
        colormap(cmap);
        colorbar;
        if opt.limits12
            caxis(opt.limits12);
        end
        plotSecLines(model, 1, [0, size(BBC, 1) + 1]);
        xlabel('Time (s)');
        ylabel('Astrocytes');
        title('Gli global (whole-cell)');
        xlim(vectimen);
        ylim([0.5, params.mastro * params.nastro + 0.5]);
        disp('Finished!');
    end

    %% ===================== ELECTRODE ANALYSIS & COMPARISONS =====================
    % These plots show electrode positions, damage overlays, and comparative analyses

    %% Part 18 - Electrodes vs Damage Overlay
    % Shows electrode positions with brain damage overlay
    if plotElectrodesVSdamage ~= 0
        disp('Started - plotElectrodesVSdamage');
        ssi = setsizee(plotElectrodesVSdamage);
        figure(); ssi();
        sgtitle('Electrodes with damage overlay');
        subplot(121);
        electrode_sensitivity_neuron_E_stim_sqr = reshape(model.electrode_sensitivity_neuron_E_stim, params.mneuro_E, params.nneuro_E, []);
        imageVar_stim = sum(electrode_sensitivity_neuron_E_stim_sqr(:,:,opt.electrode_stim), 3);
        imageVar_stim = imageVar_stim.';
        imageVar_stim_inv = 1 - imageVar_stim;
        colormap(hot);
        colrs = [255, 255, 255; 29, 67, 80; 237, 31, 30] / 255;
        locs = [0, 0.5, 1];
        cmap = makeGradient(256, colrs, locs, 0);
        imageVar_stim_inv_colorized = ind2rgb(im2uint8(imageVar_stim), cmap); % Colorize imageVar_stim_inv with the colormap
        imshow(imageVar_stim_inv_colorized);
        hold on;
        % Overlay the second image (damage) with transparency
        damage = double(model.Damage_E) / 255;
        h = imshow(damage);
        set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
        colormap(gray);
        title('Stimulating Electrodes');
        hold off;

        subplot(122);
        electrode_sensitivity_neuron_E_rec_sqr = reshape(model.electrode_sensitivity_neuron_E_rec, params.mneuro_E, params.nneuro_E, []);
        imageVar1 = sum(electrode_sensitivity_neuron_E_rec_sqr(:,:,opt.electrode_record), 3);
        imageVar01 = 1 - imageVar1;
        colormap(hot); % Apply the colormap to the current figure
        imageVar01_colorized = ind2rgb(im2uint8(imageVar1), cmap); % Colorize imageVar01 with the colormap
        imshow(imageVar01_colorized);
        hold on;
        % Overlay the second image (damage) with transparency
        damage = double(model.Damage_E) / 255;
        h = imshow(damage);
        set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
        colormap(gray);
        title('Recording Electrodes');
        hold off;
        disp('Finished!');
    end

    %% Part 19 - Electrodes vs Input Overlay
    % Shows electrode positions with input pattern overlays
    if plotElectrodesVSinputs ~= 0
        disp('Started - plotElectrodesVSinputs');
        ssi = setsizee(plotElectrodesVSinputs);
        NumimG = numel(params.learn_order);
        for imG = 1:NumimG
            imGindex = params.learn_order(imG);
            figure(); ssi();
            sgtitle(['Electrodes with input (',num2str(imGindex-1),') overlay']);
            subplot(121);
            electrode_sensitivity_neuron_E_stim_sqr = reshape(model.electrode_sensitivity_neuron_E_stim, params.mneuro_E, params.nneuro_E, []);
            imageVar_stim = sum(electrode_sensitivity_neuron_E_stim_sqr(:,:,opt.electrode_stim), 3);
            imageVar_stim_inv = 1 - imageVar_stim;
            colormap(hot); 
            imageVar_stim_inv_colorized = ind2rgb(im2uint8(imageVar_stim_inv), hot); % Colorize imageVar_stim_inv with the colormap
            imshow(imageVar_stim_inv_colorized);
            hold on;
            % Overlay the second image (inputs) with transparency
            inputs = double(cell2mat(model.images(imGindex))) / 255;
            h = imshow(inputs);
            set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
            colormap(gray); 
            title('Stimulating Electrodes');
            hold off;

            subplot(122);
            electrode_sensitivity_neuron_E_rec_sqr = reshape(model.electrode_sensitivity_neuron_E_rec, params.mneuro_E, params.nneuro_E, []);
            imageVar1 = sum(electrode_sensitivity_neuron_E_rec_sqr(:,:,opt.electrode_record), 3);
            imageVar01 = 1 - imageVar1;
            colormap(hot); 
            imageVar01_colorized = ind2rgb(im2uint8(imageVar01), hot); % Colorize imageVar01 with the colormap
            imshow(imageVar01_colorized);
            hold on;
            % Overlay the second image (inputs) with transparency
            inputs = double(cell2mat(model.images(imGindex))) / 255;
            h = imshow(inputs);
            set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
            colormap(gray);
            title('Recording Electrodes');
            hold off;
        end
        disp('Finished!');
    end

    %% Part 20 - Electrode Recording Integral Analysis
    % Integrates electrode signals over time windows with filtering
    if plotElectrodeRecord_integral ~= 0 
        disp('Started - plotElectrodeRecord_integral');
        ssi = setsizee(plotElectrodeRecord_integral);
        %-----
        if params.stimulation_mode %automatic stimulation
            electrode_neuron_E_rec = model.electrode_neuron_E_rec;
            electrode_neuron_I_rec = model.electrode_neuron_I_rec;
            electrode_neuron_total_rec ...
                = electrode_neuron_E_rec ...
                + electrode_neuron_I_rec ...
                + model.adjusted_matrix(:, 1);
            
        else %manual stimulation
            electrode_neuron_E_rec = zeros(params.n_electrode_record, params.n);
            electrode_neuron_I_rec = zeros(params.n_electrode_record, params.n);
            electrode_neuron_total_rec = zeros(params.n_electrode_record, params.n);
            for n = 1 : params.n_electrode_record
                for i = 1:params.n
                electrode_neuron_E_rec(n,i) = sum(model.electrode_sensitivity_neuron_E_rec(:,n).*model.V_line_E(:,i));
                electrode_neuron_I_rec(n,i) = sum(model.electrode_sensitivity_neuron_I_rec(:,n).*model.V_line_I(:,i));
                electrode_neuron_total_rec(n,i) = electrode_neuron_E_rec(n,i) + electrode_neuron_I_rec(n,i);
                end
            end
            %normalizing the values
            electrode_neuron_total_rec = adjust_initial_values(electrode_neuron_total_rec);
        end

        % train_patterns = size (params.learn_order,2);
        % for i = 1:train_patterns 
        %     integral_interval(i) = [model.T(i), model.T(i)];
        %     integral_electrode_rec(:,:,i) = calculate_integral(params, electrode_neuron_total_rec, integral_interval(i), [1 : params.n_electrode_record]);
        % end


        % Calculate integral and apply filter
        integral_electrode_rec = calculate_integral(params, electrode_neuron_total_rec, opt.integralTime, [1 : params.n_electrode_record]);
        integral_electrode_rec = reshape(integral_electrode_rec, [params.electrode_grid_record(1) params.electrode_grid_record(2)]);
        integral_electrode_rec = integral_electrode_rec';
        % integral_electrode_rec_filter1 = integral_electrode_rec > opt.integralThreshold1;
        top_values = (max(integral_electrode_rec(:))*50)/100;
        integral_electrode_rec_filter1 = integral_electrode_rec > top_values;

        % first interval integral
        figure('Position', [50, 250, 600, 500]);
        imagesc(integral_electrode_rec);
        colorbar; axis square; axis off;
        [rows, cols] = size(integral_electrode_rec); 
        for i = 1:rows
            for j = 1:cols
                text(j, i, sprintf('%.1f', integral_electrode_rec(i,j)), ... 
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 14);
            end
        end
        titleStr = ['Integral from ', num2str(opt.integralTime(1)), ' to ', num2str(opt.integralTime(2))];
        title(titleStr, 'FontSize', 16);

        % filter 1
        figure('Position', [500, 250, 600, 500]);
        imagesc(integral_electrode_rec_filter1);
        [rows, cols] = size(integral_electrode_rec_filter1); 
        for i = 1:rows
            for j = 1:cols
                text(j, i, sprintf('%d', integral_electrode_rec_filter1(i,j)), ... 
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 14);
            end
        end
        titleStr = ['Filter 1 from ', num2str(opt.integralTime(1)), ' to ', num2str(opt.integralTime(2))];
        title(titleStr, 'FontSize', 16);

        % comparing 1
        figure('Position', [850, 250, 600, 500]);
        imagesc(integral_electrode_rec);
        [rows, cols] = size(integral_electrode_rec_filter1);
        for i = 1:rows
            for j = 1:cols
                value = integral_electrode_rec_filter1(i,j);
                if value == 1
                    color = 'w';
                    weight = 'bold';
                else
                    color = 'black';
                    weight = 'normal';
                end
                text(j, i, sprintf('%d', value), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 14, ...
                    'Color', color, 'FontWeight', weight);
            end
        end
        titleStr = ['Integral 1 with Filter 1 numbers [', num2str(opt.integralTime(1)), ' to ', num2str(opt.integralTime(2)),']'];
        title(titleStr, 'FontSize', 16);



        % comparing 2
        % Create a new figure with a larger size
        figure('Position', [750, 150, 700, 600]);

        % Create a diverging colormap centered around the threshold
        cmap = [flipud(winter(128)); autumn(128)];
        cmin = min(integral_electrode_rec(:));
        cmax = max(integral_electrode_rec(:));
        crange = max(abs(cmin-opt.integralThreshold1), abs(cmax-opt.integralThreshold1));
        clim = [opt.integralThreshold1-crange, opt.integralThreshold1+crange];

        % Plot the matrix
        imagesc(integral_electrode_rec, clim);
        colormap(cmap);
        c = colorbar;
        c.Label.String = 'Integral Value';
        axis square; axis off;

        % Add contour for filtered areas with improved visibility
        hold on;
        [B, L] = bwboundaries(integral_electrode_rec_filter1, 'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'Color', [0.5, 0.5, 0.5, 0.7], 'LineWidth', 9);
        end

        % Add text labels in each cell
        [rows, cols] = size(integral_electrode_rec);
        for i = 1:rows
            for j = 1:cols
                value = integral_electrode_rec(i,j);
                if integral_electrode_rec_filter1(i,j)
                    color = 'k';  % White text for values above threshold
                    weight = 'bold';
                else
                    color = [201, 201, 201]/255; % 0, 239/255, 136/255];  % Black text for values below threshold
                    weight = 'normal';
                end
                text(j, i, sprintf('%.1f', value), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 14, 'Color', color, 'FontWeight', weight);
            end
        end

        titleStr = sprintf('Integral from %.2f to %.2f (Threshold: %.2f)', ...
            opt.integralTime(1), opt.integralTime(2), opt.integralThreshold1);
        title(titleStr, 'FontSize', 16);


        integral_electrode_rec_t2 = calculate_integral(params, electrode_neuron_total_rec, opt.integralTime2, [1 : params.n_electrode_record]);
        integral_electrode_rec_t2 = reshape(integral_electrode_rec_t2, [params.electrode_grid_record(1) params.electrode_grid_record(2)]);
        integral_electrode_rec_t2 = integral_electrode_rec_t2';
        integral_electrode_rec_t2 = integral_electrode_rec_t2.*integral_electrode_rec_filter1;
        integral_electrode_rec_t2_filter = (integral_electrode_rec_t2 > 0) & (integral_electrode_rec_t2 < opt.integralThreshold2);
        
        figure('Position', [100, 50, 600, 500]);
        imagesc(integral_electrode_rec_t2);
        [rows, cols] = size(integral_electrode_rec_t2); 
        for i = 1:rows
            for j = 1:cols
                text(j, i, sprintf('%.1f', integral_electrode_rec_t2(i,j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 14);
            end
        end
        titleStr = ['Integral from ', num2str(opt.integralTime2(1)), ' to ', num2str(opt.integralTime2(2))];
        title(titleStr, 'FontSize', 16);
        
        
        figure('Position', [700, 50, 600, 500]);
        imagesc(integral_electrode_rec_t2_filter);
        [rows, cols] = size(integral_electrode_rec_t2_filter);
        for i = 1:rows
            for j = 1:cols
                text(j, i, sprintf('%d', integral_electrode_rec_t2_filter(i,j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 14);
            end
        end
        titleStr = ['Filter 2 on the Integral from ', num2str(opt.integralTime2(1)), ' to ', num2str(opt.integralTime2(2))];
        title(titleStr, 'FontSize', 16);
    end
        
    %% Part 21 - Damage vs Input Comparison
    % Shows relationship between damage patterns and input patterns
    if plot_Damage_VS_Inputs ~= 0
        disp('Started - plot_Damage_VS_Inputs');
        ssi = setsizee(plot_Damage_VS_Inputs);

        NumimG = numel(params.learn_order);
        
        damage = double(model.Damage_E) / 255;
        damage = 1 - damage;
        cmap = makeGradient(256, opt.DamVSinpu_clr,opt.DamVSinpu_loc, 1);
        damage_colorized = ind2rgb(im2uint8(damage), cmap); % Colorize imageVar_stim_inv with the colormap
        

        % for imG = 1:NumimG
        %     imGindex = params.learn_order(imG);
        %     figure(); ssi();

        %     sgtitle(['input (',num2str(imG),')']);
        %     imshow(damage_colorized);
        %     hold on;
            
        %     % Overlay the second image (inputs) with transparency
        %     inputs = double(cell2mat(model.images(imGindex))) / 255;
        %     h = imshow(inputs);
        %     set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
        %     colormap(gray); 
        %     % colorbar;
            
        %     hold off;
        %     set(gcf, 'WindowState', 'maximized');
        % end

        figure();
        for imG = 1:NumimG
            imGindex = params.learn_order(imG);
            subplot(ceil(sqrt(NumimG)), ceil(NumimG / sqrt(NumimG)), imG);
            
            ssi();
            imshow(damage_colorized);
            hold on;
            
            % Overlay the second image (inputs) with transparency
            inputs = double(cell2mat(model.images(imGindex))) / 255;
            h = imshow(inputs);
            set(h, 'AlphaData', opt.Opacity); % Adjust transparency level as needed
            colormap(gray);
            
            hold off;
            title(['Input (', num2str(imG), ')'], 'FontSize', 10);
        end
        
        set(gcf, 'WindowState', 'maximized');
        disp('Finished!');
    end

    %% Part 22 - Electrode Accuracy/Similarity Analysis
    % Shows pattern recognition accuracy over time
    if plotAccuracy_elec ~= 0
        disp('Started - plotAccuracy_elec');
        if plotAccuracy_elec == 11
            patTers = params.test_order;
            ssi = setsizee(2);
        else
            patTers = params.learn_order;
            ssi = setsizee(plotAccuracy_elec);
        end
        
        figure(); ssi();
        hold on;
        plotSecAreas(model, params.step, [-200, 200]);
        
        [similarity, tLi, target] = accuracy_elec(params, model, patTers);
        
        for k = 1:length(patTers)
            similarity2(:,k) = smooth(similarity(:,k), 10);
        end
        
        for k = 1:length(patTers)
            siGnal = similarity2(:,k);
            plot(tLi*params.step*20, siGnal);
            % plot(tLi*params.step*20, siGnal, linestyle{rem(k - 1, nols) + 1}, ...
            %     'Linewidth', 0.5, ...
            %     'color', clrs(rem(k - 1, noc) + 1, :), ...
                % 'DisplayName', num2str(k));
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
        
    %% ===================== 3D VISUALIZATIONS =====================
    % These plots provide 3D perspectives and time-series snapshots of network activity

    %% Part 23 - 3D Glutamate Astrocyte Zone Visualization
    % Shows glutamate activity in astrocyte zones at different time points
    if plotGluAstroZone3D ~= 0
        disp('Started - plotGluAZMAP');

        ssi = setsizee(plotGluAstroZone3D);
        BBC = model.neuron_astrozone_activity_Glu;
        colrs = opt.colrs44;
        locs = opt.locs44;
        cmap = makeGradient(256, colrs, locs, 1);
        % maxGlu = max(max(model.neuron_astrozone_activity_Glu(:)));

        timeStamps = model.T_Iapp([numel(params.learn_order)+1:end], 2);
        
        figure();
        ssi();
        for i = 1:length(timeStamps)
            subplot(2, 4, i); 
            imagesc(BBC(:,:,timeStamps(i)));
            colormap(cmap);
            axis tight; 
            set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
            % colorbar; 
            % caxis([0, maxGlu]);
            caxis([0, 8]);
            axis equal tight; 
            title(['Time: ', num2str(double(timeStamps(i)) / 10000, '%.2f')]);
        end

        disp('Finished!');

    end

    %% Part 24 - 3D Calcium Dynamics Visualization
    % Shows calcium concentration across astrocytes at different time points
    if plotCalcium3D ~= 0
        disp('Started - plotCalcium3D');

        ssi = setsizee(plotCalcium3D);
        BBC = model.Ca;
        
        % same as calcium map colors
            % colrs = opt.colrs6;
            % locs = opt.locs6;

        % with more spaces between colors for better observation of the changes
            colrs = [255, 255, 255; 29, 67, 80; 237,31,30; 104, 0, 104]/255;
            max_ca = max(max(model.Ca(:)));
            min_ca = min(min(model.Ca(:)));
            locs = [min_ca, max_ca*0.3,max_ca*0.7, max_ca];

        cmap = makeGradient(256, colrs, locs, 0);

        % maxCa = max(max(model.Ca(:)));

        timeStamps = model.T_Iapp([numel(params.learn_order)+1:end], 2);
        
        % timeStamps = timeStamps([5;6;7;8;11;12;13;14], 1);
        % timeStamps = timeStamps([1;2], 1);

        figure();
        ssi();
        for i = 1:length(timeStamps)
            subplot(2, 4, i); 
            imagesc(BBC(:,:,timeStamps(i)));
            colormap(cmap);
            axis tight; 
            set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
            % colorbar; 
            % caxis([0, maxCa]);
            caxis([params.ca_0, 0.8]); 
            axis equal tight; 
            title(['Time: ', num2str(double(timeStamps(i)) / 10000, '%.2f')]);
        end

        disp('Finished!');

    end


    %% Part 25 - 3D Neuron-to-Astrocyte Current Training Visualization
    % Shows training dynamics of neuron-astrocyte coupling
    if plotIneuronAstro_train_3D ~= 0
        disp('Started - plotIneuronAstro_train_3D');

        ssi = setsizee(plotIneuronAstro_train_3D);
        % BBC = model.Ineuro;
        BBC = model.neuron_astrozone_activity_Glu;
        colrs = opt.colrs14;
        locs = opt.locs14;
        cmap = [1 1 1; 0.6,0.1,0.2];
        timeStamps = model.T_Iapp([1:numel(params.learn_order)], 2);
        
        figure();
        ssi();
        for i = 1:length(timeStamps)
            subplot(2, 4, i); 
            imagesc(BBC(:,:,timeStamps(i))>3.4);
            colormap(cmap);
            axis tight; 
            set(gca, 'XTick', [], 'YTick', []); % Removes tick marks
            % colorbar; 
            % caxis([0, 8]);
            axis equal tight; 
            title(['Time: ', num2str(double(timeStamps(i)) / 10000, '%.2f')]);
        end

        disp('Finished!');

    end



    %% Part 26 - 3D Gliotransmitter Global Visualization
    % Shows whole-cell gliotransmitter activity at different time points
    if plotGli_global3D ~= 0
        disp('Started - plotGli_global3D');
        ssi = setsizee(plotGli_global3D);
        BBC = model.Gli_global;

        colrs = opt.Gli_global_clr;
        locs = opt.Gli_global_loc;
        cmap = makeGradient(256, colrs, locs, 0);

        timeStamps = model.T_Iapp([numel(params.learn_order)+1:end], 2) + params.impact_astro;

        figure();
        ssi();
        for i = 1:length(timeStamps)
            subplot(2, 4, i); 
            imagesc(BBC(:,:,timeStamps(i)));
            colormap(cmap);
            axis tight; 
            set(gca, 'XTick', [], 'YTick', []); % Removes tick marks
            % colorbar; 
            % caxis([0, 8]);
            axis equal tight; 
            title(['Time: ', num2str(double(timeStamps(i)) / 10000, '%.2f')]);
        end

    disp('Finished!'); 
    end
    ress = 1;
end

%% ===================== HELPER FUNCTIONS =====================
% These functions provide utility operations for plotting and data processing


%% Functions
function ssi = setsizee(mode)

    switch mode
        case 1
            ssi = @(t) set(gcf, 'Position', [100, 100, 800, 500]);
        case 2
            ssi = @(t) set(gcf, 'Position', [100, 100, 700, 250]);
        case 3
            ssi = @(t) set(gcf, 'Position', [100, 100, 800, 300]);
        case 4
            ssi = @(t) set(gcf, 'Position', [100, 100, 800, 200]);
        otherwise
            ssi = @(t) set(gcf, 'Position', [0, 0, 3000, 3000]);
    end

end

function plotSecLines(model, timescalse, ampl)
    h1 = double(model.T_Iapp(:));

    for i = 1:length(h1)
        timeee = [h1(i) * timescalse, h1(i) * timescalse + timescalse];
        h = plot(timeee, ampl, 'k:');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

end

function plotSecAreas(model, timescalse, ampl)
    h1 = double(model.T_Iapp);

    for i = 1:size(h1, 1)
        timeee = [h1(i, 1) * timescalse, h1(i, 1) * timescalse + timescalse, ...
                      h1(i, 2) * timescalse, h1(i, 2) * timescalse + timescalse];
        h = area(timeee, [ampl(1), ampl(2), ampl(2), ampl(1)], ...
            'FaceColor', [0.96, 0.96, 0.96], ...
            'EdgeColor', [0.86, 0.86, 0.86], ...
            'LineStyle', '-.', ...
            'LineWidth', 0.7, ...
            'basevalue', -1000);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

end

% function PPo = scaleMap(HH)
%     rsip = size(HH,1)/size(HH,2);
%     nor  = floor(0.1/rsip);
%     PP   = zeros(size(HH,1)*nor, size(HH,2));
%     k    = 1;
%     for i = 1:size(HH,1)
%         for j = 1:nor
%             PP(k,:) = HH(i,:);
%             k = k + 1;
%         end
%     end
%     PPo = PP;
% end

function out = makeImage(HH, MAPcolor1)
    out(:, :, 3) = HH * MAPcolor1(3); %G
    out(:, :, 2) = HH * MAPcolor1(2); %R
    out(:, :, 1) = HH * MAPcolor1(1); %B
end

function pic = damageBackground(dat, noTime)
    dat = dat(:);
    pic = repmat(dat, 1, noTime);
end


function integral_electrode_rec = calculate_integral(params, electrode_neuron_total_rec, range, dimension)
    sampletime = params.step;
    
    % Adjust range
    if range == 0
        range = [1, size(electrode_neuron_total_rec, 2)];
    else
        range = floor(range / sampletime);
        range = max(range, 1);
        range(2) = min(range(2), size(electrode_neuron_total_rec, 2));
    end
 
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

