function [params, paramss] = model_parameters(need_set, new_params, new_paramss)
    persistent params_p paramss_p
	if nargin < 1 || ~need_set
        params = params_p;
        paramss = paramss_p;
        return;
	end
    if nargin == 1
    params = struct;
    paramss = struct; 
    % Note: Use "paramss" (double 's') to display parameters in GUI.
    %       Be careful: this may cause type changes in that parameter.
    %% Simulation Settings
    % Network configuration and input size
    params.simPattern = 3;
        % (N:71x71, A:35x35, AZ:2x2, requires 32-64G RAM)
        % 3: Enhance mode
        % 4: Disruption mode
        % 5: Augmentation mode
    
    %% Synaptic Damage Modes
    % paramss.impairmode options:
    %   0: Healthy network                                 (Experiments 1, 2, 3)
    %   1: Random damage with fixed amplitude              (Not used)
    %   2: Random damage with random amplitude             (Not used)
    %   3: Concentrated damage from damaging pattern       (Experiment 1 - part 2)
    %   4: Single neuron test mode

    % Additional settings for damage modes 1 & 2 (Not used)
    params.scenario1a2.probability  = 1; % [0-1] Higher value = more neurons damaged
    params.scenario1a2.amplitude    = 0; % [0-1] Lower value = stronger damage
                                         % In paramss.impairmode
                                         % Mode 1: All selected neurons get this strength
                                         % Mode 2: Minimum damage strength parameter
    %% Stimulation Settings
    % paramss.stimulation_mode: 0 = manual, 1 = automatic (changes stimulation & recording methods)
    % paramss.stimulation_Augmentation: 0 = disabled, 1 = enabled (requires manual mode)
        
    % Example electrode locations for manual stimulation:
        % Column 1: [1, 5, 9, 13]
        % Column 2: [2, 6, 10, 14]
        % Column 3: [3, 7, 11, 15]
        % Column 4: [4, 8, 12, 16]
        % Columns 2&3: [2, 6, 10, 14, 3, 7, 11, 15]
        % All electrodes: [1:16]
        % Example grouping: params.electrode_neuron_stimulate_G = {[6,10,14], [3,4,11,12]} <- Current default
    
    %% Experiment and Timings
	switch params.simPattern
        case 3
            paramss.simulation_mode = "Enhance";
            paramss.t_end = 5.6; 
            
            % Learning phase parameters
            params.learn_start_time         = 0.05;
            params.learn_impulse_duration   = 0.22; 
            params.learn_impulse_shift      = 0.3; 
            params.learn_order              = [1, 2, 3, 4];
            
            % Testing phase parameters
            params.test_start_time          = 1.6;
            params.test_impulse_duration    = 0.2; 
            params.test_impulse_shift       = 0.4; 
            params.test_order               = [1, 2, 3, 4];
            params.delayTests               = 0.8;
            
            % Stimulation parameters
            params.stimulation_time_astro   = 1.2; 
            paramss.stimulation_time_neuro_E = 3.2;
            params.stimulation_duration_Ex  = 0.4; 
            params.stimulation_time_neuro_E_Augmentation = [0.12, 0.42, 0.72];
            params.stimulation_time_neuro_I = 1.28; 
            params.stimulation_duration_Inh = 0.2; 

            % Stimulation control
            paramss.stimulation_mode         = 1;  % 0: manual, 1: automatic
            paramss.stimulation_Augmentation = 0;
            paramss.stimulation_neuro_E      = 1;  % excitatory stimulation enabled
            paramss.stimulation_neuro_I      = 0;  % inhibitory stimulation disabled
            params.electrode_neuron_stimulate = [15, 8, 12, 16];
            params.electrode_neuron_stimulate_G = {[6, 10, 14], [3, 4, 11, 12]};

            % Damage configuration
            paramss.impairmode = 3;  % concentrated damage based on damaging pattern
            % paramss.impairmode = 0;  % healthy network
            
        case 4
            paramss.simulation_mode = "Disruption";
            paramss.t_end = 3.5; 
            
            % Learning phase parameters
            params.learn_start_time         = 0.05;
            params.learn_impulse_duration   = 0.22;
            params.learn_impulse_shift      = 0.4; 
            params.learn_order              = [9, 3, 11, 7];

            % Testing phase parameters
            params.test_start_time          = 1.8; 
            params.test_impulse_duration    = 0.2;
            params.test_impulse_shift       = 0.4; 
            params.test_order               = [9, 3, 11, 7];

            % Disruption configuration
            params.disruption_pattern = [3, 7];  % patterns to disrupt during training
            params.disruptBothTrainTest = 1;     % 1: disrupt both phases, 0: disrupt training only

            % Stimulation parameters
            params.stimulation_time_astro   = 1;
            paramss.stimulation_time_neuro_E = 0.8;
            params.stimulation_time_neuro_E_Augmentation = [0.12, 0.42, 0.72];
            params.stimulation_duration_Ex  = 0.2;
            params.stimulation_time_neuro_I = 1.28; 
            params.stimulation_duration_Inh = 0.2; 
            
            % Stimulation control
            paramss.stimulation_mode          = 1;  % automatic mode
            paramss.stimulation_Augmentation  = 0;
            paramss.stimulation_neuro_E       = 0;  % excitatory stimulation disabled
            paramss.stimulation_neuro_I       = 1;  % inhibitory stimulation enabled
            params.electrode_neuron_stimulate = [15, 8, 12, 16];
            params.electrode_neuron_stimulate_G = {[6, 10, 14], [3, 4, 11, 12]};

            % Damage configuration
            paramss.impairmode = 0;  % healthy network
            
        case 5
            paramss.simulation_mode = "Augmentation";
            paramss.t_end = 2.6; 

            % Learning phase parameters
            params.learn_start_time         = 0.05;
            params.learn_impulse_duration   = 0.22;
            params.learn_impulse_shift      = 0.3; 
            params.learn_order              = [4, 5];

            % Testing phase parameters
            params.test_start_time          = 1;  
            params.test_impulse_duration    = 0.2;
            params.test_impulse_shift       = 0.4; 
            params.test_order               = [4, 2, 5, 1]; 

            % Disruption configuration
            params.disruption_pattern = [3, 7];

            % Stimulation parameters
            params.stimulation_time_astro   = 1;
            paramss.stimulation_time_neuro_E = 0.8;
            params.stimulation_time_neuro_E_Augmentation = [0.12, 0.42];  % timings for augmentation
            params.stimulation_duration_Ex  = 0.2; 
            params.stimulation_time_neuro_I = 1.28;
            params.stimulation_duration_Inh = 0.2;

            % Stimulation control
            paramss.stimulation_mode          = 0;  % manual mode
            paramss.stimulation_Augmentation  = 1;  % augmentation enabled
            paramss.stimulation_neuro_E       = 1;  % excitatory stimulation enabled
            paramss.stimulation_neuro_I       = 0;  % inhibitory stimulation disabled
            params.electrode_neuron_stimulate = [15, 8, 12, 16];  % selected electrodes (not used)
            params.electrode_neuron_stimulate_G = {[6, 10, 14], [3, 4, 11, 12]};  % electrode groups (used for [0.12, 0.42] timings)
            
            % Damage configuration
            paramss.impairmode = 0;  % healthy network
	end 
    
    params.step  = 0.0001; 
    params.n     = fix(paramss.t_end / params.step);

    %% Applied Pattern Current
    params.variance_learn           = 0.02;  
    params.variance_test            = 0.05;  
    params.Iapp_learn               = 12;    % µA
    params.Iapp_test                = 8;     % µA

    %% Movie
    params.after_sample_frames      = 200;
    params.before_sample_frames     = 1;

    %% Poisson Noise
    params.poisson_nu                = 1.4;    % Hz
    params.poisson_n_impulses        = 7;      
    params.poisson_impulse_duration  = fix(0.02 / params.step);  
    params.poisson_impulse_initphase = fix(1.4 / params.step);   
    params.poisson_amplitude         = 5;      % µA  20

    %% Runge-Kutta steps
    params.u2                       = params.step / 2;
    params.u6                       = params.step / 6;

    %% Network size
    params.mneuro_E = 71; 
    params.nneuro_E = 71; 
    params.mneuro_I = 35;
    params.nneuro_I = 35;
    params.quantity_neurons_E = params.mneuro_E * params.nneuro_E;
    params.quantity_neurons_I = params.mneuro_I * params.nneuro_I;
    params.mastro   = 70;
    params.nastro   = 70;

    %% Initial conditions
    params.V_0          = -71.97;
    params.U_0          = -12.88;
    params.ca_0         = 0.0724;
    params.h_0          = 0.8863;
    params.ip3_0        = 0.8202;

    %% Neuron Model
    params.aa               = 0.1;	% Time scale of recovery variable
    params.b                = 0.2;	% Recovery variable coupling
    params.c                = -65;	% After-spike reset value (mV)
    params.d                = 2;    % After-spike recovery increment
    
    params.tau_Glu          = 0.1;  % Glutamate clearance time constant
	params.r_Glu            = 450;  % Glutamate release rate (µM/s)
    params.is_active        = 25;   % Activity threshold
    params.neuron_fired_thr = 30;	% Maximum presynaptic input current (µA)
    params.I_input_thr      = 25;	% Maximum neuron voltage output (mV)

    %% Synaptic connections
    % Exponential distribution parameters (average)
    params.N_connections = 40;  % Synapses per neuron
    params.lambda        = 5;   % Average exponential distribution

    params.lambda_EE = 4;
    params.lambda_ExI = 2;
    params.lambda_InE = 3;
    
    % Synapse counts for 71x71 network
    params.N_connections_EE = 40;      % E->E connections per neuron
    params.quantity_connections_EE = params.quantity_neurons_E * params.N_connections_EE;
    params.N_connections_ExI = 10;     % E->I connections per E neuron
    params.quantity_connections_ExI = params.quantity_neurons_E * params.N_connections_ExI;
    params.N_connections_InE = 30;     % I->E connections per I neuron
    params.quantity_connections_InE = params.quantity_neurons_I * params.N_connections_InE;


    %% Synapse Parameters
    % Excitatory/Inhibitory channel properties
    params.Esyn                = 0;      % Excitatory reversal potential (mV)
    params.Esyn_               = -90;    % Inhibitory reversal potential (mV)
    params.Ssyn                = 0.2;    % Slope of synaptic activation function
    
    % Working Memory
    params.Eta_WM              = 25;     % Astrocyte effect parameter
    
    % Initial synaptic weights
    params.W_EE_0 = 0.015;  % E->E synapses
    params.W_InE_0 = 0.15;  % I->E synapses
    params.W_ExI = 0.1;     % E->I synapses
    
    %% Neuron-Astrocyte Interaction
    params.t_neuro_astro            = 0.06; % Duration of neuron-astrocyte current (s)
    params.amplitude_neuro_astro    = 5;    % Amplitude of neuron-astrocyte current (µA)
    
    %% Astrocyte Model
    params.dCa                      = 0.05; % Calcium dynamics parameter
    params.dIP3                     = 0.05; % IP3 dynamics parameter (alt: 0.1)
	params.zeta_WM					= 0.10; % whole-cell activity consuption factor for WM (alt: 0.25)
    
    % Astrocyte observation window parameters
    window_astro_watch              = 0.01;   % Observation window duration (s)
    shift_window_astro_watch        = 0.001;  % Window shift interval (s)
    
    % Astrocyte impact parameters (working memory - whole-cell activity - I_astro_neuron)
    impact_astro                    = 0.1;    % Impact duration (s)
    params.impact_astro             = fix(impact_astro / params.step);
    params.window_astro_watch       = fix(window_astro_watch / params.step);
    params.shift_window_astro_watch = fix(shift_window_astro_watch / params.step);
    
    % Astrocyte input window
	input_window                    = 0.35;   % Delay between astrocyte inputs (s)
    params.input_window             = fix(input_window / params.step);
    
    % Glutamate release from neurons and astrocytes
    params.tau_Gli                 = 0.1;    % Glutamate time constant
    params.r_Gli_global            = 15;     % Glutamate release rate

    %% Working Memory Conditions and Electrodes
    params.Glu_memorize             = 3.4;  % Glutamate threshold for memory encoding
    params.Glu_recall_global        = 1.6;  % Glutamate threshold for memory recall
    params.ca_threshold_global      = 0.15; % Calcium threshold required for whole-cell activity
    params.electrode_lambda_record  = 2.5;  % Recording electrode sensitivity
    params.electrode_lambda_stim    = 5.5;  % Stimulation electrode sensitivity

    %% Memory Performance
    params.max_spikes_thr          = 30;      % Maximum spike threshold for memory evaluation (compute_memory_performance.m)

    %% Electrodes, Recording and Stimulation
    % Stimulation electrode grid configuration
    params.electrode_grid_stim     = [4, 4];
    params.n_electrode_stim        = prod(params.electrode_grid_stim);
    
    % Recording electrode grid configuration
    params.electrode_grid_record   = [8, 8];
    params.n_electrode_record      = prod(params.electrode_grid_record);
    
    % Recording and stimulation parameters
    params.stim_neuro_amp_E        = 20;      % Excitatory neuron stimulation amplitude (µA)
    params.stim_neuro_amp_I        = 50;      % Inhibitory neuron stimulation amplitude (µA)
    
    % Inhibitory neuron stimulation timing (Disrupt)
    params.stim_start_neuro_I      = fix(params.stimulation_time_neuro_I / params.step);
    params.stim_end_neuro_I        = params.stim_start_neuro_I + fix(params.stimulation_duration_Inh / params.step);
    params.stim_neuro_I_auto_dur   = params.stimulation_duration_Inh / params.step;
    
    % Excitatory neuron stimulation timing (Enhance)
    params.stim_start_neuro_E      = fix(paramss.stimulation_time_neuro_E / params.step);
    params.stim_end_neuro_E        = params.stim_start_neuro_E + fix(params.stimulation_duration_Ex / params.step);

    % Excitatory neuron stimulation timing (augmentation groups)
    params.stim_start_neuro_E_G    = fix(params.stimulation_time_neuro_E_Augmentation / params.step);
    params.stim_end_neuro_E_G      = params.stim_start_neuro_E_G + fix(params.stimulation_duration_Ex / params.step);
    
    % Data sampling rate for online stimulation
    params.dis_sampling            = 0.01 / params.step;
    
    %% 
    else
    params = new_params;
    paramss = new_paramss;
    end
    params_p = params;
    paramss_p = paramss;
end

