function model = step_neurons(model, params, is_STDP, i)
    %% Inputs
    V_E = model.V_line_E(:, i);
    U_E = model.U_line_E(:, 1);
    V_I = model.V_line_I(:, i);
    U_I = model.U_line_I(:, 1);
    G = model.G(:, i);
    t_Iapp_met = model.T_Iapp_met(i);
    array_Iapp = model.Iapp;
    Isyn_EE = model.Isyn_line_EE;
    Isyn_InE = model.Isyn_line_InE;
    Isyn_ExI = model.Isyn_line_ExI;
    Pre_EE = model.Pre_EE;
    Post_EE = model.Post_EE;
    Pre_InE = model.Pre_InE;
    Post_InE = model.Post_InE;
    Pre_ExI = model.Pre_ExI;
    Post_ExI = model.Post_ExI;
    
    W_EE = params.W_EE_0;
    W_InE = params.W_InE_0;

    numofPost_EE = model.numofPost_EE;
    I_poisson_noise = model.I_poisson_noise(:, i);
    ST_N_EE = model.ST_N_EE;
    ST_S_EE = model.ST_S_EE;
    ST_S_InE = model.ST_S_InE;
    ST_S_ExI = model.ST_S_ExI;
    Gli_global_expand = model.Gli_global_expand;
    if params.stimulation_mode % automatic
        I_stim_E = model.auto_electrode_neuro_E_stim_timed(:, i);
        I_stim_I = model.auto_electrode_neuro_I_stim_timed(:, i);
    else % manual
        I_stim_E = model.manual_electrode_neuro_E_stim_timed(:, i);
        I_stim_I = model.manual_electrode_neuro_I_stim_timed(:, i);
    end
    % T_Iapp = model.T_Iapp;
    
    %% -------------------------------------------------------------------------------------------
    
    %% initialization
    I_poisson_noise = double(I_poisson_noise);
    % I_poisson_noise = 0;
    
    % Input image as rectangle function of applied current to neuronal layer
    if t_Iapp_met == 0
        Iapp = zeros(params.mneuro_E, params.nneuro_E, 'uint8');
    else
        Iapp = array_Iapp(:, :, t_Iapp_met); % for the timeline of applied input
    end
    Iapp_line   = Iapp(:);
    Iapp_line   = double(Iapp_line);
    Gli_global_expand_line  = Gli_global_expand(:);
    
    %% Izhikevich neuron model
    % Excitatory neurons
    fired = find(V_E >= params.neuron_fired_thr); % This limit is needed to keep the differential equations functioning.
    V_E(fired) = params.c; % After-spike reset value of the membrane potential
    U_E(fired) = U_E(fired) + params.d; % After-spike reset value of the recovery variable
    if params.stimulation_neuro_E == 1
        I_sum_E  = Iapp_line + Isyn_EE + Isyn_InE + ST_N_EE' .* I_poisson_noise + I_stim_E; 
    else
        I_sum_E  = Iapp_line + Isyn_EE + Isyn_InE + ST_N_EE' .* I_poisson_noise; 
    end
    I_sum_E = min(I_sum_E, params.I_input_thr);    % ISUM = I_sum_E; % used for plot
    V_E    = V_E + params.step .* 1000 .* (0.04 .* (V_E .^ 2) + (5 .* V_E) + 140 + I_sum_E - U_E);
    U_E    = U_E + params.step .* 1000 .* params.aa .* (params.b .* V_E - U_E);
    V_E    = min(V_E, params.neuron_fired_thr);
    
    % Inhibitory neurons
    fired = find(V_I >= params.neuron_fired_thr);
    V_I(fired) = params.c;
    U_I(fired) = U_I(fired) + params.d;
    if params.stimulation_neuro_I == 1
        I_sum_I = Isyn_ExI + I_stim_I;
    else
        I_sum_I = Isyn_ExI;
    end
    V_I = V_I + params.step .* 1000 .* (0.04 .* (V_I .^ 2) + (5 .* V_I) + 140 + I_sum_I - U_I); 
    U_I = U_I + params.step .* 1000 .* params.aa .* (params.b .* V_I - U_I);
    V_I = min(V_I, params.neuron_fired_thr);

    %% Neuron synaptic currents
    V_line_E_mask = zeros(params.quantity_neurons_E, 1, 'logical');
    V_line_E_mask(V_E > params.neuron_fired_thr - 1) = 1;
    
    Isyn_EE = zeros(params.quantity_neurons_E, 1, 'double');
    Isyn_ExI = zeros(params.quantity_neurons_I, 1, 'double'); 
    Isyn_InE = zeros(params.quantity_neurons_E, 1, 'double'); 
    S_E = 1 ./ (1 + exp(( - V_E ./ params.Ssyn)));
    S_I = 1 ./ (1 + exp(( - V_I ./ params.Ssyn)));

    %% Isyn_EE
    numofIsyn_EE = zeros(params.quantity_neurons_E, 1, 'int8');
    W_EE_1 = ST_S_EE' .* W_EE .* (1 + Gli_global_expand_line(Post_EE) .* params.Eta_WM);
    Isync_EE = W_EE_1 .* S_E(Pre_EE) .* (0 - V_E(Post_EE));
    for q = 1 : params.quantity_connections_EE
        Isyn_EE(Post_EE(q)) = Isyn_EE(Post_EE(q)) + Isync_EE(q);
        numofIsyn_EE(Post_EE(q)) = numofIsyn_EE(Post_EE(q)) + int8((Isync_EE(q) > 0.001));
    end
    Isyn_EE(Isyn_EE < 0) = 0;
    
    %% Isyn_InE
    Isync_InE = ST_S_InE' .* W_InE .* S_I(Pre_InE) .* (params.Esyn_ - V_E(Post_InE));
    for q = 1 : params.quantity_connections_InE
        Isyn_InE(Post_InE(q)) = Isyn_InE(Post_InE(q)) + Isync_InE(q);
    end
    
    %% Isyn_ExI
    Isync_ExI = ST_S_ExI' .* params.W_ExI .* S_E(Pre_ExI) .* (0 - V_I(Post_ExI));
    for q = 1 : params.quantity_connections_ExI
        Isyn_ExI(Post_ExI(q)) = Isyn_ExI(Post_ExI(q)) + Isync_ExI(q);
    end    

    %% Glutamate (neurotransmitter model)
    SoI = double(numofIsyn_EE) ./ double(numofPost_EE); 
    G = G + params.step .* (params.r_Glu .* SoI - G ./ params.tau_Glu) .* (ST_N_EE)';
    
    %% -------------------------------------------------------------------------------------------
    
    %% outputs
    model.V_line_E(:, i + 1) = V_E;
    model.U_line_E(:, 1) = U_E;
    model.V_line_I(:, i + 1) = V_I;
    model.U_line_I(:, 1) = U_I;
    model.G(:, i + 1) = G;
    
    model.Isyn_line_EE = Isyn_EE;
    model.Isyn_line_InE = Isyn_InE;
    model.Isyn_line_ExI = Isyn_ExI;

    model.V_line_E_mask = V_line_E_mask;
    model.Iapp_v_full(:, :, i) = Iapp; % for video
end
