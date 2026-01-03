function [Ca_expanded, Gli_global_expand] = ...
    expand_astrocytes(Ca, Gli_global)
    % Expands astrocyte calcium concentration and electrical currents to connected neurons
    params = model_parameters();
    Ca_expanded = zeros(params.mneuro_E, params.nneuro_E);
    Gli_global_expand = zeros(params.mneuro_E, params.nneuro_E);

    for J = 1:params.mastro
        for K = 1:params.nastro
            % Extract the corresponding value for a 2x2 block
            Ca_block = Ca(J, K);
            Gli_global_block = Gli_global(J, K);
    
            % Apply the max operation to the corresponding 2x2 region 
            Ca_expanded(J:J+1, K:K+1) = max(Ca_expanded(J:J+1, K:K+1), Ca_block);
            Gli_global_expand(J:J+1, K:K+1) = max(Gli_global_expand(J:J+1, K:K+1), Gli_global_block);
        end
    end
end
