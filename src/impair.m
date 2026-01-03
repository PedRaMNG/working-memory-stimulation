function [ST_S, ST_N] = impair(n, params, Damage_pattern, numofPost, Post)
    sel = params.impairmode;
    switch sel
        case 0 % Healthy
            ST_S = ones(1, n);
			ST_N = ones(1, params.quantity_neurons_E);

        case 1 % Random damage to the synapses with fixed amplitude

            numofper = floor(params.scenario1a2.probability * n);
            h1 = [params.scenario1a2.amplitude*ones(1, numofper), ones(1, n - numofper)];
            ST_S = h1(randperm(n));
            
            ST_N = ones(1, params.quantity_neurons_E);
            count = 0;
            for i = 1:numel(numofPost)
                ST_N(i) = sum(ST_S(count + (1:numofPost(i))))/numofPost(i);
                count = count + numofPost(i);
            end
			% numofper = floor(params.scenario1a2.probability * params.quantity_neurons_E);
            % h1 = [params.scenario1a2.amplitude*ones(1, numofper), ones(1, params.quantity_neurons_E - numofper)];
            % ST_N =  h1(randperm(params.quantity_neurons_E));

        case 2 % Random damage to the synapses with random amplitude

            numofper = floor(params.scenario1a2.probability * n);
            h1 = [(1-params.scenario1a2.amplitude)*rand(1, numofper) + params.scenario1a2.amplitude, ...
                         ones(1, n - numofper)];
            ST_S = h1(randperm(n));
            
			numofper = floor(params.scenario1a2.probability * params.quantity_neurons_E);
            h1 = [(1-params.scenario1a2.amplitude)*rand(1, numofper) + params.scenario1a2.amplitude, ...
                         ones(1, params.quantity_neurons_E - numofper)];
            ST_N = h1(randperm(params.quantity_neurons_E));

        case 3 % Damage with regard to the damaging pattern

            % Damage_pattern = double(Damage_pattern)./255;
            % ST_S = ones(1,n);
			% cumul= 0;
			% for nn = 1 : params.quantity_neurons_E
            %     ST_S(cumul + (1:numofPost(nn))) = Damage_pattern(nn)*ones(1, floor(numofPost(nn)));
			% 	cumul = cumul + numofPost(nn);
			% end
			% [~, indi] = sort(Pre2);
			% ST_S = ST_S(indi);
            ST_S = double(Post);
            tumVec = double(Damage_pattern(:))/255;
            for i = 1:numel(tumVec)
                ST_S(ST_S==i) = tumVec(i);
            end
            ST_N = double(Damage_pattern(:))'/255;
		case 31
			ST_S = double(Damage_pattern(:))'/255;
            ST_N = ST_S;
    end
end
