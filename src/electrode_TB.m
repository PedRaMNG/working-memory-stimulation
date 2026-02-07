% clear;
% close all;
% clc;

% mneuro_E = 26;
% nneuro_E = 26;
% electrode_lambda = 0.4;
% electrode_grid = [16,16];


%---------------------------------
% mneuro_E = 43;
% nneuro_E = 43;
% mastro = 14;
% nastro = 14;
% az = 3;

% electrode_lambda = 3; %3.5 without offset %1.4 %7
% electrode_grid = [4,4];

% electrode_lambda = 1.45; %1.6 without offset
% electrode_grid = [8,8];
%---------------------------------
colrs = [255, 255, 255; 29, 67, 80; 237, 31, 30] / 255;
locs = [0, 0.5, 1];
pro.cmap = makeGradient(256, colrs, locs, 0);
%---------------------------------
mneuro_E = 71; %71 35
nneuro_E = 71; %71 35
mneuro_I = 35; 
nneuro_I = 35; 
mastro = 35;
nastro = 35;
az = 2;

scale_x_E = max(mneuro_E, mneuro_I) / mneuro_E;
scale_y_E = max(nneuro_E, nneuro_I) / nneuro_E;
scale_x_I = max(mneuro_E, mneuro_I) / mneuro_I;
scale_y_I = max(nneuro_E, nneuro_I) / nneuro_I;

electrode_lambda = 2.5; %1.4 %7
electrode_grid = [8,8];

electrode_lambda = 5.5; %1.4 %7
electrode_grid = [4,4];
%---------------------------------

% Output: [a matrix in n*2 which n is the number of electrodes]
% electrode_grid = [2,2] OR [6,3] OR ...


electrode_positions = zeros(prod(electrode_grid), 2);

% Define the offset values (you can adjust these as needed)
x_offset = 2.5;  % example offset in x direction
y_offset = 2.5;  % example offset in y direction

% Calculate the spacing between electrodes
x_sub_net = (mneuro_E - 2*x_offset) / electrode_grid(1);
y_sub_net = (nneuro_E - 2*y_offset) / electrode_grid(2);

% Compute the electrode positions with the specified offsets
for i = 1:electrode_grid(1) % x
    for j = 1:electrode_grid(2) % y
        electrode_positions((i-1)*electrode_grid(2) + j, 1) = x_offset + (i - 0.5) * x_sub_net;
        electrode_positions((i-1)*electrode_grid(2) + j, 2) = y_offset + (j - 0.5) * y_sub_net;
    end
end
    


electrode_sensitivity_neuron_E = zeros(mneuro_E, nneuro_E, size(electrode_positions, 1));
electrode_sensitivity_neuron_I = zeros(mneuro_I, nneuro_I, size(electrode_positions, 1));
for electrode = 1:size(electrode_positions, 1)
    for i = 1:mneuro_E
        for j = 1:nneuro_E
        d = sqrt((i - electrode_positions(electrode, 1))^2 + (j - electrode_positions(electrode, 2))^2);
        % Let electrode_lambda = 0.05 for 43*43 network
        % electrode_sensitivity_neuron_E(i, j, electrode) = exp(-0.5*d/electrode_lambda);
        electrode_sensitivity_neuron_E(i, j, electrode) = exp(-0.5*(d/electrode_lambda)^2);
        end
    end
end
% Calculate the sensitivity for the inhibitory network
for electrode = 1:size(electrode_positions, 1)
    for i = 1:mneuro_I
        for j = 1:nneuro_I
            d = sqrt((i * scale_x_I - electrode_positions(electrode, 1))^2 + (j * scale_y_I - electrode_positions(electrode, 2))^2);
            electrode_sensitivity_neuron_I(i, j, electrode) = exp(-0.5 * (d / electrode_lambda)^2);
        end
    end
end

hFigSIngle = figure();
% surf(sum(electrode_sensitivity_neuron_E(:,:,[1,electrode_grid(1),electrode_grid(1)*electrode_grid(2)-electrode_grid(1)+1,electrode_grid(1)*electrode_grid(2)]),3));
surf(sum(electrode_sensitivity_neuron_E(:,:,[5,6,7, 9,10,11]),3));
colormap(pro.cmap);
newPosition = [750, 300, 600, 400];
set(hFigSIngle, 'Position', newPosition);
% view(90,90);

%%
I_app_astro = zeros(params.mastro, params.nastro);
I_stim = sum(electrode_sensitivity_neuron_E(:,:,[5,6,7, 9,10,11]),3);


for j = 1:params.mastro
    for k = 1:params.nastro
        
        I_app_astro(j,k) = sum((I_stim(j:j+1, k:k+1) > 0.2), 'all')>=2;

    end
end
figure();
imagesc (I_app_astro);
%%

hFig1 = figure();
variablesToSum = [1,5];
% variablesToSum = [5:8];
% surf(sum(electrode_sensitivity_neuron_E(:,:,variablesToSum), 3)); % selective
% surf(sum(electrode_sensitivity_neuron_E, 3)); % all
%     axis off; % Remove axis and labels
%     view(3); % Set 3D view
%     grid off; % Turn off grid
%     shading interp;
%     colormap(pro.cmap);

imagesc(sum(electrode_sensitivity_neuron_E, 3)); % Plot the data as an image
    set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
    colormap(pro.cmap); % Apply the colormap
    axis equal tight;

% view(90,90);
colorbar;
newPosition = [150, 300, 600, 400]; % [left, bottom, width, height]
set(hFig1, 'Position', newPosition);

hFig2 = figure();
surf(sum(electrode_sensitivity_neuron_E, 3)); % all
    axis off; % Remove axis and labels
    view(3); % Set 3D view
    grid off; % Turn off grid
    shading interp;
    colormap(pro.cmap);
% imagesc(sum(electrode_sensitivity_neuron_I, 3)); % Plot the data as an image
%     set(gca, 'XTick', [], 'YTick', []); % Remove tick marks
%     colormap(pro.cmap); % Apply the colormap
%     axis equal tight;
newPosition = [150, 300, 600, 400]; % [left, bottom, width, height]
set(hFig2, 'Position', newPosition);


% imwrite(sum(electrode_sensitivity_neuron_E, 3), 'SavedData/electrode_record_71x71_offset.jpg');

%% astrocyte mode 1
sum_electrode_sensitivities = sum(electrode_sensitivity_neuron_E, 3);
sum_electrode_sensitivities_astro = zeros(mastro, nastro);
sj = 0;

for j = 1 : az : (mneuro_E - az)
   sk = 0;
   for k = 1 : az : (nneuro_E - az)
       % Number of neurons that passed the glu threshold
       sum_electrode_sensitivities_astro(j - sj, k - sk) = ... 
           sum(sum_electrode_sensitivities(j : j + az, k : k + az), 'all');
        if az == 2
            sk = sk + 1;
        else
            sk = sk + 2;
        end
    end
    if az == 2
        sj = sj + 1;
    else
        sj = sj + 2;
    end
end


% hFig1 = figure();
% surf(sum_electrode_sensitivities_astro);
% colormap(hsv);
% % view(90,90);
% colorbar;
% newPosition = [150, 300, 600, 400]; % [left, bottom, width, height]
% set(hFig1, 'Position', newPosition)
% newPosition = [750, 300, 600, 400]; % [left, bottom, width, height]
% set(hFig1, 'Position', newPosition)

%% astrocyte mode 2 - main one
electrode_sensitivity_astro = zeros(mastro, nastro, size(electrode_positions, 1));

for electrode = 1:size(electrode_positions, 1)
    sj = 0;
    for j = 1 : az : (mneuro_E - az)
       sk = 0;
       for k = 1 : az : (nneuro_E - az)
            electrode_sensitivity_astro(j - sj, k - sk, electrode) = ... 
                sum(electrode_sensitivity_neuron_E(j : j + az, k : k + az, electrode), 'all');
            if az == 2
                sk = sk + 1;
            else
                sk = sk + 2;
            end
        end
        if az == 2
            sj = sj + 1;
        else
            sj = sj + 2;
        end
    end
end

% hFig2 = figure();
% surf(sum(electrode_sensitivity_astro, 3));
% colormap(hsv);
% view(90,90);
% colorbar;
% newPosition = [750, 300, 600, 400]; % [left, bottom, width, height]
% set(hFig2, 'Position', newPosition);


