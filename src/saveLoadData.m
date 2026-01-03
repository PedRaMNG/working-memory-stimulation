%% save data
% Generate a folder name with today's date and time
folderName = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
folderPath = fullfile('SavedData', char(folderName));  % Convert to char for folder path
mkdir(folderPath);

% Get the names of all variables in the 'model' struct
varss = fieldnames(model);

% Loop over each variable and save it
for i = 1:numel(varss)
    eval([varss{i}, ' = model.', varss{i}, ';']);
    save(fullfile(folderPath, [varss{i}, '.mat']), varss{i}, '-v7.3');
    disp(varss{i});
end

% Save the 'params' variable
save(fullfile(folderPath, 'params.mat'), 'params', '-v7.3');
disp('params');
save(fullfile(folderPath, 'memory_performance.mat'), 'memory_performance', '-v7.3');
disp('memory_performance');

disp('------ Finished Saving ------');

%% Import data
% cd(fileparts(mfilename('fullpath')))
clearvars;
folderPath = 'SavedData/2026-01-01_20-35-23/'; % put your saved data here and load them from this folder
fileList = dir(fullfile(folderPath, '*.mat'));
numFiles = numel(fileList);
fprintf('Locating files...\n');
for i = 1:numFiles
    filename = fullfile(folderPath, fileList(i).name);
    load(filename);  
    clc; fprintf('Loading %d / %d\r', i, numFiles);
end
fprintf('\nDone loading.\n');
clear filename folderPath
% load data
workspaceVariables = who; % Get a list of variable names in the workspace
workspaceVariables = workspaceVariables(~strcmp(workspaceVariables, 'params'));
workspaceVariables = workspaceVariables(~strcmp(workspaceVariables, 'memory_performance'));
Mymodel = struct(); % Create an empty structure named 'model'
for i = 1:numel(workspaceVariables)
    variableName = workspaceVariables{i};
    Mymodel.(variableName) = eval(workspaceVariables{i});
end
clear(workspaceVariables{:});
eval('model = Mymodel;');
clear Mymodel i variableName workspaceVariables;
clc;










