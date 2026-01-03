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