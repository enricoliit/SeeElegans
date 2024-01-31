% loading stack
% load('wholestack.mat');

% loading tracking results
% load('neurons_reconstructed_max_2.mat');

% specifying output folder
% outputpath='.\output\';

% specifying voxel size
voxel_size=[0.267 0.267 2];

[folders, files] = getFileStructureMat('E:\PanNeuronale\Acquisizioni\20230912');

%%

for i = 1 : 1%numel(files)
    % running manual correction GUI
    [folders{i} '\' files{i}{1}]
    load([folders{i} '\' files{i}{1}]);
    outputpath=[folders{i} '\output\'];
    parameter_selection_and_tracking_batch(wholestack,outputpath,voxel_size);
    drawnow
    uiwait(gcf);
end

%%

for i = 2 : numel(files)
    % running manual correction GUI
    [folders{i} '\' files{i}{1}]
    load([folders{i} '\' files{i}{1}]);
    outputpath=[folders{i} '\output\'];
    tracking_step(wholestack,outputpath);
    drawnow
    close all
end