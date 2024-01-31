% loading stack
load('wholestack.mat');

% loading tracking results
load('neurons_reconstructed_max_2.mat');

% specifying output folder
outputpath='.\output\';

% specifying voxel size
voxel_size=[0.267 0.267 2];

% running manual correction GUI
parameter_selection_and_tracking(wholestack,outputpath,voxel_size);