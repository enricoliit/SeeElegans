% loading stack
load('.\..\wholestack.mat');

% loading tracking results
load('.\..\output\neurons_reconstructed_max_2.mat');

% specifying output folder
outputpath='.\..\output\';

% running manual correction GUI
remove_and_add_neurons(wholestack,neurons_reconstructed_max_2,outputpath);