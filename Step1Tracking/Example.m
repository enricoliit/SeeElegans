% loading stack
load('.\..\wholestack.mat');

% specifying output folder
outputpath='.\..\output\';

% running parameter selection step
SE_parameter_selection(wholestack, outputpath, [0.267 0.267 2]);
