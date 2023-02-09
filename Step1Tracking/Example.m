% loading stack
load('.\..\wholestack.mat');

% specifying output folder
outputpath='.\..\output\';

% running parameter selection step
SE_parameter_selection(wholestack, outputpath);

% running tracking step
SE_tracking_step(wholestack, outputpath);