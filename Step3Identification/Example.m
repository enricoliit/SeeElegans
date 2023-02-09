% loading stack
load('.\..\wholestack.mat');

% loading cleaned tracks
load('.\neurons_cleaned.mat');

% specifying output folder
outputpath='.\..\output\';

% running the identification step
identification_step(wholestack, neurons_cleaned, outputpath,[0.267 0.267 2]);
