% loading stack
% load('.\..\wholestack.mat');

% loading cleaned tracks
% load('.\neurons_cleaned.mat');

% specifying output folder
outputpath='E:\GNAO experiments\E246K\v1\';
% outputpath='.\';

% running the identification step
identification_step(wholestack, neurons_cleaned, outputpath, [0.267 0.267 2]);

% MA IL V4 E IL V2 SONO IDENTICI? E246K


