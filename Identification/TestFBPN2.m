warning('off', 'all');

m = [];

%% Dataset 4

load('.\Dataset4.mat');
fname = 'D4'; fixed_time = 132; voxel_size = [0.267 0.267 2]; current_neuron_arangement = pca_reconstruction(neurons_cleaned,20); % prima sembrava giusto [29.7563 76.6338]; [319.0453  115.6386]; 
chosen_point = [40 70]; 
[correspondences,directions,guess,ruler]=automatic_identity(current_neuron_arangement,chosen_point,voxel_size,fixed_time,2400);

BW1 = [44 51 10 33];
BW2 = [66 20 9 14];

checkId(fname, guess(2,:), BW1, BW2);
disp(['CG1: ' num2str(guess(1,:)) ';   R: ' num2str(BW1)]);
disp(['CG2: ' num2str(guess(2,:)) ';   R: ' num2str(BW2)]);
disp(['r = ' num2str(ruler)])

title(fname)

[tp, ~] = match_neurons_with_correspondences(neurons_identified, correspondences);
m = [m tp];
%% Dataset 5

load('.\Dataset5.mat');
fname = 'D5'; fixed_time = 132; voxel_size = [0.267 0.267 2]; current_neuron_arangement = pca_reconstruction(neurons_cleaned,20); 
chosen_point = [30 70];
[correspondences,directions,guess,ruler]=automatic_identity(current_neuron_arangement,chosen_point,voxel_size,fixed_time,2400); % t = 132

BW1 = [41 65 61 64];
BW2 = [55 57 70];

checkId(fname, guess(2,:), BW1, BW2);
disp(['CG1: ' num2str(guess(1,:)) ';   R: ' num2str(BW1)]); 
disp(['CG2: ' num2str(guess(2,:)) ';   R: ' num2str(BW2)]); 
disp(['r = ' num2str(ruler)])

title(fname)

[tp, ~] = match_neurons_with_correspondences(neurons_identified, correspondences);
m = [m tp];

%% Dataset 6

load('.\Dataset6.mat');
fname = 'D6'; fixed_time = 132; voxel_size = [0.267 0.267 2]; current_neuron_arangement = pca_reconstruction(neurons_cleaned,20); 
chosen_point = [29 76];
[correspondences,directions,guess,ruler]=automatic_identity(current_neuron_arangement,chosen_point,voxel_size,fixed_time,2400); % t = 132

BW1 = [25 17 74 63];
BW2 = [9 106 58 48];

checkId(fname, guess(2,:), BW1, BW2);
disp(['CG1: ' num2str(guess(1,:)) ';   R: ' num2str(BW1)]); disp(['CG2: ' num2str(guess(2,:)) ';   R: ' num2str(BW2)]); 
disp(['r = ' num2str(ruler)])

title(fname)

[tp, ~] = match_neurons_with_correspondences(neurons_identified, correspondences);
m = [m tp];

%% Dataset 2

load('.\Dataset2.mat');

fname = 'D2'; fixed_time = 132; voxel_size = [0.356 0.356 2]; current_neuron_arangement = pca_reconstruction(neurons_cleaned,20); 
chosen_point = [380 80];
[correspondences,directions,guess,ruler]=automatic_identity(current_neuron_arangement,chosen_point,voxel_size,fixed_time,1100); % t = 132

BW1 = [77 68 76 71];
BW2 = [40 89 90 111];

checkId(fname, guess(2,:), BW1, BW2);
disp(['CG1: ' num2str(guess(1,:)) ';   R: ' num2str(BW1)]); disp(['CG2: ' num2str(guess(2,:)) ';   R: ' num2str(BW2)]); 
disp(['r = ' num2str(ruler)])

title('D2');

[tp, ~] = match_neurons_with_correspondences(neurons_identified, correspondences);
m = [m tp];

%% Dataset 1

load('.\Dataset1.mat');

fname = 'D1'; fixed_time = 132; voxel_size = [0.267 0.267 2]; current_neuron_arangement = pca_reconstruction(neurons_cleaned,20); 
chosen_point = [30 70];
[correspondences,directions,guess,ruler]=automatic_identity(current_neuron_arangement,chosen_point,voxel_size,fixed_time,2400); % t = 132

BW1 = [52 94 126 82];
BW2 = [75 128 54 40]; 

checkId(fname, guess(2,:), BW1, BW2);
disp(['CG1: ' num2str(guess(1,:)) ';   R: ' num2str(BW1)]); disp(['CG2: ' num2str(guess(2,:)) ';   R: ' num2str(BW2)]); 
disp(['r = ' num2str(ruler)])


[tp, ~] = match_neurons_with_correspondences(neurons_identified, correspondences);
m = [m tp];

disp(["TPr = " num2str(mean(m))]);

%%
function [check1, check2]=checkId(fname, guess, BW1, BW2)
okmessagge = ['<strong>' fname ' OK!</strong> \n'];
errormessage = ['<strong>' fname ' Errore</strong> \n'];
check = 0;
for j = 1:size(guess,2)
    if ismember(guess(1,j),BW1)
        check = check + 1;
    end
end
if check == 3
    fprintf(okmessagge);
else
    check1 = 0;
    for j = 1:size(guess,2)
        if ismember(guess(1,j),BW2)
            check1 = check1 + 1;
        end
    end
    if check1 == 3
        fprintf(okmessagge);
    else
        fprintf(errormessage);
    end
end
end