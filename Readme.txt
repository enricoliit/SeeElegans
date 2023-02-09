Requirements:
Thi software run on MATLAB and requires Image Processing Toolbox
It has been tested on MATLAB 9.12 Image Processing Toolbox 11.5

How to use:
1) run SE_parameter_selection(input_stack, output_folder) in folder Step1Tracking
   Give as input the image stack (input_stack) and the path to the output folder (output_folder). 
   The stack requires the shape (x,y,z,t). 

   The function will open a window to inspect the stack and to set the parameters. 
   After running, the output folder will contain 4 files named: parameters.mat, parameters2.mat, parameters3.mat, 
   and parameters4.mat. 
   The set of parameters includes:
	a) the size and sigma of the log filter and the intensity threshold of the spots found through log filtering
	b) the selection of a video crop to test the tracking parameters
	c) the maximum time gap (in frames) and the maximum distance (in pixels) to link spots, and the minimum length to 
	   retain tracks
	d) the distance maximum gap (in pixels) to link spots based on the preservation of the distances among spots

2) run SE_tracking_step(input_stack, output_folder) in folder Step1Tracking
   Give as input the image stack (input_stack) and the path to the output folder (output_folder). 
   The stack requires the shape (x,y,z,t). The output folder must contain the files parameter.mat, parameter2.mat, 
   parameter3.mat, and parameter4.mat. 
   
   The function will perform the tracking step and output the variable 
   neurons_reconstructe_max_2.mat. This will contain the final sets of tracked spots. The function also opens a figure to show
   the resulting tracks and the unfiltered set of signals.

3) run remove_and_add_neurons(input_stack, neurons_reconstructed_max_2.mat, outputfolder, inversion)  in folder Step2Correction
   Give as input the image stack (input_stack), the final set of tracked spots (neurons_reconstructed_max2.mat) 
   and the path to the output folder (outputfolder). You may also set the inversion parameter that inverts the x and y 
   coordinates, useful in case the data come from a different tracking procedure.
   
   The function will open a GUI to inspect the tracking results. A mouse scroll will change the visualized time, 
   while zooming or ctrl+mouse scroll will change the visualized plane. By clicking on a spot, it is possible to visualize the 
   corresponding signal and to include it in the neurons to be removed by clicking on "Delete". Neurons will be removed only when
   clicking on the "Save" button will delete the neuron from the set of tracked spots. By clicking on "Add" button, it is
   possible to click on the image and by clicking on "Track" it is possible to add a spot by reconstructing its track through its
   neighbouring ones. By clicking Save the new spot will be displayed along with the other.
   Every time the button "Save" is clicked, the program checks for the existence of files called 
   "neurons_reconstructed_max_*.mat" in the output folder and will add a new file called "neurons_reconstructed_max_[N+1].mat",
   where N is the number of files with that name. So, the file with the highest number is the most recent one. The newly saved
   files will contain a variable called "neurons_verified", indicating that the original set of track has been modified. 
 
4) run identification_step(input_stack, neurons_cleaned, outputpath, voxel_size) in folder Step3Identification
   Give as input the image stack (input_stack), the cleaned set of tracked spots (neurons_cleaned.mat), the path to the output 
   folder (outputfolder) and the voxel size in um as [x, y, z]; if the voxel size is not specified the default value is 
   [0.267 0.267 2].
   
   The function will open a GUI to assign identities. A mouse scroll will change the visualized time, while zooming or ctrl+mouse
   scroll will change the visualized plane. The lower part will show a slice, the upper part is a 3D plot that can be moved
   By clicking on a spot, it is possible to visualize the corresponding signal. You can specify the direction of the head by clicking
   on the button "head point" and then clicking on the slice on a pixel close to the tip of the animal. This will identify the
   direction of the head of the animal. By clicking on auto-assign an identity hypothesis will be generated for a subset of neurons
   according to the signal coherence between neurons and their arrangement. By clicking on "Assign" you can choose a spot and assign
   an identity with the popup menu on the right part of the GUI. If the spot already has an identity it can be deleted by clicking on
   the "Delete" button. By clicking on "Store" the plots will be updated. By clicking on "Save" the current identification hypothesis
   will be saved. Every time the button "Save" is clicked, the program checks for the existence of files called 
   "neurons_reconstructed_max_id_*.mat" in the output folder and will add a new file called "neurons_reconstructed_max_id_[N+1].mat",
   where N is the number of files with that name. So, the file with the highest number is the most recent one. The newly saved
   files will contain a variable called "neurons_identified", reporting the current identification hypothesis in a cell.The cell will
   contain ordered ids in the "id" variable and names in the "name" variable ordered in the same way.
   



Example of use: the "wholestack.mat" file available at https://figshare.com/articles/media/wholestack_mat/22063073 contains a variable
called wholestack.mat that represents a calcium imaging 4D recording of a C. elegans with pan-neuronal GCaMP expression. If it is 
moved in the parent folder of the project the following lines may be run.

Step 1) commands to run the identification step when in ./Step1Tracking folder

	% loading stack
	load('.\..\wholestack.mat');

	% specifying output folder
	outputpath='.\..\output\';

	% running parameter selection step
	SE_parameter_selection(wholestack, outputpath);

	% running tracking step
	SE_tracking_step(wholestack, outputpath);


Step 2) commands to run the validation step when in ./Step2Correction folder

	% loading stack
	load('.\..\wholestack.mat');

	% loading tracking results
	load('.\..\output\neurons_reconstructed_max_2.mat');

	% specifying output folder
	outputpath='.\..\output\';

	% running manual correction GUI
	remove_and_add_neurons(wholestack,neurons_reconstructed_max_2,outputpath);



Step 3) commands to run the identification step when in ./Step3Identification folder

	% loading stack
	load('.\..\wholestack.mat');

	% loading cleaned tracks
	load('.\neurons_cleaned.mat');

	% specifying output folder
	outputpath='.\..\output\';

	% running the identification step
	identification_step(wholestack, neurons_cleaned, outputpath,[0.267 0.267 2]);



