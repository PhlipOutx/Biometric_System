Biometric System
v1.4
7/9/12

Change Log
----------
@v1.4
- Integrated libsvm
@v1.3
- Fixed bugs associated with SIFT and SURF when 0 keypoints were found
- Added v1.2 functionality for LDA
- Added <select_type> and <ellipse_type> subfields
- Added ICA as a feature type. 
@v1.2
- Added functionality to handle class ids by file so that the system will
  work with the Pinellas data
@v1.1
- Fixed situational errors with displaying results
- Improved output of variable_selection.m
- Allowed SIFT and SURF to store distance matrix incrementally so that it
  works with the Palmetto Cluster
- Store .mat files in -v7.3 format
- Save .mat file of results structure



Table of Contents:
Section 1.0		Available Functions
Section 2.0		Experiment Files
Section 3.0		Outputs
Section 4.0		Usage
Section 5.0		Features To Add



Section 1.0   -   Available Functions
-------------------------------------
%%%%
function install_bprl()
Adds all subdirectories to the path. This function must be called with each
new instance of the Matlab environment.

%%%%
function results = biometric_system(input)
% input
input(struct):		Matlab structure containing fields that hold variables
                    used in a biometric experiment. See Section 2.0 for
                    list of required fields.
                    
input(string):		Filepath to xml file that will be parsed into a struct.

%ouput
results:			See Section 3.0 for list of outputs.

%%%%
function estimate_memory(input)
%input
input(struct):		Matlab structure containing fields that hold variables
                    used in a biometric experiment. See Section 2.0 for
                    list of required fields.
                    
input(string):		Filepath to xml file that will be parsed into a struct.

%output
Prints estimated memory to be used during each part of the experiments.
Ensure that there is enough memory on the machine running the experiment to
accomidate the experiment. 

%%%%
function results = twoset_fusion(input)
% input
input(struct):		Matlab structure containing fields that hold variables
                    used in a fusion of two experiments. See Section 2.0
                    for list of required fields.
                    
input(string):		Filepath to xml file that will be parsed into a struct.

%ouput
results:			See Section 3.0 for list of outputs.

%%%%
function variable_selection(input)
% input
input(struct):		Matlab structure containing fields that hold variables
                    used in a biometric experiment. See Section 2.0 for
                    list of required fields.
                    
input(string):		Filepath to xml file that will be parsed into a struct.

%output
Prints out two 9x5 matrices. One contains Rank-1 recognition rates and the
other contains equal error rates. This function tests an experiment with
variable image and patch size. The available patch sizes (columns) are 
[15, 20, 25, 30, 35] and the available image sizes (rows) are
[4;5;6;7;8;9;10;11;12]*patch_size



Section 2.0   -   Experiment Files
----------------------------------
Each of the functions in this system use XML files as input. These XML
files are found in the Experiments/ folder. Variables for an experiments
are placed inside an <experiment> field. This field has 4 required
subfields:
    <id>
    <variables>
    <input>
    <output>
    
Section 2.1 - <id> field
--------------------------
The <id> field must contain an alphanumeric identifier that is unique for
each experiment. I suggest that the identifier start with a letter because
the system will truncate any leading zeros. 

Section 2.2 - <variables> field
---------------------------------
The <variables> field will have different requirements depending on the
function called. biometric_system(input), estimate_memory(input), and
variable_selection(input) have the following required subfields:
    <nameprefix> % A number representing how many digits are used in the
                   subject ID portion of an image filename.
    <distance> % A string representing the desired distance function to
                 use. Options are eucdist, sqdist, dotprod, nrmcorr,
                 corrdist, angle, cityblk, maxdiff, mindiff, intersect,
                 intersectdis, chisq, kldiv, jeffrey. See
                 Similarity/slmetric_pw.m for definitions of each.
    <printing> % Number of desired printing method of iterations.
                 0 - do not print iterations of loops.
                 1 - print in a manner handled by the matlab environment.
                 2 - print in a manner handled by the command line.
    <select_type> % Only required when running variable_selection
                 1 - Use Rank-1 to decide best variables.
                 2 - Use EER to decide best variables.
    <feature> % Contains variables for chosen feature extraction method.
                Has subfield <name> which will require different addition
                subfields of <feature> depending on the value of name.
        If <name>LBP</name> include:
            <radius>
            <samples>
            <type>
            % See Features/Local/lbp.m and Util/getmapping.m for
              description of fields.
        If <name>HOG</name> include:
            <bins> % number of histogram bins
        If <name>LPQ</name> include:
            <winSize>
            <decorr>
            <freqestim> 
            <mode>
            % See Features/Local/lpq.m for description of fields.
        If <name>WLD</name> include:
            <num_orientation>
            <num_excitation>
            % See Features/Local/wld.m for description of fields.
        If <name>PCA</name> include:
            <vectors> % number of eigen vectors to retain
        If <name>LDA</name> include:
            <vectors> % number of eigen vectors to retain
        If <name>ICA</name> include:
            <vectors> % number of eigen vectors to retain
        If <name>KPCA</name> include:
            <vectors> % number of eigen vectors to retain
        If <name>KFA</name> include:
            <vectors> % number of eigen vectors to retain
        If <name>SIFT</name> or <name>SURF</name> there are no
        required fields.
    <patches> % Has subfields for defining block/patch configuration when
                using certain features. 
        <type> % Only possible value is "nonover" at the moment.
        <size_x> % number of pixels in the x direction
        <size_y> % number of pixels in the y direction
    <data> % Has subfields for defining preprocessing steps on data.
        <channel> % color channel to use
        <max_x> % number of columns to resize all of the data to
        <max_y> % number of rows to resize all of the data to
        <histeq> % 0 or 1 to perform histogram equalization
        <ellipse_x> % x loction of gray color ellipse overlain on image
        <ellipse_y> % y loction of gray color ellipse overlain on image
        <ellipse_a> % major axis of ellipse
        <ellipse_b> % minor axis of ellipse
        <ellipse_t> % theta of ellipse
        <ellipse_type> % 1 for shading inside ellispe and 2 for outside
        
twoset_fusion(input) have the following required subfields:
    <nameprefix> % A number representing how many digits are used in the
                   subject ID portion of an image filename.
    <distance> % A string representing the desired distance function to
                 use. Options are eucdist, sqdist, dotprod, nrmcorr,
                 corrdist, angle, cityblk, maxdiff, mindiff, intersect,
                 intersectdis, chisq, kldiv, jeffrey. See
                 Similarity/slmetric_pw.m for definitions of each.
    <printing> % Number of desired printing method of iterations.
                 0 - do not print iterations of loops.
                 1 - print in a manner handled by the matlab environment.
                 2 - print in a manner handled by the command line.
    <id1> % ID of an experiment to fuse
    <id2> % ID of an experiment to fuse
    <fusion> % type of fusion to perform
               Possible values: "feature", "score"
    <feature> % Has only 1 subfield <name> which contains type of feature
                originally used in the experiments.

Section 2.3 - <input> field
---------------------------
The <input> field contains filepaths to imagelists and directories where
the images are located. The subfields are as follows:
    <training> % filepath to imagelist of training data
    <training_id> % filepath to list of class ids of training data
    <probe> % filepath to imagelist of probe data
    <probe_id> % filepath to list of class ids of probe data
    <gallery> % filepath to imagelist of gallery data
    <gallery_id> % filepath to list of class ids of gallery data
    <datadir> % directory containing images in the imagelists
                Must have trailing slash.

Section 2.4 - <output> field
----------------------------
The <output> field contains variables to direct the output of the system.
The subfields are as follows:
    <resultsdir> % directory to store output
                   Must have trailing slash.
    <training_mat> % filename of matrix to store training features in
    <probe_mat> % filename of matrix to store probe features in
    <gallery_mat> % filename of matrix to store gallery features in
    <distance_mat> % filename of matrix to store distances in
    <results_mat> % filename of matrix to store results in
    <roc_resolution> % number of point in the ROC and DET curves
    <make_files> % 0 or 1 to save points on curve to file
    <make_figures> % 0 or 1 to save figures
    <plot_color> % color of line in figure
    <plot_linewidth> % width of line in figure



Section 3.0   -   Outputs
-------------------------
biometric_system(input) and twoset_fusion(input) return a struct containing
the following fields:
    CMC_rec_rates - A vector containing the recognition rate at the Rank of
                    index.
    ROC_ver_rate - A vector of size 1 x (resolution+1) containing the 
                   verification rates of the ROC curve evaluated at 
                   the corresponding treshold values; the values in
                   the vector stand for 1-FRR; where FRR is the
                   false rejection rate
    ROC_miss_rate - A vector of size 1 x (resolution+1) containing the miss 
                    verification rates of the ROC curve evaluated at 
                    the corresponding treshold values; the values in
                    the vector stand for 1-FAR; where FAR is the
                    false acceptance rate
    ROC_rates_and_threshs - A struct with data computed at characteristic
                            operating points on the ROC curve; the struct
                            contains the following fields: 
                         
        .EER_er       - equal error rate - the value where FAR=FRR
        .EER_tr       - the threshold needed to obtain the EER 
        .EER_frr      - the value of the FRR at the EER 
        .EER_ver      - the value of the verification rate at the EER                      
        .EER_far      - the value of the FAR at the EER 

        .FRR_01FAR_er  - the value of the half total error rate at the 
                      ROC operating point that ensures that the false 
                      acceptance rate is 10 times higher than the 
                      false rejection rate; i.e., FAR=0.1FRR
        .FRR_01FAR_tr  - the threshold needed to obtain FAR=0.1FRR 
        .FRR_01FAR_frr - the value of the FRR at FAR=0.1FRR 
        .FRR_01FAR_ver - the value of the verification rate at FAR=0.1FRR                      
        .FRR_01FAR_far - the value of the FAR at FAR=0.1FRR 

        .FRR_10FAR_er  - the value of the half total error rate at the 
                      ROC operating point that ensures that the false 
                      rejection rate is 10 times higher than the 
                      false acceptance rate; i.e., FAR=10FRR
        .FRR_10FAR_tr  - the threshold needed to obtain FAR=10FRR 
        .FRR_10FAR_frr - the value of the FRR at FAR=10FRR 
        .FRR_10FAR_ver - the value of the verification rate at FAR=10FRR                      
        .FRR_10FAR_far - the value of the FAR at FAR=10FRR

        .VER_001FAR_er  - the value of the half total error rate at the 
                       ROC operating point where the FAR equals 0.01% 
        .VER_001FAR_tr  - the threshold needed to a FAR of 0.01%
        .VER_001FAR_frr - the value of the FRR a FAR of 0.01% 
        .VER_001FAR_ver - the value of the verification rate at a FAR of 
                       0.01%                      
        .VER_001FAR_far - the value of the FAR at a FAR of 0.01% (note 
                       that this is the actual value of the FAR, since 
                       there might not be anough data to obtain an FAR 
                       of exactly 0.01%)

        .VER_01FAR_er  - the value of the half total error rate at the 
                      ROC operating point where the FAR equals 0.1% 
        .VER_01FAR_tr  - the threshold needed to a FAR of 0.1%
        .VER_01FAR_frr - the value of the FRR a FAR of 0.1% 
        .VER_01FAR_ver - the value of the verification rate at a FAR of 
                      0.1%                      
        .VER_01FAR_far - the value of the FAR at a FAR of 0.1% (note 
                      that this is the actual value of the FAR, since 
                      there might not be anough data to obtain an FAR 
                      of exactly 0.1%)

        .VER_1FAR_er  - the value of the half total error rate at the 
                      ROC operating point where the FAR equals 1% 
        .VER_1FAR_tr  - the threshold needed to a FAR of 1%
        .VER_1FAR_frr - the value of the FRR a FAR of 1% 
        .VER_1FAR_ver - the value of the verification rate at a FAR of 
                      1%                      
        .VER_1FAR_far - the value of the FAR at a FAR of 1% (note 
                      that this is the actual value of the FAR, since 
                      there might not be anough data to obtain an FAR 
                      of exactly 1%)

        .dprime - The d' value of the distributions



Section 4.0   -   Usage
-----------------------
The source files are placed in subfolders to maintain an order to the
file structure of the system. This forces me to include the subfolders
in the Matlab path before using any of the high level functions.
install_bprl() must be called with each new instance of Matlab.

This, in addition to the fact that Matlab functions can't be called from
the command line with input variables, makes the functions unusable as is
from the command line. We can remedy this by compiling the functions. Use
the command "mcc -m biometric_system.m" to compile the Matlab function.
The created executable will be able to be run with input parameters from
the command line.

Running the executable from the command line using SIFT or SURF features
will fail without an additional step. Run the executable once and notice
the error message. In it you will find a path to a location of a missing
file. Copy the Features/Keypoint/libvl.so (for Unix) or 
Features/Keypoint/vl.dll (for Windows) file to the path from the error
message. For some reason the operating system will not copy these files
on its own.



Section 5.0   -   Features To Add
---------------------------------
- Function to combine Matlab .fig files
- Add functionality to do filtering like DCT or FFT on image preprocessing



