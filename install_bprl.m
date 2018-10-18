%% Get current directory and add all subdirectories to path
function install_bprl()
root = pwd;
addpath(root);
addpath(fullfile(root,'Evaluation'));
addpath(fullfile(root,'Experiments'));
addpath(fullfile(root,'Features'));
addpath(fullfile(root,'Features','Keypoint'));
addpath(fullfile(root,'Features','Local'));
addpath(fullfile(root,'Features','Subspace'));
addpath(fullfile(root,'Features','Filter'));
addpath(fullfile(root,'Plotting'));
addpath(fullfile(root,'Similarity'));
addpath(fullfile(root,'Util'));
addpath(fullfile(root,'SideProjects'));
fprintf('Finished installing the BPRL System\n');

%savepath
