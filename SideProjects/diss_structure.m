% this file is for running the experiments that go in the structure chapter of my dissertation
% it displays the rank-1 recognition rate of experiments comparing features from the different structural areas of the periocular region


function results = diss_structure(input)
% input(struct):	Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%					list of required fields.
%					
% input(string):	Filepath to xml file that will be parsed into a struct.
%
% results:			See Section 3.0 for list of outputs.

	start_t = tic;	
	if isstruct(input)
		experiment = input;
	else
		% parse xml experiment file
		experiment = xml_read(input);
	end
    
    % extract gallery features
    t1 = tic;
    fprintf('##Extracting gallery features\n');
    fprintf('#From images in %s\n', experiment.input.gallery);
    fprintf('#With training set in %s\n', experiment.input.training);
	
	features = load([experiment.output.resultsdir experiment.id experiment.output.gallery_mat], 'features');	
	gf1 = features.features;
	fprintf('#Loaded gallery features [%d X %d]\n', size(gf1,1), size(gf1,2));

    toc(t1)
	fprintf('\n');
    
    % extract probe features
    t1 = tic;
    fprintf('##Extracting probe features\n');
    fprintf('#From images in %s\n', experiment.input.probe);
    fprintf('#With training set in %s\n', experiment.input.training);

	features = load([experiment.output.resultsdir experiment.id experiment.output.probe_mat], 'features');
	pf1 = features.features;
	fprintf('#Loaded probe features [%d X %d]\n', size(pf1,1), size(pf1,2));
	
    toc(t1)
	fprintf('\n');
	
	loc = [1,7;8,12;13,16;17,20;21,28;29,34;35,74];
	
	%get flength from feature name
	flength = 0;
	if strcmp(experiment.variables.feature.name,'LBP') == 1
		flength = 59;
	elseif strcmp(experiment.variables.feature.name,'HOG') == 1
		flength = experiment.variables.feature.bins;
	elseif strcmp(experiment.variables.feature.name,'LPQ') == 1
		flength = 256;
	end
	
	% performance based
	vrs = zeros(7,1);
	distances = cell(7,1);
	for i = 1:7
		% get subset of features matricies
		gallery_features =  gf1(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
		probe_features =  pf1(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);

		d1 = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
		vrs(i,1) = CMC_rec_rates(1)*100;
		distances{i} = d1;
		printing(experiment.variables.printing, 'performance based', i, 7);
	end
	
	% saving this to be used with diss_multiscale_fusion.m. This can probably be commented out for any other use.
	save([experiment.output.resultsdir experiment.id 'distance.mat'], 'distances', '-v7.3');

	
    fprintf('\nPerformance based scores (higher better)\n');
	for i = 1:7
		fprintf('%7.4f,',vrs(i));
	end
	fprintf('\n');
end

function ret = checkSameList(gallery, probe)
	%%% Scan in image filenames %%%
	fid = fopen(gallery);
	glist = textscan(fid, '%s');
	fclose(fid);
	glist = glist{1};
	gnum = size(glist,1);
	
	%%% Scan in image filenames %%%
	fid = fopen(probe);
	plist = textscan(fid, '%s');
	fclose(fid);
	plist = plist{1};
	pnum = size(plist,1);
	
	if pnum ~= gnum
		ret=0;
	elseif (sum(strcmp(plist, glist)) == gnum)
		ret=1;
	else
		ret=0;
	end
end
