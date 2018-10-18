% this file is for running the experiments that go in the structure chapter of my dissertation
% it displays the rank-1 recognition rate of experiments comparing features from the different structural areas of the periocular region
% some features are masked out based on the feature footprint for their given sub-region


function results = diss_structure_mask(input)
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
	
	features = load([experiment.output.resultsdir experiment.id experiment.output.gallery_mat], 'features');	gf1 = features.features;

    toc(t1)
	fprintf('\n');
    
    % extract probe features
    t1 = tic;
    fprintf('##Extracting probe features\n');
    fprintf('#From images in %s\n', experiment.input.probe);
    fprintf('#With training set in %s\n', experiment.input.training);

	features = load([experiment.output.resultsdir experiment.id experiment.output.probe_mat], 'features');	pf1 = features.features;
	
    toc(t1)
	fprintf('\n');
	
	loc = [1,7;8,12;13,16;17,20;21,28;29,34;35,74];
	
	%get flength from feature name
	flength = 0;
	if strcmp(experiment.variables.feature.name,'LBP') == 1
		flength = 59;
	elseif strcmp(experiment.variables.feature.name,'HOG') == 1
		flength = 12;
	elseif strcmp(experiment.variables.feature.name,'LPQ') == 1
		flength = 256;
	end
	
	% performance based
	vrs = zeros(7,1);
	for i = 1:7
		
		% get features to use
		mask = getMask(experiment,i);
	
		% get subset of features matricies
		gallery_features =  gf1(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
		probe_features =  pf1(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
		
		mask = repmat(mask,size(gallery_features,2),(loc(i,2)-loc(i,1)+1));
		mask = mask';
		
		gallery_features = gallery_features .* mask;
		probe_features = probe_features .* mask;

		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		vrs(i,1) = CMC_rec_rates(1)*100;
		printing(experiment.variables.printing, 'performance based', i, 7);
	end
    fprintf('\nPerformance based scores (higher better)\n');
	disp(vrs);
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


function mask = getMask(experiment,i)
	% get features to use
		if 2 == 1 % means
			if strcmp(experiment.variables.feature.name,'LBP') == 1
				if i == 1
					mask = [0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1];
				elseif i == 2
					mask = [0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,1,1];
				elseif i == 3
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1];
					else
						mask = [0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1];
					end
				elseif i == 4
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1];
					else
						mask = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1];
					end
				elseif i == 5
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1];
					else
						mask = [0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1];
					end
				elseif i == 6
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1];
					else
						mask = [0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1];
					end
				elseif i == 7
					mask = [0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1];
				end
			end
		elseif 2 == 1 %dp above mean
			if strcmp(experiment.variables.feature.name,'LBP') == 1
				if i == 1
					mask = [1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1];
				elseif i == 2
					mask = [1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1];
				elseif i == 3
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,1,0,1,1];
					else
						mask = [1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,1,1,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,0,0,0,1,0];
					end
				elseif i == 4
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1];
					else
						mask = [1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,1,0,0,1,1,0,1,0];
					end
				elseif i == 5
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [0,1,0,0,0,0,1,0,0,0,1,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1];
					else
						mask = [0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,1,1,1,1,0,0,1,1,1,1,1,0,0,0,0,0,1,1];
					end
				elseif i == 6
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,1,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,1,0,0,0,0,1,1,1];
					else
						mask = [1,1,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0];
					end
				elseif i == 7
					mask = [1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,1,0,1,0];
				end
			end
		elseif 1 == 1 % top 3/4 dp
			if strcmp(experiment.variables.feature.name,'LBP') == 1
				if i == 1
					mask = [1,1,0,1,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1];
				elseif i == 2
					mask = [1,1,0,1,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1];
				elseif i == 3
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,1,0,1,0,1,1,0,0,1,1,1,0,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,1,1,1,1,1];
					else
						mask = [1,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1];
					end
				elseif i == 4
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,1,0,1,0,0,1,0,0,1,1,1,0,1,0,0,0,0,1,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,1,1,1];
					else
						mask = [1,1,0,0,0,0,1,0,0,1,1,1,0,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1];
					end
				elseif i == 5
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,1,0,0,1,0,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,0,0,1,1,1,1];
					else
						mask = [1,1,0,0,1,0,1,0,0,1,1,1,0,1,1,1,0,0,0,0,1,1,1,0,1,1,0,0,0,0,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,0,0,1,1];
					end
				elseif i == 6
					if strcmp(experiment.variables.patches.side,'left') == 1
						mask = [1,1,0,0,1,0,1,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,1,1,1,1];
					else
						mask = [1,1,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,1];
					end
				elseif i == 7
					mask = [1,1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,1,0,1,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1];
				end
			end
		end
end