% this file is for running the experiments that go in the structure chapter of my dissertation
% it displays the rank-1 recognition rate of experiments comparing features from the different structural areas of the periocular region


function results = diss_structure_feature_variance(input)
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
	for i = 1:7
		% get subset of features matricies
		gallery_features =  gf1(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
		probe_features =  pf1(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
		%fprintf('### gallery_features ###\n');
		%disp(mean(gallery_features,2)')
		%fprintf('### gallery_features ###\n');
		
		gf_mat = zeros(flength, size(gallery_features, 2));
		pf_mat = zeros(flength, size(probe_features, 2));
		for x = 1:(loc(i,2) - loc(i,1) + 1)
			gf_mat = gf_mat + gallery_features( ((x-1)*flength)+1:x*flength , : );
			pf_mat = pf_mat + probe_features( ((x-1)*flength)+1:x*flength , : );
		end
		gf_mat = gf_mat / (loc(i,2) - loc(i,1) + 1); % [59, 466] for lbp frgc
		pf_mat = pf_mat / (loc(i,2) - loc(i,1) + 1);
		
		dprime = zeros(flength,1);
		for y = 1:flength
			ford_gf = gf_mat(y,:);
			ford_pf = pf_mat(y,:)';
			ford = bsxfun(@minus, ford_gf, ford_pf); % [466, 466] for frgc
			ford = abs(ford);
			%disp(size(ford))
			%disp(ford_gf(1:5))
			%disp(ford_pf(1:5))
			%disp(ford(1:5,1:5))
			
			true_scores = ford(1:length(ford)+1:numel(ford));
			false_scores = ford(setdiff(1:numel(ford), 1:length(ford)+1:numel(ford)));
			
			%disp(size(true_scores))
			%disp(size(false_scores))
			
			dprime(y) = abs(mean(true_scores) - mean(false_scores))/sqrt((std(true_scores)^2+std(false_scores)^2)/2);
		end
		
		gf_mat = (gf_mat + pf_mat) / 2;
		avg_f = mean(gf_mat,2); % [59, 1] for lbp frgc
		std_f = std(gf_mat,0,2); % [59, 1] for lbp frgc
		
		for y = 1:length(avg_f)
			fprintf('%f,', avg_f(y));
		end
		fprintf('\n');
		for y = 1:length(std_f)
			fprintf('%f,', std_f(y));
		end
		fprintf('\n');
		for y = 1:length(dprime)
			fprintf('%f,', dprime(y));
		end
		fprintf('\n\n');
		
	end

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
