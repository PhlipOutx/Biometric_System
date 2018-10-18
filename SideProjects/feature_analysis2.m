function results = feature_analysis2(input)
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
	
	%features = load([experiment.output.resultsdir experiment.LBPriu2 experiment.output.gallery_mat], 'features');	gf1 = features.features;
	%features = load([experiment.output.resultsdir experiment.LBPri experiment.output.gallery_mat], 'features');	gf2 = features.features;
	features = load([experiment.output.resultsdir experiment.LBPu2 experiment.output.gallery_mat], 'features');	gf3 = features.features;
	features = load([experiment.output.resultsdir experiment.LBPall experiment.output.gallery_mat], 'features');	gf4 = features.features;
	%features = load([experiment.output.resultsdir experiment.HOG8 experiment.output.gallery_mat], 'features');	gf5 = features.features;
	features = load([experiment.output.resultsdir experiment.HOG12 experiment.output.gallery_mat], 'features');	gf6 = features.features;
	features = load([experiment.output.resultsdir experiment.LPQ experiment.output.gallery_mat], 'features');	gf7 = features.features;
	
	%either = load([experiment.output.resultsdir experiment.LBPriu2 'fa.mat'], 'either');	e1 = either.either;
	%either = load([experiment.output.resultsdir experiment.LBPri 'fa.mat'], 'either');	e2 = either.either;
	either = load([experiment.output.resultsdir experiment.LBPu2 'fa.mat'], 'either');	e3 = either.either;
	either = load([experiment.output.resultsdir experiment.LBPall 'fa.mat'], 'either');	e4 = either.either;
	%either = load([experiment.output.resultsdir experiment.HOG8 'fa.mat'], 'either');	e5 = either.either;
	either = load([experiment.output.resultsdir experiment.HOG12 'fa.mat'], 'either');	e6 = either.either;
	either = load([experiment.output.resultsdir experiment.LPQ 'fa.mat'], 'either');	e7 = either.either;

    toc(t1)
	fprintf('\n');
    
    % extract probe features
    t1 = tic;
    fprintf('##Extracting probe features\n');
    fprintf('#From images in %s\n', experiment.input.probe);
    fprintf('#With training set in %s\n', experiment.input.training);

	%features = load([experiment.output.resultsdir experiment.LBPriu2 experiment.output.probe_mat], 'features');	pf1 = features.features;
	%features = load([experiment.output.resultsdir experiment.LBPri experiment.output.probe_mat], 'features');	pf2 = features.features;
	features = load([experiment.output.resultsdir experiment.LBPu2 experiment.output.probe_mat], 'features');	pf3 = features.features;
	features = load([experiment.output.resultsdir experiment.LBPall experiment.output.probe_mat], 'features');	pf4 = features.features;
	%features = load([experiment.output.resultsdir experiment.HOG8 experiment.output.probe_mat], 'features');	pf5 = features.features;
	features = load([experiment.output.resultsdir experiment.HOG12 experiment.output.probe_mat], 'features');	pf6 = features.features;
	features = load([experiment.output.resultsdir experiment.LPQ experiment.output.probe_mat], 'features');	pf7 = features.features;
	
    toc(t1)
	fprintf('\n');
	
	flength = [10 36 59 256 8 12 256];
	if strcmp(experiment.variables.patches.type,'miller');
		loc = [1,11;12,22;23,26;27,30;31,48;49,63;64,84];
	elseif strcmp(experiment.variables.patches.type,'miller2');
		loc = [1,7;8,12;13,16;17,20;21,28;29,34;35,74];
	end
	
	% correlation based
    printing(experiment.variables.printing, 'correlation based', 1, 49);
    r1m = zeros(7,7);
	for i = 1:7
		for j = [3 4 6 7]
			if j == 1
				gallery_features =  gf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 2
				gallery_features =  gf2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 3
				gallery_features =  gf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 4
				gallery_features =  gf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 5
				gallery_features =  gf5(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf5(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e5(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 6
				gallery_features =  gf6(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf6(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e6(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 7
				gallery_features =  gf7(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf7(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e7(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			end
			
			r1=zeros(size(gallery_features,1),1);
			for k=1:size(gallery_features,1)
				[r1(k),~,~,~,~,~,~] = ICC([gallery_features(k,:);probe_features(k,:)]', 'A-1');
            end
            temp = ~isnan(r1);
            r1 = r1(temp);
            
            features = [gallery_features,probe_features];
            r2 = corrcoef(features');
            r2 = r2(temp,temp);
            r2 = ((sum(r2,2)-1)/size(r2,1));
            
            r1m(i,j) = mean(abs(r1)./sqrt(abs(r2)));
            printing(experiment.variables.printing, 'correlation based', (i-1)*7+j, 49);
		end
    end
    clearvars features temp r1 r2;
	fprintf('\nCorrelation based scores (higher better)\n');
	disp(r1m);
	
	
	% performance based
    printing(experiment.variables.printing, 'performance based', 1, 49);
	vrs = zeros(7,7);
	for i = 1:7
		for j = [3 4 6 7]
			if j == 1
				gallery_features =  gf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 2
				gallery_features =  gf2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 3
				gallery_features =  gf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 4
				gallery_features =  gf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 5
				gallery_features =  gf5(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf5(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e5(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 6
				gallery_features =  gf6(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf6(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e6(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			elseif j == 7
				gallery_features =  gf7(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				probe_features =  pf7(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
				ef =  e7(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j));
				gallery_features = gallery_features(ef,:);
				probe_features = probe_features(ef,:);
			end
			distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
			CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
			vrs(i,j) = CMC_rec_rates(1)*100;
			%[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
			%[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
			%vrs(i,j) = ROC_rates_and_threshs.VER_01FAR_ver*100;
			printing(experiment.variables.printing, 'performance based', (i-1)*7+j, 49);
		end
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
