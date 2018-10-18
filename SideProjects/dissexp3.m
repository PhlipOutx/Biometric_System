function results = dissexp3(input)
% input(struct):	Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%					list of required fields.
%					
% input(string):	Filepath to xml file that will be parsed into a struct.
%
% results:			See Section 3.0 for list of outputs.



%%%%%%%%%%%%%%%%%%%%%
%%%%	unweighted score based fusion - Miller dissertation patches - Miller disseration feature extraction method for untrained classifier - on periocular using FRGC data
%%%%%%%%%%%%%%%%%%%%%


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
	
	features = load([experiment.output.resultsdir experiment.LBPu2 experiment.output.gallery_mat], 'features');		gf1 = features.features;
	features = load([experiment.output.resultsdir experiment.LBPall experiment.output.gallery_mat], 'features');	gf2 = features.features;
	features = load([experiment.output.resultsdir experiment.HOG12 experiment.output.gallery_mat], 'features');		gf3 = features.features;
	features = load([experiment.output.resultsdir experiment.LPQ experiment.output.gallery_mat], 'features');		gf4 = features.features;
	clearvars features;

    toc(t1)
	fprintf('\n');
    
    % extract probe features
    t1 = tic;
    fprintf('##Extracting probe features\n');
    fprintf('#From images in %s\n', experiment.input.probe);
    fprintf('#With training set in %s\n', experiment.input.training);

	%features = load([experiment.output.resultsdir experiment.LBPu2 experiment.output.probe_mat], 'features');	pf1 = features.features;
	%features = load([experiment.output.resultsdir experiment.LBPall experiment.output.probe_mat], 'features');	pf2 = features.features;
	%features = load([experiment.output.resultsdir experiment.HOG12 experiment.output.probe_mat], 'features');	pf3 = features.features;
	%features = load([experiment.output.resultsdir experiment.LPQ experiment.output.probe_mat], 'features');		pf4 = features.features;
	clearvars features;
	
    toc(t1)
	fprintf('\n');
	
	% extract probe features
    t1 = tic;
    fprintf('##Computing partial distance matrices\n');

	flength = [59 256 12 256];
	loc = [1,11;12,22;23,26;27,30;31,48;49,63;64,84];
	
	i=1;j=3;
	gallery_features =  gf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d1 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d1 = (d1-min(min(d1)))/(max(max(d1))-min(min(d1)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	i=2;j=1;
	gallery_features =  gf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d2 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d2 = (d2-min(min(d2)))/(max(max(d2))-min(min(d2)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	i=3;j=4;
	gallery_features =  gf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d3 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d3 = (d3-min(min(d3)))/(max(max(d3))-min(min(d3)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	i=4;j=4;
	gallery_features =  gf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d4 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d4 = (d4-min(min(d4)))/(max(max(d4))-min(min(d4)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	i=5;j=1;
	gallery_features =  gf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf1(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d5 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d5 = (d5-min(min(d5)))/(max(max(d5))-min(min(d5)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	i=6;j=2;
	gallery_features =  gf2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf2(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d6 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d6 = (d6-min(min(d6)))/(max(max(d6))-min(min(d6)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	i=7;j=4;
	gallery_features =  gf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	%probe_features =  pf4(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	d7 = compute_distances(gallery_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	d7 = (d7-min(min(d7)))/(max(max(d7))-min(min(d7)));
	fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	% i=8;j=3;
	% gallery_features =  gf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	% probe_features =  pf3(((loc(i,1)-1)*flength(j))+1:loc(i,2)*flength(j),:);
	% d8 = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	% d8 = (d8-min(min(d8)))/(max(max(d8))-min(min(d8)));
	% fprintf('Finished [i=%d][j=%d]\n',i,j);
	
	clearvars gf1 gf2 gf3 gf4 pf1 pf2 pf3 pf4 gallery_features probe_features;
	
	toc(t1)
	fprintf('\n');
	
	% compute distance matrix
    t1 = tic;
    fprintf('##Computing distance matrix with %s\n', experiment.variables.distance);
	distances = d1+d2+d3+d4+d5+d6+d7;
	s=whos('distances');
	fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
    toc(t1)
	fprintf('\n');
	clearvars d1 d2 d3 d4 d5 d6 d7;
	
	% evaluate the experiment
	t1 = tic;
	fprintf('##Evaluating the experiment\n');
	% compute CMC
	fprintf('1. Computing CMC\n');
	t2 = tic;
	results.CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
	toc(t2)
	
	fprintf('2. Split true and false scores\n');
	t2 = tic;
	[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
	clearvars distances;
	toc(t2)
	
	% compute ROC
	fprintf('3. Computing ROC\n');
	t2 = tic;
	[results.ROC_ver_rate, results.ROC_miss_rate, results.ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
	toc(t2)
	
	% print points on each curve to a file
	if experiment.output.make_files == 1
		fprintf('5. Print points to file\n');
		t2 = tic;
		file = fopen([experiment.output.resultsdir experiment.id '_det.txt'], 'w');
		for i = 1:length(results.ROC_ver_rate)
		   fprintf(file, '%f,%f\n', results.ROC_miss_rate(i), 1-results.ROC_ver_rate(i));
		end
		fclose(file);

		file = fopen([experiment.output.resultsdir experiment.id '_roc.txt'], 'w');
		for i = 1:length(results.ROC_ver_rate)
		   fprintf(file, '%f,%f\n', results.ROC_miss_rate(i), results.ROC_ver_rate(i));
		end
		fclose(file);

		file = fopen([experiment.output.resultsdir experiment.id '_cmc.txt'], 'w');
		for i = 1:length(results.CMC_rec_rates)
		   fprintf(file, '%f\n', results.CMC_rec_rates(i));
		end
		fclose(file);
		toc(t2)	
	end
	
	toc(t1)
	fprintf('\n');
	
	% plotting results
	if experiment.output.make_figures == 1
		t1 = tic;
		fprintf('##Plotting the results\n');
		
		hdet = plot_det(1-results.ROC_ver_rate, results.ROC_miss_rate, results.ROC_rates_and_threshs.EER_er, experiment.variables.feature.name, experiment.output.plot_color, experiment.output.plot_linewidth, 1);
		saveas(hdet,[experiment.output.resultsdir experiment.id '_det.fig']);
		close(hdet);
		
		hcmc = plot_cmc(results.CMC_rec_rates, experiment.variables.feature.name, experiment.output.plot_color, experiment.output.plot_linewidth, 2);
		saveas(hcmc,[experiment.output.resultsdir experiment.id '_cmc.fig']);
		close(hcmc);
		
		hroc = plot_roc(results.ROC_ver_rate, results.ROC_miss_rate, results.ROC_rates_and_threshs.EER_er, experiment.variables.feature.name, experiment.output.plot_color, experiment.output.plot_linewidth, 3);
		saveas(hroc,[experiment.output.resultsdir experiment.id '_roc.fig']);
		close(hroc);
		
		toc(t1)
		fprintf('\n');
	end
    end_t = toc(start_t);
	
	% save results matrix
	if ~isempty(experiment.output.results_mat)
		if(~exist(experiment.output.resultsdir, 'dir'))
		   mkdir(experiment.output.resultsdir);
		end
		save([experiment.output.resultsdir experiment.id experiment.output.results_mat], 'results');
	end
    
    % print out results
	fprintf('\nExperimental Results: %s\n', experiment.id);
    fprintf('=============================================================\n');
    fprintf('Identification experiments:\n');
    fprintf('Rank-1   recognition rate: %7.4f%%\n', results.CMC_rec_rates(1)*100);
	if length(results.CMC_rec_rates) > 10
		fprintf('Rank-10  recognition rate: %7.4f%%\n', results.CMC_rec_rates(10)*100);
	else
		fprintf('Rank-%d recognition rate: %7.4f%%\n', length(results.CMC_rec_rates), results.CMC_rec_rates(length(results.CMC_rec_rates))*100);
	end
	if length(results.CMC_rec_rates) > 100
		fprintf('Rank-100 recognition rate: %7.4f%%\n', results.CMC_rec_rates(100)*100);
	else
		fprintf('Rank-%d recognition rate: %7.4f%%\n', length(results.CMC_rec_rates), results.CMC_rec_rates(length(results.CMC_rec_rates))*100);
	end
    fprintf('Verification/Authentication experiments:\n');
    fprintf('EER: %7.4f%%\n', results.ROC_rates_and_threshs.EER_er*100);
	fprintf('Dprime: %7.4f\n', results.ROC_rates_and_threshs.dprime);
    fprintf('VR at 1%% FAR:    %7.4f%%\n', results.ROC_rates_and_threshs.VER_1FAR_ver*100);
    fprintf('VR at 0.1%% FAR:  %7.4f%%\n', results.ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('VR at 0.01%% FAR: %7.4f%%\n', results.ROC_rates_and_threshs.VER_001FAR_ver*100);
	fprintf('Recap:\n');
	fprintf('%7.4f\t%7.4f\t%7.4f\n',results.CMC_rec_rates(1)*100,results.ROC_rates_and_threshs.EER_er*100,results.ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('=============================================================\n');
    fprintf('Elapsed time: %f seconds\n', end_t);
	

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
