function results = dissexp12(input)
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
    
    % experiment.input.probe = '\\skynet\users\pemille\imagelists\2012_03_Dissertation\frgc_face_exp1_neutral.srt';
    % experiment.input.gallery = '\\skynet\users\pemille\imagelists\2012_03_Dissertation\frgc_face_exp1_neutral.srt';
    % experiment.input.datadir = '\\skynet\static\Processed_Images\FRGC\FRGC_ALIGNEDFACE\';
    % experiment.output.resultsdir = '\\skynet\users\pemille\results\Dissertation\';
    % experiment.variables.printing = 1;
    
	fprintf('###Beginning %s\n', experiment.id);
	fprintf('Image size: %d X %d\n', experiment.variables.data.max_y, experiment.variables.data.max_x);
	fprintf('Patch size: %d X %d\n', experiment.variables.patches.size_y, experiment.variables.patches.size_x);
	fprintf('Patch type: %s\n', experiment.variables.patches.type);
	fprintf('Distance:   %s\n', experiment.variables.distance);
	fprintf('Feature:    %s\n', experiment.variables.feature.name);

	% ids
	% w1id = 'diss41000001';
	% w3id = 'diss41000004';
	% w5id = 'diss41000007';
	% w7id = 'diss41000010';
	% w9id = 'diss41000013';
	
	% compute distance matrix
	weights = [1.0 5.0 1.0 0.0 0.0 0.0];
    t1 = tic;
    fprintf('##Computing distance matrix with %s\n', experiment.variables.distance);
	results = load([experiment.output.resultsdir experiment.w1id 'distance.mat'], 'distances'); 
		distances = (weights(1)*results.distances);
	results = load([experiment.output.resultsdir experiment.w3id 'distance.mat'], 'distances');
		distances = distances + (weights(2)*results.distances);
	results = load([experiment.output.resultsdir experiment.w5id 'distance.mat'], 'distances');
		distances = distances + (weights(3)*results.distances);
	results = load([experiment.output.resultsdir experiment.w7id 'distance.mat'], 'distances');
		distances = distances + (weights(4)*results.distances);
	%results = load([experiment.output.resultsdir experiment.w9id 'distance.mat'], 'distances'); 
	%	distances = distances + (weights(5)*results.distances);
	%results = load([experiment.output.resultsdir experiment.w11id 'distance.mat'], 'distances'); 
	%	distances = distances + (weights(6)*results.distances);
	clearvars results;
	s=whos('distances');
	fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
	if ~isempty(experiment.output.distance_mat)
		save([experiment.output.resultsdir experiment.id experiment.output.distance_mat], 'distances', '-v7.3');
	end
    toc(t1)
	fprintf('\n');
	
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
	fprintf('%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n',results.CMC_rec_rates(1)*100,results.ROC_rates_and_threshs.EER_er*100,results.ROC_rates_and_threshs.VER_01FAR_ver*100,results.ROC_rates_and_threshs.VER_001FAR_frr*100,results.ROC_rates_and_threshs.VER_01FAR_frr*100,results.ROC_rates_and_threshs.VER_1FAR_frr*100);
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
