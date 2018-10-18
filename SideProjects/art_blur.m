function results = art_blur(input)
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
    
	fprintf('###Beginning %s\n', experiment.id);
	
    % extract gallery features
    t1 = tic;
    fprintf('##Extracting gallery features\n');
    fprintf('#From images in %s\n', experiment.input.gallery);
    fprintf('#With training set in %s\n', experiment.input.training);
	experiment.variables.quality.blur_size = 5;
	experiment.variables.quality.blur_sigma = 1;
	experiment.variables.quality.blur_num = 0;
    gallery_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.gallery, experiment.input.datadir, experiment.output.gallery_mat, experiment.output.resultsdir);
	s=whos('gallery_features');
	fprintf('Produced gallery feature matrix [%d features x %d images] [%6.2f MB]\n', size(gallery_features,1), size(gallery_features,2), s.bytes/1000000);
    toc(t1)
	fprintf('\n');
	
	blur_values = [0, 1, 3, 5, 10, 15, 20, 30, 40];
	for index = 1:9
    
		% extract probe features
		t1 = tic;
		fprintf('##Extracting probe features\n');
		fprintf('#From images in %s\n', experiment.input.probe);
		fprintf('#With training set in %s\n', experiment.input.training);
		experiment.variables.quality.blur_num = blur_values(index);
		probe_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.probe, experiment.input.datadir, experiment.output.probe_mat, experiment.output.resultsdir);
		s=whos('probe_features');
		fprintf('Produced probe feature matrix [%d features x %d images] [%6.2f MB]\n', size(probe_features,1), size(probe_features,2), s.bytes/1000000);
		toc(t1)
		fprintf('\n');
		
		% compute distance matrix
		t1 = tic;
		fprintf('##Computing distance matrix with %s\n', experiment.variables.distance);
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		s=whos('distances');
		fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
		toc(t1)
		fprintf('\n');
		clearvars probe_features;
		
		% evaluate the experiment
		t1 = tic;
		fprintf('##Evaluating the experiment\n');
		
		fprintf('2. Split true and false scores\n');
		t2 = tic;
		[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		results(index).num_ts = numel(true_scores);
		results(index).num_fs = numel(false_scores);
		clearvars distances;
		toc(t2)
		
		% compute ROC
		fprintf('3. Computing ROC\n');
		t2 = tic;
		[results(index).ROC_ver_rate, results(index).ROC_miss_rate, results(index).ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		toc(t2)
	
	end
	
    end_t = toc(start_t);
    
    % print out results
	fprintf('\nExperimental Results: %s\n', experiment.id);
    fprintf('=============================================================\n');
    fprintf('Verification/Authentication experiments:\n');
    fprintf('EER:\n');
	fprintf('%7.4f%%  ', results(1).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(2).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(3).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(4).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(5).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(6).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(7).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%  ', results(8).ROC_rates_and_threshs.EER_er*100);
	fprintf('%7.4f%%\n', results(9).ROC_rates_and_threshs.EER_er*100);
    fprintf('VR at 0.1%% FAR:\n');
    fprintf('%7.4f%%  ', results(1).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(2).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(3).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(4).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(5).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(6).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(7).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%  ', results(8).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f%%\n', results(9).ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('Number of true scores\n');
	fprintf('%d  ', results(1).num_ts);
	fprintf('%d  ', results(2).num_ts);
	fprintf('%d  ', results(3).num_ts);
	fprintf('%d  ', results(4).num_ts);
	fprintf('%d  ', results(5).num_ts);
	fprintf('%d  ', results(6).num_ts);
	fprintf('%d  ', results(7).num_ts);
	fprintf('%d  ', results(8).num_ts);
	fprintf('%d\n', results(9).num_ts);
	fprintf('Number of false scores\n');
	fprintf('%d  ', results(1).num_fs);
	fprintf('%d  ', results(2).num_fs);
	fprintf('%d  ', results(3).num_fs);
	fprintf('%d  ', results(4).num_fs);
	fprintf('%d  ', results(5).num_fs);
	fprintf('%d  ', results(6).num_fs);
	fprintf('%d  ', results(7).num_fs);
	fprintf('%d  ', results(8).num_fs);
	fprintf('%d\n', results(9).num_fs);
    fprintf('=============================================================\n');
    fprintf('Elapsed time: %f seconds\n', end_t);
end
