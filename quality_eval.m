function results = quality_eval(input)
% input(struct):    Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%                   list of required fields.
%					
% input(string):    Filepath to xml file that will be parsed into a struct.
%
% results:          See Section 3.0 for list of outputs.

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
    gallery_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.gallery, experiment.input.datadir, experiment.output.gallery_mat, experiment.output.resultsdir);
    s=whos('gallery_features');
    fprintf('Produced gallery feature matrix [%d features x %d images] [%6.2f MB]\n', size(gallery_features,1), size(gallery_features,2), s.bytes/1000000);
    toc(t1)
    fprintf('\n');
    
    % extract probe features
    t1 = tic;
    fprintf('##Extracting probe features\n');
    fprintf('#From images in %s\n', experiment.input.probe);
    fprintf('#With training set in %s\n', experiment.input.training);
    if checkSameList(experiment.input.gallery, experiment.input.probe) == 1
        fprintf('Lists are the same\n');
        probe_features = gallery_features;
    else
        probe_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.probe, experiment.input.datadir, experiment.output.probe_mat, experiment.output.resultsdir);
    end
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
    clearvars gallery_features probe_features;
    
    % evaluate the experiment
    fprintf('##Evaluating the experiment\n');
    
    for i=1:10
        % compute CMC
        fprintf('1. Computing CMC\n');
        t2 = tic;
        %results.CMC_rec_rates = blur_compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
        toc(t2)
        
        fprintf('2. Split true and false scores\n');
        t2 = tic;
        [true_scores false_scores range] = blur_split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output, i);
        results(i).num_ts = numel(true_scores);
        results(i).num_fs = numel(false_scores);
        results(i).range = range;
        toc(t2)
        
        % compute ROC
        fprintf('3. Computing ROC\n');
        t2 = tic;
        [results(i).ROC_ver_rate, results(i).ROC_miss_rate, results(i).ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
        toc(t2)
    end
    clearvars distances;	
    end_t = toc(start_t);
    
    % save results matrix
    if ~isempty(experiment.output.results_mat)
        if(~exist(experiment.output.resultsdir, 'dir'))
           mkdir(experiment.output.resultsdir);
        end
        save([experiment.output.resultsdir experiment.id experiment.output.results_mat], 'results');
    end
    
    % print out results
    fprintf('\nExperimental Results\n');
    fprintf('=============================================================\n');
    %fprintf('Identification experiments:\n');
    %fprintf('Rank-1   recognition rate: %7.4f%%\n', results.CMC_rec_rates(1)*100);
    %if length(results.CMC_rec_rates) > 10
    %	fprintf('Rank-10  recognition rate: %7.4f%%\n', results.CMC_rec_rates(10)*100);
    %else
    %	fprintf('Rank-%d recognition rate: %7.4f%%\n', length(results.CMC_rec_rates), results.CMC_rec_rates(length(results.CMC_rec_rates))*100);
    %end
    %if length(results.CMC_rec_rates) > 100
    %	fprintf('Rank-100 recognition rate: %7.4f%%\n', results.CMC_rec_rates(100)*100);
    %else
    %	fprintf('Rank-%d recognition rate: %7.4f%%\n', length(results.CMC_rec_rates), results.CMC_rec_rates(length(results.CMC_rec_rates))*100);
    %end
    fprintf('Verification/Authentication experiments:  1\n');
    fprintf('EER:\n');
    fprintf('%7.4f%%  ', results(1).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(2).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(3).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(4).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(5).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(6).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(7).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(8).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%  ', results(9).ROC_rates_and_threshs.EER_er*100);
    fprintf('%7.4f%%\n', results(10).ROC_rates_and_threshs.EER_er*100);
    fprintf('VR at 0.1%% FAR:\n');
    fprintf('%7.4f%%  ', results(1).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(2).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(3).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(4).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(5).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(6).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(7).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(8).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%  ', results(9).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('%7.4f%%\n', results(10).ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('Number of true scores\n');
    fprintf('%d  ', results(1).num_ts);
    fprintf('%d  ', results(2).num_ts);
    fprintf('%d  ', results(3).num_ts);
    fprintf('%d  ', results(4).num_ts);
    fprintf('%d  ', results(5).num_ts);
    fprintf('%d  ', results(6).num_ts);
    fprintf('%d  ', results(7).num_ts);
    fprintf('%d  ', results(8).num_ts);
    fprintf('%d  ', results(9).num_ts);
    fprintf('%d\n', results(10).num_ts);
    fprintf('Number of false scores\n');
    fprintf('%d  ', results(1).num_fs);
    fprintf('%d  ', results(2).num_fs);
    fprintf('%d  ', results(3).num_fs);
    fprintf('%d  ', results(4).num_fs);
    fprintf('%d  ', results(5).num_fs);
    fprintf('%d  ', results(6).num_fs);
    fprintf('%d  ', results(7).num_fs);
    fprintf('%d  ', results(8).num_fs);
    fprintf('%d  ', results(9).num_fs);
    fprintf('%d\n', results(10).num_fs);
    fprintf('Maximum value of range\n');
    fprintf('%7.4f  ', results(1).range);
    fprintf('%7.4f  ', results(2).range);
    fprintf('%7.4f  ', results(3).range);
    fprintf('%7.4f  ', results(4).range);
    fprintf('%7.4f  ', results(5).range);
    fprintf('%7.4f  ', results(6).range);
    fprintf('%7.4f  ', results(7).range);
    fprintf('%7.4f  ', results(8).range);
    fprintf('%7.4f  ', results(9).range);
    fprintf('%7.4f\n', results(10).range);
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

