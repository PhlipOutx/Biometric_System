function results = experiment_for_every_patch(input)
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
    
    experiment.output.distance_mat = '';
    experiment.output.results_mat = '';
    
    fprintf('###Beginning %s\n', experiment.id);
    
    %%% setup variables %%%
    if strcmp(experiment.variables.feature.name,'LBP') &&  strcmp(experiment.variables.feature.type,'all')
        flength = 256;
    elseif strcmp(experiment.variables.feature.name,'LBP')
        mapping = getmapping(experiment.variables.feature.samples, experiment.variables.feature.type);
        flength = mapping.num;
    elseif strcmp(experiment.variables.feature.name,'HOG')
        flength = experiment.variables.feature.bins;
    elseif strcmp(experiment.variables.feature.name,'LPQ')
        flength = 256;
    elseif strcmp(experiment.variables.feature.name,'WLD')
        flength = experiment.variables.feature.num_orientation*experiment.variables.feature.num_excitation;
    elseif strcmp(experiment.variables.feature.name,'LCH')
        flength = 16;
    end
    patches = getpatches(experiment.variables);

    % allocate space
    patch_results = zeros(patches.num,4);
    fprintf('#Finding individual results for %d patches\n', patches.num);
    
    
    if exist([experiment.output.resultsdir 'temp' experiment.id 'efep.mat'],'file')
        temp = load([experiment.output.resultsdir 'temp' experiment.id 'efep.mat'], 'cur_p', 'patch_results');
        s_cur_p = temp.cur_p+1;
        patch_results = temp.patch_results;
        clearvars temp;
    else
        s_cur_p = 0;
    end
    for cur_p = s_cur_p:patches.num-1
        fprintf('#Starting patch %d of %d\n', cur_p+1, patches.num);
        
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
        
        gf_temp = gallery_features((cur_p*flength)+1:((cur_p+1)*flength),:);
        
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
    
        pf_temp = probe_features((cur_p*flength)+1:((cur_p+1)*flength),:);
        clearvars gallery_features probe_features;
        
        % compute distance matrix
        t1 = tic;
        fprintf('##Computing distance matrix with %s\n', experiment.variables.distance);
        distances = compute_distances(gf_temp, pf_temp, experiment.id, experiment.variables, experiment.input, experiment.output);
        s=whos('distances');
        fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
        toc(t1)
        fprintf('\n');
        clearvars gf_temp pf_temp;
        
        % evaluate the experiment
        t1 = tic;
        fprintf('##Evaluating the experiment\n');
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
        
        toc(t1)
        fprintf('\n');
        
        patch_results(cur_p+1, 1) = results.CMC_rec_rates(1)*100;
        patch_results(cur_p+1, 2) = results.ROC_rates_and_threshs.EER_er*100;
        patch_results(cur_p+1, 3) = results.ROC_rates_and_threshs.VER_01FAR_ver*100;
        patch_results(cur_p+1, 4) = results.ROC_rates_and_threshs.dprime;
        fprintf('Rank1:%7.4f  EER:%7.4f  VR:%7.4f  Dp:%7.4f\n', patch_results(cur_p+1, 1), patch_results(cur_p+1, 2), patch_results(cur_p+1, 3), patch_results(cur_p+1, 4));
        clearvars results true_scores false_scores;
        save([experiment.output.resultsdir 'temp' experiment.id 'efep.mat'], 'cur_p', 'patch_results', '-v7.3');
    
    end
    end_t = toc(start_t);
    
    % print out results
    fprintf('\nExperimental Results: %s\n', experiment.id);
    fprintf('=============================================================\n');
    fprintf('Rank1:\n');
    for i=1:patches.num
        fprintf('%7.4f,', patch_results(i,1));
    end
    fprintf('\nEER:\n');
    for i=1:patches.num
        fprintf('%7.4f,', patch_results(i,2));
    end
    fprintf('\nVR at 0.1 FAR:\n');
    for i=1:patches.num
        fprintf('%7.4f,', patch_results(i,3));
    end
    fprintf('\nDp:\n');
    for i=1:patches.num
        fprintf('%7.4f,', patch_results(i,4));
    end
    fprintf('\n');
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
