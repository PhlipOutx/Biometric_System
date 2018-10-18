function results = twoset_fusion(input)
%input
% input(struct):	Matlab structure containing fields that hold variables
%                   used in a fusion of two experiments. See Section 2.0
%					for list of required fields.
%					
% input(string):	Filepath to xml file that will be parsed into a struct.
%
%ouput
% results:			See Section 3.0 for list of outputs.


    start_t = tic;	
    if isstruct(input)
        experiment = input;
    else
        % parse xml experiment file
        experiment = xml_read(input);
    end
    
    if strcmp(experiment.variables.fusion,'feature')
        fprintf('##Loading gallery features\n');
        % load gallery features from id1
        assert(exist([experiment.output.resultsdir experiment.variables.id1 experiment.output.gallery_mat],'file')==2, 'gallery features from id1 do not exist');
        gf1 = load([experiment.output.resultsdir experiment.variables.id1 experiment.output.gallery_mat], 'features');
        gf1 = gf1.features;
        % load gallery features from id2
        assert(exist([experiment.output.resultsdir experiment.variables.id2 experiment.output.gallery_mat],'file')==2, 'gallery features from id2 do not exist');
        gf2 = load([experiment.output.resultsdir experiment.variables.id2 experiment.output.gallery_mat], 'features');
        gf2 = gf2.features;

        gallery_features = cat(1, gf1, gf2);
        clearvars gf1 gf2;		
        
        fprintf('##Loading probe features\n');
        % load probe features from id1
        assert(exist([experiment.output.resultsdir experiment.variables.id1 experiment.output.probe_mat],'file')==2, 'probe features from id1 do not exist');
        pf1 = load([experiment.output.resultsdir experiment.variables.id1 experiment.output.probe_mat], 'features');
        pf1 = pf1.features;
        % load probe features from id2
        assert(exist([experiment.output.resultsdir experiment.variables.id2 experiment.output.probe_mat],'file')==2, 'probe features from id2 do not exist');
        pf2 = load([experiment.output.resultsdir experiment.variables.id2 experiment.output.probe_mat], 'features');
        pf2 = pf2.features;
        
        probe_features = cat(1, pf1, pf2);
        clearvars pf1 pf2;
        
        % compute distance matrix
        t1 = tic;
        fprintf('##Computing distance matrix with %s\n', experiment.variables.distance);
        distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
        s=whos('distances');
        fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
        toc(t1)
        fprintf('\n');
        clearvars gallery_features probe_features;
        
        % save distances matrix
        if ~isempty(experiment.output.distance_mat)
            if(~exist(experiment.output.resultsdir, 'dir'))
               mkdir(experiment.output.resultsdir);
            end
            save([experiment.output.resultsdir experiment.id experiment.output.distance_mat], 'distances', '-v7.3');
        end
        
    elseif strcmp(experiment.variables.fusion,'score')
        fprintf('##Loading distances from id1\n');
        % load gallery features from id1
        assert(exist([experiment.output.resultsdir experiment.variables.id1 experiment.output.distance_mat],'file')==2, 'distances from id1 do not exist');
        d1 = load([experiment.output.resultsdir experiment.variables.id1 experiment.output.distance_mat], 'distances');
        d1 = d1.distances;
        
        fprintf('##Loading distances from id2\n');
        % load gallery features from id2
        assert(exist([experiment.output.resultsdir experiment.variables.id2 experiment.output.distance_mat],'file')==2, 'distances from id2 do not exist');
        d2 = load([experiment.output.resultsdir experiment.variables.id2 experiment.output.distance_mat], 'distances');
        d2 = d2.distances;

        distances = d1+d2;
        clearvars d1 d2;
        
        % save distances matrix
        if ~isempty(experiment.output.distance_mat)
            if(~exist(experiment.output.resultsdir, 'dir'))
               mkdir(experiment.output.resultsdir);
            end
            save([experiment.output.resultsdir experiment.id experiment.output.distance_mat], 'distances', '-v7.3');
        end
    end
    
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
    fprintf('=============================================================\n');
    fprintf('Elapsed time: %f seconds\n', end_t);
end



