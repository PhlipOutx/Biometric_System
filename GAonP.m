function results = GAonP(input)
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
    
    % extract gallery features
    t1 = tic;
    fprintf('##Extracting gallery features\n');
    fprintf('#From images in %s\n', experiment.input.gallery);
    gallery_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.gallery, experiment.input.datadir, experiment.output.gallery_mat, experiment.output.resultsdir);
    s=whos('gallery_features');
    fprintf('Produced gallery feature matrix [%d features x %d images] [%6.2f MB]\n', size(gallery_features,1), size(gallery_features,2), s.bytes/1000000);
    toc(t1)
    fprintf('\n');
    
    % extract probe features
    t1 = tic;
    fprintf('##Extracting probe features\n');
    fprintf('#From images in %s\n', experiment.input.probe);
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
    
    %% initialize GA variables
    numgenes = (experiment.variables.data.max_x/experiment.variables.patches.size_x)*(experiment.variables.data.max_y/experiment.variables.patches.size_y);
    numchrom = 32;
    mutation_rate = .1;
    numiter = 50000;
    actiter = 1;
    since = 0;
    
    %%% Scan in image filenames %%%
    fid = fopen(experiment.input.probe);
    probe_names = textscan(fid, '%s');
    fclose(fid);
    probe_names = probe_names{1};
    pnum = size(probe_names,1);
    fid = fopen(experiment.input.gallery);
    gallery_names = textscan(fid, '%s');
    fclose(fid);
    gallery_names = gallery_names{1};
    gnum = size(gallery_names,1);
    
    %%% find id for each image %%%
    if isempty(experiment.input.probe_id) && isempty(experiment.input.gallery_id)
        psub = cellfun(@(x) x(1:experiment.variables.nameprefix), probe_names, 'UniformOutput', false);
        gsub = cellfun(@(x) x(1:experiment.variables.nameprefix), gallery_names, 'UniformOutput', false);

        asub=[psub;gsub];
        [avalues, aloc, aids] = unique(asub);
        
        [trash, psubloc] = ismember(psub,asub);
        pids = aids(psubloc); 
        
        [trash, gsubloc] = ismember(gsub,asub);
        gids = aids(gsubloc);
        
        clearvars avalues aloc trash psubloc gsubloc aids gsub psub;
    else
        fid = fopen(experiment.input.probe_id);
        pids = textscan(fid, '%d');
        fclose(fid);
        pids = pids{1};
        assert(pnum == size(pids,1), 'IDs does not match imagelists: Probe count mismatch');
        
        fid = fopen(experiment.input.gallery_id);
        gids = textscan(fid, '%d');
        fclose(fid);
        gids = gids{1};
        assert(gnum == size(gids,1), 'IDs does not match imagelists: Gallery count mismatch');
    end
        
    % find location of ts and fs
    nts = 0;nfs = 0;
    temp = zeros(gnum,pnum);
    for i = 1:pnum
        s1 = pids(i) == gids;
        s2 = strcmp(probe_names(i), gallery_names);
        nts = nts + sum(xor(s1,s2));
        nfs = nfs + sum(~s1);
        temp(:,i) = xor(s1,s2);
    end
    experiment.variables.ts_mat=find(temp==1);
    for i = 1:pnum
        s1 = pids(i) == gids;
        s2 = strcmp(probe_names(i), gallery_names);
        temp(:,i) = ~s1;
    end
    experiment.variables.fs_mat=find(temp==1);
    experiment.variables.gids = gids;
    experiment.variables.pids = pids;
    experiment.variables.nts = nts;
    experiment.variables.nfs = nfs;
    experiment.variables.gnum = gnum;
    experiment.variables.pnum = pnum;
    experiment.variables.probe_names = probe_names;
    experiment.variables.gallery_names = gallery_names;
    clearvars s1 s2 i nfs nts pids gids gallery_names probe_names pnum gnum temp; 
    
    
    % find score for each patch
    t1 = tic;
    fprintf('##Finding score for each patch\n');
    if ~isempty(experiment.output.pscores_mat) && exist([experiment.output.resultsdir experiment.id experiment.output.pscores_mat],'file')
        fprintf('#Already done\n');
    elseif ~isempty(experiment.output.pscores_mat)
        pscores = zeros(numgenes,1);
        for i=1:numgenes
            chromosome = zeros(numgenes,1);
            chromosome(i) = 1;
            pscores(i) = getScore(chromosome, probe_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
            printing(experiment.variables.printing, 'finding score for each patch', i, numgenes);
        end
        fprintf('\n');
        save ([experiment.output.resultsdir experiment.id experiment.output.pscores_mat], 'pscores');
    end
    toc(t1)
    
    % randomly initilize population
    t1 = tic;
    fprintf('##Randomly initilize population\n');
    if ~isempty(experiment.output.population_mat) && exist([experiment.output.resultsdir experiment.id experiment.output.population_mat],'file')
        fprintf('#Load population\n');
        load ([experiment.output.resultsdir experiment.id experiment.output.population_mat], 'population', 'scores', 'baseline', 'actiter', 'since');
    else
        population = round(rand(numgenes,numchrom));
        scores = zeros(numchrom,1);
        baseline = getScore(ones(numgenes,1), probe_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
        for i=1:numchrom
            chromosome = population(:,i);
            scores(i) = getScore(chromosome, probe_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
            printing(experiment.variables.printing, 'randomly initilizing popluation', i, numchrom);
        end
        fprintf('\n');
        if ~isempty(experiment.output.population_mat)
            save ([experiment.output.resultsdir experiment.id experiment.output.population_mat], 'population', 'scores', 'baseline', 'actiter', 'since');
        end
    end
    fprintf('actiter:   %d\n', actiter);
    fprintf('baseline:  %f\n', baseline);
    fprintf('max score: %f\n', max(scores));
    fprintf('min score: %f\n', min(scores));
    toc(t1)
    
    % run iterations
    fprintf('##Running iterations of genetic algorithm\n');
    for actiter = actiter:numiter
        t1 = tic;
        winners = tournament_selection(scores);
        child = uniform_crossover(population(:,winners(1)), population(:,winners(2)));
        child = random_mutation(child, mutation_rate);

        % see if child is already in population
        i = 1;
        same = 0;
        while i <= length(scores)
            temp = sum(child==population(:,i));
            if temp == size(population,1)
                same = 1;
                break;
            end
            i=i+1;
        end
        
        % if child is not in population
        if same == 0
            score = getScore(child, probe_features, gallery_features, experiment.id, experiment.variables, experiment.input, experiment.output);
            [population scores change] = replace_worst(population, scores, child, score);
            if change == 1
                since = 0;
            else
                since = since + 1;
            end
            time = toc(t1);
            stemp = sprintf('[%f seconds   %7.4f score   %7.4f max   %d since change] running iteration',time,score,max(scores),since);
            printing(experiment.variables.printing, stemp, actiter, numiter);
            save ([experiment.output.resultsdir experiment.id experiment.output.population_mat], 'population', 'scores', 'baseline', 'actiter', 'since');
        end
        
        if since > 1000
            fprintf('\nIt has been over 1000 iterations since a change has been made to the population. Terminating loop.\n');
            break;
        end
    end
    fprintf('\n#Iterations complete\n');
    [best_score,I] = max(scores);
    fprintf('Final Score: %7.4f\n', best_score);
    best_chromosome = population(:,I);
    disp(best_chromosome);
    

    
    
    
    
    
    
    
    
        
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
    
    % choose patches from the best chromosome
    lenfeat = size(probe_features,1)/((experiment.variables.data.max_x/experiment.variables.patches.size_x)*(experiment.variables.data.max_y/experiment.variables.patches.size_y));
    f = find(best_chromosome==1);
    chrom = [];
    for i=1:size(f,1)
        chrom=[chrom ((f(i)-1)*lenfeat)+1:f(i)*lenfeat];
    end
    probe_features=probe_features(chrom,:);
    gallery_features=gallery_features(chrom,:);
    
    
    % compute distance matrix
    t1 = tic;
    fprintf('##Computing distance matrix with %s\n', experiment.variables.distance);
    distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
    s=whos('distances');
    fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
    % save distances matrix
    save([experiment.output.resultsdir experiment.id 'distance.mat'], 'distances', '-v7.3');
    toc(t1)
    fprintf('\n');
    clearvars gallery_features probe_features;
    
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

function score = getScore(chromosome, probe_features, gallery_features, id, variables, input, output)
    lenfeat = size(probe_features,1)/((variables.data.max_x/variables.patches.size_x)*(variables.data.max_y/variables.patches.size_y));
    f = find(chromosome==1);
    chrom = [];
    for i=1:size(f,1)
        chrom=[chrom ((f(i)-1)*lenfeat)+1:f(i)*lenfeat];
    end
    probe_features=probe_features(chrom,:);
    gallery_features=gallery_features(chrom,:);
    distances = slmetric_pw(gallery_features, probe_features, variables.distance);
    % split true false scores
    true_scores = distances(variables.ts_mat);
    false_scores = distances(variables.fs_mat);
    clearvars s1 s2 i distances;
    %score = paep_roc(true_scores, false_scores, output.roc_resolution)*100;
    [~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, output.roc_resolution);
    score = ROC_rates_and_threshs.EER_er*100;
end

function winners = tournament_selection(scores)

    rand_scores = randperm(length(scores));
    rand_scores = rand_scores(1:length(scores)/2);
    while length(rand_scores) > 2
        temp = [];
        for x=1:2:length(rand_scores)
            if scores(rand_scores(x)) < scores(rand_scores(x+1))
                temp = [temp;rand_scores(x)];
            else
                temp = [temp;rand_scores(x+1)];
            end
        end
        rand_scores = temp;
    end
    winners = rand_scores;
end

function child = uniform_crossover(parent1, parent2)

    map = round(rand(length(parent1),1));
    child = zeros(length(parent1),1);
    for i=1:length(parent1)
        if (map(i) == 1)
            child(i,1) = parent1(i);
        else
            child(i,1) = parent2(i);
        end
    end
end
    
function child = random_mutation(child, mutation_rate)

    for i=1:length(child)
        if rand() < mutation_rate
            if (child(i) == 0)
                child(i) = 1;
            else
                child(i) = 0;
            end
        end
    end
end

function [population scores change] = replace_worst(population, scores, child, score)

    [worst loc] = min(scores);
    change = 0;
    i = 1;
    same = 0;
    while i <= length(scores)
        temp = sum(child==population(:,i));
        if temp == size(population,1)
            same = 1;
            break;
        end
        i=i+1;
    end
    if score < worst & same == 0
        population(:,loc) = child;
        scores(loc) = score;
        change = 1;
    end
end

