function results = demographic_classifier(input)
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
    
    % extract training features
    t1 = tic;
    fprintf('##Extracting training features\n');
    fprintf('#From images in %s\n', experiment.input.training);
    training_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.training, experiment.input.datadir, experiment.output.training_mat, experiment.output.resultsdir);
    s=whos('training_features');
    fprintf('Produced experiment feature matrix [%d features x %d images] [%6.2f MB]\n', size(training_features,1), size(training_features,2), s.bytes/1000000);
    toc(t1)
    fprintf('\n');
    
    % split training lists up by class
    [tids uids] = getClassIDs(experiment);
    fprintf('There are %d unique classes.\n', numel(uids));

    
    class_features = cell(numel(uids),1);
    class_num = zeros(numel(uids),1);
    for i=1:numel(uids)
        class_features{i,1} = training_features(:,tids==uids(i));
        class_num(i) = sum(tids==uids(i));
        fprintf('%d,', class_num(i));
    end
    fprintf('\n');
    clearvars training_features;
    save([experiment.output.resultsdir experiment.id 'temp.mat'], 'class_features', '-v7.3');
    fprintf('Finished splitting features by class.\n');
    
    % experimental loop
    fprintf('Starting experimental loop.\n');
    for i=1:5
        % split list into fifths
        % compute average from 80% of list
        % build test list from remaining 20%
        % compute NN for each
        % compute accuracy
        
        test_set = [];
        test_ids = [];
        train_set = [];
        train_ids = [];
        class_features = load([experiment.output.resultsdir experiment.id 'temp.mat'], 'class_features');
        class_features = class_features.class_features;
        for j=1:numel(uids)
            itrain = 1:class_num(j);
            itest = (i-1)*floor(class_num(j)/5)+1:i*floor(class_num(j)/5);
            itrain(itest) = [];
            
            %matrix = class_features{j};
            train_avg(:,j) = mean(class_features{j}(:,itrain),2);
            train_set = [train_set class_features{j}(:,itrain)];
            train_ids = [train_ids; ones(length(itrain),1)*j];
            test_set = [test_set class_features{j}(:,itest)];
            test_ids = [test_ids; ones(length(itest),1)*j];
        end
        clearvars class_features;
        fprintf('-Finished splitting train and test set.\n');
        
        if strcmp(experiment.variables.distance,'libsvm')
            fprintf('1. Train SVM Classifier with %d images\n',size(train_ids,1));
            t2 = tic;
            model = svmtrain(train_ids, train_set', '-q -s 0 -t 0');
            toc(t2)
            fprintf('2. Predict Classification with SVM\n');
            t2 = tic;
            pred=zeros(size(test_ids,1),1);
            dec_values=zeros(size(test_ids,1),model.nr_class*(model.nr_class-1)/2);
            for k=1:size(test_ids,1)
                [pred(k,1), ~, dec_values(k,:)] = svmpredict(test_ids(k), test_set(:,k)', model);
            end
            toc(t2)
            fprintf('3. build votes array\n');
            t2 = tic;
            distances = zeros(size(test_ids,1),model.nr_class);
            index = 1;
            for k=1:size(test_ids,1)
                index=1;
                for t=1:model.nr_class
                    for j=t+1:model.nr_class
                        if dec_values(k,index) > 0
                            distances(k,t)=distances(k,t)+1;
                        else
                            distances(k,j)=distances(k,j)+1;
                        end
                        index=index+1;
                    end 
                end
            end
            toc(t2)
            save([experiment.output.resultsdir experiment.id 'distance' num2str(i) '.mat'], 'distances');
        
            %experiment.variables.pids = test_ids;
            %experiment.variables.gids = train_ids;
            %distances = compute_distances(train_set, test_set, experiment.id, experiment.variables, experiment.input, experiment.output);
            %[~,pred] = max(distances,[],2);
        elseif strcmp(experiment.variables.distance,'bayesian')
            O1 = NaiveBayes.fit(train_set',train_ids,'Distribution','normal');
            pred = O1.predict(test_set');			
        else
            distances = compute_distances(train_avg, test_set, experiment.id, experiment.variables, experiment.input, experiment.output);
            [~,pred] = min(distances,[],2);
        end
        fprintf('-Finished classifying.\n');
        
        accuracy(i) = sum(test_ids == pred)/length(test_ids);
        prediction(:,i) = pred;
        truth(:,i) = test_ids;
        for j=1:numel(uids)
            test1 = test_ids(test_ids == j);
            pred1 = pred(test_ids == j);
            class_accuracy(j,i) = sum(test1 == pred1)/length(test1);
            for k=1:numel(uids)
                class_inaccuracy(k,j,i) = sum(k == pred1)/length(test1);
            end
        end
        fprintf('-Finished computing results.\n');
        
        fprintf('Completed loop %d of %d.\n',i,5);
    end
    
    results.accuracy = accuracy;
    results.class_accuracy = class_accuracy;
    results.class_inaccuracy = class_inaccuracy;
    results.mean_accuracy = mean(accuracy);
    results.mean_class_accuracy = mean(class_accuracy,2);
    results.mean_class_inaccuracy = mean(class_inaccuracy,3);
    results.prediction = prediction;
    results.truth = truth;
    delete([experiment.output.resultsdir experiment.id 'temp.mat']);
    
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
    fprintf('Classification experiments:\n');
    fprintf('Mean Accuracy: %7.4f%%\n', results.mean_accuracy*100);
    for i=1:length(results.mean_class_accuracy)
        fprintf('Mean Accuracy for Class %d: %7.4f%%', i, results.mean_class_accuracy(i)*100);
        fprintf('   Mean Inaccuracy: ');
        for j=1:size(results.mean_class_inaccuracy,1)
            if j ~= i
                fprintf('%7.2f%%(%d) ', results.mean_class_inaccuracy(j,i)*100, j);
            end
        end
        fprintf('\n');
    end
    fprintf('=============================================================\n');
    fprintf('Elapsed time: %f seconds\n', end_t);
    
    fprintf('\n');
    for i=1:size(results.mean_class_inaccuracy,2)
        for j=1:size(results.mean_class_inaccuracy,1)
            fprintf('%7.4f,', results.mean_class_inaccuracy(j,i)*100);
        end
        fprintf('\n');
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

function [tids uids] = getClassIDs(experiment)
    %%% Scan in image filenames %%%
    fid = fopen(experiment.input.training);
    training_names = textscan(fid, '%s');
    fclose(fid);
    training_names = training_names{1};
    tnum = size(training_names,1);
    
    %%% find id for each image %%%
    assert(isempty(experiment.input.training_id) ~= 1, 'Must provide class list');
    fid = fopen(experiment.input.training_id);
    tids = textscan(fid, '%d');
    fclose(fid);
    tids = tids{1};
    assert(tnum == size(tids,1), 'IDs does not match imagelists: training count mismatch');
    uids = unique(tids);
end
