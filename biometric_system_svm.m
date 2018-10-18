function results = biometric_system_svm(input)
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
    
    % create probe and gallery lists
    fid = fopen(experiment.input.probe_id);
    pids = textscan(fid, '%d');
    fclose(fid);
    pids = double(pids{1});
    
    fid = fopen(experiment.input.gallery_id);
    gids = textscan(fid, '%d');
    fclose(fid);
    gids = double(gids{1});
    pnum = size(pids,1);
    gnum = size(gids,1);

    t1 = tic;
    
    
        fprintf('1. Train SVM Classifier with %d images\n',gnum);
        t2 = tic;
        % load or create the SVM model
        if ~isempty(experiment.output.distance_mat) && exist([experiment.output.resultsdir experiment.id 'model.mat'],'file')
            model = load([experiment.output.resultsdir experiment.id 'model.mat'], 'model');
            model = model.model;
        else
            model = svmtrain(gids, gallery_features', '-q -s 0 -t 0');
            if ~isempty(experiment.output.distance_mat)
                if(~exist(experiment.output.resultsdir, 'dir'))
                   mkdir(experiment.output.resultsdir);
                end
                save([experiment.output.resultsdir experiment.id 'model.mat'], 'model', '-v7.3');
            end
        end
        clearvars gallery_features;
        s=whos('model');
        fprintf('Produced SVM model [%6.2f MB]\n', s.bytes/1000000);
        toc(t2)
        
        fprintf('2. Predict Classification with SVM\n');
        t2 = tic;
        % predict classification with SVM
        if ~isempty(experiment.output.distance_mat) && exist([experiment.output.resultsdir experiment.id 'predicted.mat'],'file')
            predicted = load([experiment.output.resultsdir experiment.id 'predicted.mat'], 'predicted');
            predicted = predicted.predicted;
        else
            predicted = zeros(pnum,1);
            for k=1:pnum
                [predict_label, ~, ~] = svmpredict(pids(k), probe_features(:,k)', model);
                if predict_label == pids(k)
                    predicted(k) = 1;
                else
                    predicted(k) = 0;
                end
                printing(experiment.variables.printing, sprintf('predict label [%07d] actual [%07d]', predict_label, pids(k)), k, pnum);
            end
            if ~isempty(experiment.output.distance_mat)
                if(~exist(experiment.output.resultsdir, 'dir'))
                   mkdir(experiment.output.resultsdir);
                end
                save([experiment.output.resultsdir experiment.id 'predicted.mat'], 'predicted', '-v7.3');
            end
        end
        clearvars probe_features model;
        fprintf('\n');
    
    
    
    
    toc(t1)
    fprintf('\n');
    clearvars gallery_features probe_features;
    
    end_t = toc(start_t);

    % print out results
    fprintf('\nExperimental Results: %s\n', experiment.id);
    fprintf('=============================================================\n');
    fprintf('SVM Accuracy: %7.4f%%\n', (sum(predicted)/length(predicted))*100);
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
