function results = demographic_classifier_cascade(input)
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
	[tids gids eids uids] = getClassIDs(experiment);
	fprintf('There are %d unique classes.\n', numel(uids));
	
	class_features = cell(numel(uids),1);
	class_num = zeros(numel(uids),1);
	for i=1:numel(uids)
		class_features{i,1} = training_features(:,tids==uids(i));
		class_num(i) = sum(tids==uids(i));
	end
	clearvars training_features;
	save([experiment.output.resultsdir experiment.id 'temp.mat'], 'class_features', '-v7.3');
	fprintf('Finished splitting features by class.\n');
	
	% experimental loop
	fprintf('Starting experimental loop.\n');
	for i=1:5
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
		% fix gender ids
		gtrain_ids = train_ids;
		gtrain_ids(train_ids==1 | train_ids==3 | train_ids==5 | train_ids==7 | train_ids==9) = 1;
		gtrain_ids(train_ids==2 | train_ids==4 | train_ids==6 | train_ids==8 | train_ids==10) = 2;
		gtest_ids = test_ids;
		gtest_ids(test_ids==1 | test_ids==3 | test_ids==5 | test_ids==7 | test_ids==9) = 1;
		gtest_ids(test_ids==2 | test_ids==4 | test_ids==6 | test_ids==8 | test_ids==10) = 2;
		
		clearvars class_features;
		fprintf('-Finished splitting train and test set.\n');
        
        if strcmp(experiment.variables.distance,'libsvm')
			fprintf('1. Train SVM Classifier with %d images\n',size(gtrain_ids,1));
			t2 = tic;
			model = svmtrain(gtrain_ids, train_set', '-q -s 0 -t 2');
			toc(t2)
			fprintf('2. Predict Classification with SVM\n');
			t2 = tic;
            pred=zeros(size(gtest_ids,1),1);
			for k=1:size(gtest_ids,1)
				[pred(k,1), ~, ~] = svmpredict(gtest_ids(k), test_set(:,k)', model);
			end
			toc(t2)
		elseif strcmp(experiment.variables.distance,'bayesian')
			O1 = NaiveBayes.fit(train_set',train_ids,'Distribution','normal');
			pred = O1.predict(test_set')';			
        else
            distances = compute_distances(train_avg, test_set, experiment.id, experiment.variables, experiment.input, experiment.output);
            [~,pred] = min(distances);
        end
		fprintf('-Finished 1st order classifying.\n');
		
		
		%%% second order classifying %%%
		gpred = pred;
		% split training set
		male_train_set = train_set(:,gtrain_ids == 1);
		male_train_ids = train_ids(gtrain_ids == 1);
		female_train_set = train_set(:,gtrain_ids == 2);
		female_train_ids = train_ids(gtrain_ids == 2);
		% split test set
		male_test_set = test_set(:,gpred == 1);
		male_test_ids = test_ids(gpred == 1);
		female_test_set = test_set(:,gpred == 2);
		female_test_ids = test_ids(gpred == 2);
		clearvars train_set train_ids test_set test_ids;
		
		% classify males
		if strcmp(experiment.variables.distance,'libsvm')
			fprintf('1. Train SVM Classifier with %d images\n',size(male_train_ids,1));
			t2 = tic;
			model = svmtrain(male_train_ids, male_train_set', '-q -s 0 -t 2');
			toc(t2)
			fprintf('2. Predict Classification with SVM\n');
			t2 = tic;
            mpred=zeros(size(male_test_ids,1),1);
			dec_values=zeros(size(male_test_ids,1),model.nr_class*(model.nr_class-1)/2);
			for k=1:size(male_test_ids,1)
				[mpred(k,1), ~, dec_values(k,:)] = svmpredict(male_test_ids(k), male_test_set(:,k)', model);
			end
			toc(t2)
			fprintf('3. build votes array\n');
			t2 = tic;
			mdistances = zeros(size(male_test_ids,1),model.nr_class);
			index = 1;
			for k=1:size(male_test_ids,1)
			    index=1;
			    for t=1:model.nr_class
			        for j=t+1:model.nr_class
			            if dec_values(k,index) > 0
			                mdistances(k,t)=mdistances(k,t)+1;
			            else
			                mdistances(k,j)=mdistances(k,j)+1;
			            end
			            index=index+1;
			        end 
                end
			end
			toc(t2)
		elseif strcmp(experiment.variables.distance,'bayesian')
			O1 = NaiveBayes.fit(train_set',train_ids,'Distribution','normal');
			pred = O1.predict(test_set')';			
        else
            distances = compute_distances(train_avg, test_set, experiment.id, experiment.variables, experiment.input, experiment.output);
            [~,pred] = min(distances);
        end
		fprintf('-Finished 2nd order classifying for males.\n');
		
		% classify females
		if strcmp(experiment.variables.distance,'libsvm')
			fprintf('1. Train SVM Classifier with %d images\n',size(female_train_ids,1));
			t2 = tic;
			model = svmtrain(female_train_ids, female_train_set', '-q -s 0 -t 2');
			toc(t2)
			fprintf('2. Predict Classification with SVM\n');
			t2 = tic;
            fpred=zeros(size(female_test_ids,1),1);
			dec_values=zeros(size(female_test_ids,1),model.nr_class*(model.nr_class-1)/2);
			for k=1:size(female_test_ids,1)
				[fpred(k,1), ~, dec_values(k,:)] = svmpredict(female_test_ids(k), female_test_set(:,k)', model);
			end
			toc(t2)
			fprintf('3. build votes array\n');
			t2 = tic;
			fdistances = zeros(size(female_test_ids,1),model.nr_class);
			index = 1;
			for k=1:size(female_test_ids,1)
			    index=1;
			    for t=1:model.nr_class
			        for j=t+1:model.nr_class
			            if dec_values(k,index) > 0
			                fdistances(k,t)=fdistances(k,t)+1;
			            else
			                fdistances(k,j)=fdistances(k,j)+1;
			            end
			            index=index+1;
			        end 
                end
			end
			toc(t2)
			
		elseif strcmp(experiment.variables.distance,'bayesian')
			O1 = NaiveBayes.fit(train_set',train_ids,'Distribution','normal');
			pred = O1.predict(test_set')';			
        else
            distances = compute_distances(train_avg, test_set, experiment.id, experiment.variables, experiment.input, experiment.output);
            [~,pred] = min(distances);
        end
		fprintf('-Finished 2nd order classifying for females.\n');
		
		refit_test_ids = [male_test_ids; female_test_ids];
		refit_pred = [mpred; fpred];
		
		distances = [mdistances; fdistances]; clearvars mdistances fdistances;
        save([experiment.output.resultsdir experiment.id 'distance' num2str(i) '.mat'], 'distances');
		
        accuracy(i) = sum(refit_test_ids == refit_pred)/length(refit_test_ids);
		prediction(:,i) = refit_pred;
		truth(:,i) = refit_test_ids;
        for j=1:numel(uids)
            test1 = refit_test_ids(refit_test_ids == j);
            pred1 = refit_pred(refit_test_ids == j);
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
		fprintf('Mean Accuracy for Class %d: %7.4f%%\n', i, results.mean_class_accuracy(i)*100);
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

function [tids gids eids uids] = getClassIDs(experiment)
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
	
	gids = tids;
	gids(tids==1 | tids==3 | tids==5 | tids==7 | tids==9) = 1;
	gids(tids==2 | tids==4 | tids==6 | tids==8 | tids==10) = 2;
	
	eids = tids;
	eids(tids==1 | tids==2) = 1;
	eids(tids==3 | tids==4) = 2;
	eids(tids==5 | tids==6) = 3;
	eids(tids==7 | tids==8) = 4;
	eids(tids==9 | tids==10) = 5;
end
