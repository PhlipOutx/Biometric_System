function results = demographic_classifier_scorefusion(input)
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
    
	% split training lists up by class
	[tids uids] = getClassIDs(experiment);
	fprintf('There are %d unique classes.\n', numel(uids));
	
	results = load([experiment.output.resultsdir experiment.LBP59id experiment.output.results_mat], 'results'); LBP59results = results.results; clearvars results;
	
	
	% experimental loop
	fprintf('Starting experimental loop.\n');
	for i=1:5
		% loading distance files
		results = load([experiment.output.resultsdir experiment.LBP59id 'distance' num2str(i) '.mat'], 'distances'); distances = results.distances; clearvars results;
		results = load([experiment.output.resultsdir experiment.LBP256id 'distance' num2str(i) '.mat'], 'distances'); distances = distances + results.distances; clearvars results;
		results = load([experiment.output.resultsdir experiment.HOG12id 'distance' num2str(i) '.mat'], 'distances'); distances = distances + results.distances; clearvars results;
		results = load([experiment.output.resultsdir experiment.LPQ256id 'distance' num2str(i) '.mat'], 'distances'); distances = distances + results.distances; clearvars results;
		results = load([experiment.output.resultsdir experiment.LCH16id 'distance' num2str(i) '.mat'], 'distances'); distances = distances + results.distances; clearvars results;
		test_ids =  LBP59results.truth(:,i);
		
		[~,pred] = max(distances,[],2);
		
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
