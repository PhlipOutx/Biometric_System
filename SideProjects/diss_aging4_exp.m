function results = diss_aging4_exp(input)
% input(struct):	Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%					list of required fields.
%					
% input(string):	Filepath to xml file that will be parsed into a struct.
%
% results:			See Section 3.0 for list of outputs.

	
	if isstruct(input)
		experiment = input;
	else
		% parse xml experiment file
		experiment = xml_read(input);
	end
    
    experiment.input.datadir = '\\skynet\static\Processed_Images\Pinellas\Pinellas_ALIGNEDFACE\';
    experiment.output.resultsdir = '\\skynet\users\pemille\results\Dissertation\';
    experiment.variables.printing = 1;
    
	fprintf('###Beginning %s\n', experiment.id);
	fprintf('Image size: %d X %d\n', experiment.variables.data.max_y, experiment.variables.data.max_x);
	fprintf('Patch size: %d X %d\n', experiment.variables.patches.size_y, experiment.variables.patches.size_x);
	fprintf('Patch type: %s\n', experiment.variables.patches.type);
	fprintf('Distance:   %s\n', experiment.variables.distance);
	fprintf('Feature:    %s\n', experiment.variables.feature.name);
	
	% ids
	w1id = 'diss41000001';
	w3id = 'diss41000004';
	w5id = 'diss41000007';
	w7id = 'diss41000010';
	w9id = 'diss41000013';
	w11id = 'diss41000016';
	
	% LBP
	LBPrank1 = zeros(1,16);
	LBPeer = zeros(1,16);
	for j = 0:15
		experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG' num2str(j) '.srt'];
		experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year' num2str(j) '.srt'];
		
		features = load([experiment.output.resultsdir w3id 'g' num2str(j) 'LBP.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w3id 'y' num2str(j) 'LBP.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		
		features = load([experiment.output.resultsdir w5id 'g' num2str(j) 'LBP.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w5id 'y' num2str(j) 'LBP.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = distances + (5*compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output));
		
		features = load([experiment.output.resultsdir w7id 'g' num2str(j) 'LBP.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w7id 'y' num2str(j) 'LBP.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = distances + compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		clearvars gallery_features probe_features;
		
		CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		LBPrank1(j+1) = CMC_rec_rates(1)*100;
		[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		clearvars distances;
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		LBPeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		printing(experiment.variables.printing, 'LBP experiments', j+1, 16);
	end
    fprintf('\nLBP Results: Rank-1\n');
	for j = 1:16 
		fprintf('%f,',LBPrank1(1,j));
	end
	fprintf('\n');
	fprintf('\nLBP Results: EER\n');
	for j = 1:16 
		fprintf('%f,',LBPeer(1,j));
	end
	fprintf('\n\n');
	
	
	
	% HOG
	HOGrank1 = zeros(1,16);
	HOGeer = zeros(1,16);
	for j = 0:15
		experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG' num2str(j) '.srt'];
		experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year' num2str(j) '.srt'];
		
		features = load([experiment.output.resultsdir w3id 'g' num2str(j) 'HOG.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w3id 'y' num2str(j) 'HOG.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		
		features = load([experiment.output.resultsdir w5id 'g' num2str(j) 'HOG.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w5id 'y' num2str(j) 'HOG.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = distances + (5*compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output));
		
		features = load([experiment.output.resultsdir w7id 'g' num2str(j) 'HOG.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w7id 'y' num2str(j) 'HOG.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = distances + compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		clearvars gallery_features probe_features;
		
		
		CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		HOGrank1(j+1) = CMC_rec_rates(1)*100;
		[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
		clearvars distances;
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		HOGeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		printing(experiment.variables.printing, 'HOG experiments', j+1, 16);
	end
    fprintf('\nHOG Results: Rank-1\n');
	for j = 1:16 
		fprintf('%f,',HOGrank1(1,j));
	end
	fprintf('\n');
	fprintf('\nHOG Results: EER\n');
	for j = 1:16 
		fprintf('%f,',HOGeer(1,j));
	end
	fprintf('\n\n');
	
	
	% LPQ
	LPQrank1 = zeros(1,16);
	LPQeer = zeros(1,16);
	for j = 0:15
		experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG' num2str(j) '.srt'];
		experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year' num2str(j) '.srt'];
		
		features = load([experiment.output.resultsdir w7id 'g' num2str(j) 'LPQ.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w7id 'y' num2str(j) 'LPQ.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		
		features = load([experiment.output.resultsdir w9id 'g' num2str(j) 'LPQ.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w9id 'y' num2str(j) 'LPQ.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = distances + (5*compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output));

		features = load([experiment.output.resultsdir w11id 'g' num2str(j) 'LPQ.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir w11id 'y' num2str(j) 'LPQ.mat'], 'features'); probe_features = features.features; clearvars features;
		distances = distances + compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		clearvars gallery_features probe_features;
		
		CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		LPQrank1(j+1) = CMC_rec_rates(1)*100;
		[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
		clearvars distances;
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		LPQeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		printing(experiment.variables.printing, 'LPQ experiments', j+1, 16);
	end
    fprintf('\nLPQ Results: Rank-1\n');
	for j = 1:16 
		fprintf('%f,',LPQrank1(1,j));
	end
	fprintf('\n');
	fprintf('\nLPQ Results: EER\n');
	for j = 1:16 
		fprintf('%f,',LPQeer(1,j));
	end
	fprintf('\n\n');
	

 
end

function extractTheseFeatures(experiment, feat1, list1, name1)

	% feature extraction
	if(~exist([experiment.output.resultsdir experiment.id name1 '.mat'], 'file'))
		experiment.variables.feature.name = feat1;
		
		% extract gallery features
		t1 = tic;
		features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, list1, experiment.input.datadir, '', '');
		s=whos('features');
		fprintf('Produced feature matrix [%d features x %d images] [%6.2f MB]\n', size(features,1), size(features,2), s.bytes/1000000);
		toc(t1)
		fprintf('\n');
		
		save([experiment.output.resultsdir experiment.id name1 '.mat'], 'features', '-v7.3'); clearvars features;
	end

end
