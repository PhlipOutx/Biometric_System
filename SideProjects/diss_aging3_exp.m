function results = diss_aging3_exp(input)
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
	
	
	%extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG.srt', 'gfLBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year0.srt', 'y0LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year1.srt', 'y1LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year2.srt', 'y2LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year3.srt', 'y3LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year4.srt', 'y4LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year5.srt', 'y5LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year6.srt', 'y6LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year7.srt', 'y7LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year8.srt', 'y8LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year9.srt', 'y9LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year10.srt', 'y10LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year11.srt', 'y11LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year12.srt', 'y12LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year13.srt', 'y13LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year14.srt', 'y14LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year15.srt', 'y15LBP');
	
	%extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG.srt', 'gfHOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year0.srt', 'y0HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year1.srt', 'y1HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year2.srt', 'y2HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year3.srt', 'y3HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year4.srt', 'y4HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year5.srt', 'y5HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year6.srt', 'y6HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year7.srt', 'y7HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year8.srt', 'y8HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year9.srt', 'y9HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year10.srt', 'y10HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year11.srt', 'y11HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year12.srt', 'y12HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year13.srt', 'y13HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year14.srt', 'y14HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year15.srt', 'y15HOG');
	
	%extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG.srt', 'gfLPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year0.srt', 'y0LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year1.srt', 'y1LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year2.srt', 'y2LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year3.srt', 'y3LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year4.srt', 'y4LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year5.srt', 'y5LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year6.srt', 'y6LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year7.srt', 'y7LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year8.srt', 'y8LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year9.srt', 'y9LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year10.srt', 'y10LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year11.srt', 'y11LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year12.srt', 'y12LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year13.srt', 'y13LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year14.srt', 'y14LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year15.srt', 'y15LPQ');
	
	
	
	
	
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG0.srt', 'g0LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG1.srt', 'g1LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG2.srt', 'g2LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG3.srt', 'g3LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG4.srt', 'g4LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG5.srt', 'g5LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG6.srt', 'g6LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG7.srt', 'g7LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG8.srt', 'g8LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG9.srt', 'g9LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG10.srt', 'g10LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG11.srt', 'g11LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG12.srt', 'g12LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG13.srt', 'g13LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG14.srt', 'g14LBP');
	% extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG15.srt', 'g15LBP');
	
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG0.srt', 'g0HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG1.srt', 'g1HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG2.srt', 'g2HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG3.srt', 'g3HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG4.srt', 'g4HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG5.srt', 'g5HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG6.srt', 'g6HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG7.srt', 'g7HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG8.srt', 'g8HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG9.srt', 'g9HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG10.srt', 'g10HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG11.srt', 'g11HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG12.srt', 'g12HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG13.srt', 'g13HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG14.srt', 'g14HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG15.srt', 'g15HOG');
	
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG0.srt', 'g0LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG1.srt', 'g1LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG2.srt', 'g2LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG3.srt', 'g3LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG4.srt', 'g4LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG5.srt', 'g5LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG6.srt', 'g6LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG7.srt', 'g7LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG8.srt', 'g8LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG9.srt', 'g9LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG10.srt', 'g10LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG11.srt', 'g11LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG12.srt', 'g12LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG13.srt', 'g13LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG14.srt', 'g14LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG15.srt', 'g15LPQ');

	
	
	
	% LBP
	% LBPrank1 = zeros(1,16);
	% LBPeer = zeros(1,16);
	% for j = 0:15
		% features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'LBP.mat'], 'features'); gallery_features = features.features; clearvars features;
		% features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'LBP.mat'], 'features'); probe_features = features.features; clearvars features;
		
		% experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		% experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		% experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG' num2str(j) '.srt'];
		% experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year' num2str(j) '.srt'];
		% distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		% clearvars gallery_features probe_features;
		% CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% LBPrank1(j+1) = CMC_rec_rates(1)*100;
		% [true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% clearvars distances;
		% [~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		% LBPeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		% printing(experiment.variables.printing, 'LBP experiments', j+1, 16);
	% end
    % fprintf('\nLBP Results: Rank-1\n');
	% for j = 1:16 
		% fprintf('%f,',LBPrank1(1,j));
	% end
	% fprintf('\n');
	% fprintf('\nLBP Results: EER\n');
	% for j = 1:16 
		% fprintf('%f,',LBPeer(1,j));
	% end
	% fprintf('\n\n');
	
	
	
	% HOG
	HOGrank1 = zeros(1,16);
	HOGeer = zeros(1,16);
	for j = 0:15
		features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'HOG.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'HOG.mat'], 'features'); probe_features = features.features; clearvars features;
		
		experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG' num2str(j) '.srt'];
		experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year' num2str(j) '.srt'];
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
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
		features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'LPQ.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'LPQ.mat'], 'features'); probe_features = features.features; clearvars features;
		
		experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_yearG' num2str(j) '.srt'];
		experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_face_year' num2str(j) '.srt'];
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
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
		features = extract_features_dissexp10(experiment.id, experiment.variables, experiment.input, experiment.output, list1, experiment.input.datadir, '', '');
		s=whos('features');
		fprintf('Produced feature matrix [%d features x %d images] [%6.2f MB]\n', size(features,1), size(features,2), s.bytes/1000000);
		toc(t1)
		fprintf('\n');
		
		save([experiment.output.resultsdir experiment.id name1 '.mat'], 'features', '-v7.3'); clearvars features;
	end

end
