function results = diss_quality2_exp(input)
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
    
	fprintf('###Beginning %s\n', experiment.id);
	
	
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_lefteye.srt', 'gfLBP');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_lefteye.srt', 'gfHOG');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_lefteye.srt', 'gfLPQ');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_lefteye.srt', 'gfSIFT');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_lefteye.srt', 'gfGABOR');
	
	
	% LBP
	features = load([experiment.output.resultsdir experiment.id 'gfLBP.mat'], 'features'); gfLBP = features.features; clearvars features;
	distances = compute_distances(gfLBP, gfLBP, experiment.id, experiment.variables, experiment.input, experiment.output); clearvars gfLBP;
	%split true false scores
	[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); clearvars distances;
	
	% find quality difference
	fid = fopen('\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_qvalue.srt');
	qvalues = textscan(fid, '%f');
	fclose(fid);
	qvalues = qvalues{1};
	qdist = triu(abs(bsxfun(@minus,repmat(qvalues,1,length(qvalues)),qvalues')),1);
	
	% split true flase scores
	[qtrue_scores qfalse_scores] = split_true_false_scores(experiment.id, qdist, experiment.variables, experiment.input, experiment.output); clearvars qdist;
	true_scores = true_scores(qtrue_scores~=0);
	false_scores = false_scores(qfalse_scores~=0);
	qtrue_scores = qtrue_scores(qtrue_scores~=0);
	qfalse_scores = qfalse_scores(qfalse_scores~=0);
	
	% generate results
	for j = 1:10
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores(qtrue_scores>((j-1)*5) & qtrue_scores<=(j*5)), false_scores(qfalse_scores>((j-1)*5) & qfalse_scores<=(j*5)), experiment.output.roc_resolution);
		LBPeer(j) = ROC_rates_and_threshs.EER_er*100;
	end
	fprintf('Finished LBP experiment\n');
	
	
	% HOG
	features = load([experiment.output.resultsdir experiment.id 'gfHOG.mat'], 'features'); gfHOG = features.features; clearvars features;
	distances = compute_distances(gfHOG, gfHOG, experiment.id, experiment.variables, experiment.input, experiment.output); clearvars gfHOG;
	%split true false scores
	[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); clearvars distances;
	
	% find quality difference
	fid = fopen('\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_qvalue.srt');
	qvalues = textscan(fid, '%f');
	fclose(fid);
	qvalues = qvalues{1};
	qdist = triu(abs(bsxfun(@minus,repmat(qvalues,1,length(qvalues)),qvalues')),1);
	
	% split true flase scores
	[qtrue_scores qfalse_scores] = split_true_false_scores(experiment.id, qdist, experiment.variables, experiment.input, experiment.output); clearvars qdist;
	true_scores = true_scores(qtrue_scores~=0);
	false_scores = false_scores(qfalse_scores~=0);
	qtrue_scores = qtrue_scores(qtrue_scores~=0);
	qfalse_scores = qfalse_scores(qfalse_scores~=0);
	
	% generate results
	for j = 1:10
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores(qtrue_scores>((j-1)*5) & qtrue_scores<=(j*5)), false_scores(qfalse_scores>((j-1)*5) & qfalse_scores<=(j*5)), experiment.output.roc_resolution);
		HOGeer(j) = ROC_rates_and_threshs.EER_er*100;
	end
	fprintf('Finished HOG experiment\n');

	
	
	
	% LPQ
	features = load([experiment.output.resultsdir experiment.id 'gfLPQ.mat'], 'features'); gfLPQ = features.features; clearvars features;
	distances = compute_distances(gfLPQ, gfLPQ, experiment.id, experiment.variables, experiment.input, experiment.output); clearvars gfLPQ;
	%split true false scores
	[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); clearvars distances;
	
	% find quality difference
	fid = fopen('\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_qvalue.srt');
	qvalues = textscan(fid, '%f');
	fclose(fid);
	qvalues = qvalues{1};
	qdist = triu(abs(bsxfun(@minus,repmat(qvalues,1,length(qvalues)),qvalues')),1);
	
	% split true flase scores
	[qtrue_scores qfalse_scores] = split_true_false_scores(experiment.id, qdist, experiment.variables, experiment.input, experiment.output); clearvars qdist;
	true_scores = true_scores(qtrue_scores~=0);
	false_scores = false_scores(qfalse_scores~=0);
	qtrue_scores = qtrue_scores(qtrue_scores~=0);
	qfalse_scores = qfalse_scores(qfalse_scores~=0);
	
	% generate results
	for j = 1:10
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores(qtrue_scores>((j-1)*5) & qtrue_scores<=(j*5)), false_scores(qfalse_scores>((j-1)*5) & qfalse_scores<=(j*5)), experiment.output.roc_resolution);
		LPQeer(j) = ROC_rates_and_threshs.EER_er*100;
	end
	fprintf('Finished LPQ experiment\n');
	
		
	% SIFT
	experiment.variables.feature.name = 'SIFT';
	features = load([experiment.output.resultsdir experiment.id 'gfSIFT.mat'], 'features'); gfSIFT = features.features; clearvars features;
	distances = compute_distances(gfSIFT, gfSIFT, experiment.id, experiment.variables, experiment.input, experiment.output); clearvars gfSIFT;
	%split true false scores
	[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); clearvars distances;
	
	% find quality difference
	fid = fopen('\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_qvalue.srt');
	qvalues = textscan(fid, '%f');
	fclose(fid);
	qvalues = qvalues{1};
	qdist = triu(abs(bsxfun(@minus,repmat(qvalues,1,length(qvalues)),qvalues')),1);
	
	% split true flase scores
	[qtrue_scores qfalse_scores] = split_true_false_scores(experiment.id, qdist, experiment.variables, experiment.input, experiment.output); clearvars qdist;
	true_scores = true_scores(qtrue_scores~=0);
	false_scores = false_scores(qfalse_scores~=0);
	qtrue_scores = qtrue_scores(qtrue_scores~=0);
	qfalse_scores = qfalse_scores(qfalse_scores~=0);
	
	% generate results
	for j = 1:10
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores(qtrue_scores>((j-1)*5) & qtrue_scores<=(j*5)), false_scores(qfalse_scores>((j-1)*5) & qfalse_scores<=(j*5)), experiment.output.roc_resolution);
		SIFTeer(j) = ROC_rates_and_threshs.EER_er*100;
	end
	fprintf('Finished SIFT experiment\n');
	
	
	% GABOR
	experiment.variables.feature.name = 'GABOR';
	features = load([experiment.output.resultsdir experiment.id 'gfGABOR.mat'], 'features'); gfGABOR = features.features; clearvars features;
	distances = compute_distances(gfGABOR, gfGABOR, experiment.id, experiment.variables, experiment.input, experiment.output); clearvars gfGABOR;
	%split true false scores
	[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); clearvars distances;
	
	% find quality difference
	fid = fopen('\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_quality_qvalue.srt');
	qvalues = textscan(fid, '%f');
	fclose(fid);
	qvalues = qvalues{1};
	qdist = triu(abs(bsxfun(@minus,repmat(qvalues,1,length(qvalues)),qvalues')),1);
	
	% split true flase scores
	[qtrue_scores qfalse_scores] = split_true_false_scores(experiment.id, qdist, experiment.variables, experiment.input, experiment.output); clearvars qdist;
	true_scores = true_scores(qtrue_scores~=0);
	false_scores = false_scores(qfalse_scores~=0);
	qtrue_scores = qtrue_scores(qtrue_scores~=0);
	qfalse_scores = qfalse_scores(qfalse_scores~=0);
	
	% generate results
	for j = 1:10
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores(qtrue_scores>((j-1)*5) & qtrue_scores<=(j*5)), false_scores(qfalse_scores>((j-1)*5) & qfalse_scores<=(j*5)), experiment.output.roc_resolution);
		GABOReer(j) = ROC_rates_and_threshs.EER_er*100;
	end
	fprintf('Finished GABOR experiment\n');
	
	
	
	% Print Results
	for i=1:10
		fprintf('%f,',LBPeer(i));
	end
	fprintf('\n');
	for i=1:10
		fprintf('%f,',HOGeer(i));
	end
	fprintf('\n');
	for i=1:10
		fprintf('%f,',LPQeer(i));
	end
	fprintf('\n');
	for i=1:10
		fprintf('%f,',SIFTeer(i));
	end
	fprintf('\n');
	for i=1:10
		fprintf('%f,',GABOReer(i));
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

function [probe_names,gallery_names,pids,gids] = getIDsAndNames(experiment)

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
	
end
