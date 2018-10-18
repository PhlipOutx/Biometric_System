function results = feature_analysis(input)
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
	else
		probe_features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, experiment.input.probe, experiment.input.datadir, experiment.output.probe_mat, experiment.output.resultsdir);
		features = [gallery_features,probe_features];
	end
	s=whos('features');
	fprintf('Produced probe feature matrix [%d features x %d images] [%6.2f MB]\n', size(features,1), size(features,2), s.bytes/1000000);
    toc(t1)
	fprintf('\n');

	
	%%% Scan in image filenames %%%
	[probe_names gallery_names pnum gnum pids gids] = getClassIDs(experiment);


	
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
	end
	
	if strcmp(experiment.variables.patches.type,'miller');
		loc = [1,11;12,22;23,26;27,30;31,48;49,63;64,84];
	elseif strcmp(experiment.variables.patches.type,'miller2');
		loc = [1,7;8,12;13,16;17,20;21,28;29,34;35,74];
	end

	
	%%% Correlation based approach %%%
	r1=zeros(size(gallery_features,1),1);
	for i=1:size(gallery_features,1)
		[r1(i),~,~,~,~,~,~] = ICC([gallery_features(i,:);probe_features(i,:)]', 'A-1');
    end
    r1(isnan(r1)) = 0;
	r1r = reshape(r1,size(gallery_features,1)/flength,flength)';
	
	temp1 = corrcoef(features');
	temp1(isnan(temp1)) = 0;
    r2 = ((sum(temp1,2)-1)/size(features,1));
    r3 = abs(r1)./sqrt(abs(r2));
	r3r = reshape(r3,size(gallery_features,1)/flength,flength)';
	
	r1a = r3 > mean(r3);
    for i=1:size(loc,1)
        r3rm(i,:) = mean(r3r(:,loc(i,1):loc(i,2)),2)';
    end
	r3rm_m = mean(r3rm,2);
    r3rm_50 = r3rm > repmat(r3rm_m,1,flength);
	r3rm_50_3 = sum(r3rm_50) > 3;
	red1 = repmat(r3rm_50(1,:),1,7);
    red1 = [red1 repmat(r3rm_50(2,:),1,5)];
    red1 = [red1 repmat(r3rm_50(3,:),1,4)];
    red1 = [red1 repmat(r3rm_50(4,:),1,4)];
    red1 = [red1 repmat(r3rm_50(5,:),1,8)];
    red1 = [red1 repmat(r3rm_50(6,:),1,6)];
    red1 = [red1 repmat(r3rm_50(7,:),1,40)];
	
	fprintf('Correlation based: \n');
	for i=1:size(r3rm_50,1)
        fprintf('%d: ',i);
        for j=1:size(r3rm_50,2)
            if r3rm_50(i,j) == 1, fprintf('1 '); else fprintf('0 '); end
        end 
        fprintf('\n');
    end
    fprintf('\n');
	fprintf('feature length reduced from %d to %d (%f%%)\n\n',size(gallery_features,1),sum(red1),(sum(red1)/size(gallery_features,1))*100);
	featslen(1) = size(gallery_features,1);
	featslen(2) = sum(red1);
	
	%%% performance based approach %%%
	for i=1:7;
		gf1 = gallery_features(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
		pf1 = probe_features(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
        d1 = compute_distances(gf1, pf1, experiment.id, experiment.variables, experiment.input, experiment.output);
        CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
        score_t(i,1) = CMC_rec_rates(1)*100;
            
		feats = 0:flength:size(gf1,1)-2;
		for j=1:flength
			gf2 = gf1(feats+j,:);
			pf2 = pf1(feats+j,:);
			d1 = compute_distances(gf2, pf2, experiment.id, experiment.variables, experiment.input, experiment.output);
			CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
			score(i,j) = CMC_rec_rates(1)*100;
		end
    end
    
    score_m = mean(score,2);
    score_50 = score > repmat(score_m,1,flength);
    score_50_3 = sum(score_50) > 3;
	red2 = repmat(score_50(1,:),1,7);
    red2 = [red2 repmat(score_50(2,:),1,5)];
    red2 = [red2 repmat(score_50(3,:),1,4)];
    red2 = [red2 repmat(score_50(4,:),1,4)];
    red2 = [red2 repmat(score_50(5,:),1,8)];
    red2 = [red2 repmat(score_50(6,:),1,6)];
    red2 = [red2 repmat(score_50(7,:),1,40)];

	fprintf('Performance based: \n');
    for i=1:size(score_50,1)
        fprintf('%d: ',i);
        for j=1:size(score_50,2)
            if score_50(i,j) == 1, fprintf('1 '); else fprintf('0 '); end
        end 
        fprintf('\n');
    end
    fprintf('\n');
	fprintf('feature length reduced from %d to %d (%f%%)\n\n',size(gallery_features,1),sum(red2),(sum(red2)/size(gallery_features,1))*100);
	featslen(3) = sum(red2);
	
	%combine red1 and red2 creating logical matricies for using 1s and not using 1s
	both = (red1==1)&(red2==1);
    fprintf('feature length using bins where both chosen reduced from %d to %d (%f%%)\n\n',size(gallery_features,1),sum(both),(sum(both)/size(gallery_features,1))*100);
    either = (red1==1)|(red2==1);
    fprintf('feature length using bins where either chosen reduced from %d to %d (%f%%)\n\n',size(gallery_features,1),sum(either),(sum(either)/size(gallery_features,1))*100);
    
	%compute rank1 for both methods
    d1 = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
	CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
	fprintf('Rank-1 using bins where all chosen [%7.4f%%]\n',CMC_rec_rates(1)*100);
	perf(1) = CMC_rec_rates(1)*100;
    
    gfc = gallery_features(red1,:);
	pfc = probe_features(red1,:);
    d1 = compute_distances(gfc, pfc, experiment.id, experiment.variables, experiment.input, experiment.output);
	CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
	fprintf('Rank-1 using bins where correlation chosen [%7.4f%%]\n',CMC_rec_rates(1)*100);
    clearvars pfc gfc;
	perf(2) = CMC_rec_rates(1)*100;
    
    gfp = gallery_features(red2,:);
	pfp = probe_features(red2,:);
    d1 = compute_distances(gfp, pfp, experiment.id, experiment.variables, experiment.input, experiment.output);
	CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
	fprintf('Rank-1 using bins where performance chosen [%7.4f%%]\n',CMC_rec_rates(1)*100);
    clearvars pfp gfp;
	perf(3) = CMC_rec_rates(1)*100;
    
    gfb = gallery_features(both,:);
	pfb = probe_features(both,:);
	featslen(4) = size(gfb,1);
    d1 = compute_distances(gfb, pfb, experiment.id, experiment.variables, experiment.input, experiment.output);
	CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
	fprintf('Rank-1 using bins where both chosen [%7.4f%%]\n',CMC_rec_rates(1)*100);
    clearvars gfb pfb;
	perf(4) = CMC_rec_rates(1)*100;
    
    gfe = gallery_features(either,:);
	pfe = probe_features(either,:);
	featslen(5) = size(gfe,1);
    d1 = compute_distances(gfe, pfe, experiment.id, experiment.variables, experiment.input, experiment.output);
	CMC_rec_rates = compute_cmc(experiment.id, d1, experiment.variables, experiment.input, experiment.output);
	fprintf('Rank-1 using bins where either chosen [%7.4f%%]\n',CMC_rec_rates(1)*100);
    clearvars pfe gfe;
	perf(5) = CMC_rec_rates(1)*100;
	
	fprintf('Recap:\n');
	fprintf('%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n', perf(1),perf(2),perf(3),perf(4),perf(5));
	fprintf('%d\t%d\t%d\t%d\t%d\n', featslen(1),featslen(2),featslen(3),featslen(4),featslen(5));
	
	save([experiment.output.resultsdir experiment.id 'fa.mat'], 'either');
    
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

function [probe_names gallery_names pnum gnum pids gids] = getClassIDs(experiment)
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
	if isempty(experiment.input.probe_id) && isempty(experiment.input.training_id)
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

