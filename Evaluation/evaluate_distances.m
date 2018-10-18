function results = evaluate_distances(id, distances, variables, input, output)

    %%% Scan in image filenames %%%
    fid = fopen(input.probe);
    probe_names = textscan(fid, '%s');
    fclose(fid);
    probe_names = probe_names{1};
    pnum = size(probe_names,1);
    
    fid = fopen(input.gallery);
    gallery_names = textscan(fid, '%s');
    fclose(fid);
    gallery_names = gallery_names{1};
    gnum = size(gallery_names,1);

    %%% check that imagelists match distance matrix %%%
    [gnum1 pnum1] = size(distances);
    assert(pnum == pnum1, 'Distance matrix does not match imagelists: Probe count mismatch');
    assert(gnum == gnum1, 'Distance matrix does not match imagelists: Gallery count mismatch');
    clearvars pnum1 gnum1;

	%%% create list of tm and fm scores %%%
	maxd = max(max(distances));
    true_scores = [];
    false_scores = [];
    for i = 1:pnum
		temp1 = probe_names(i);
		% fill true_scores and false_scores
		s1 = strncmp(temp1, gallery_names, variables.nameprefix);
		s2 = strcmp(temp1, gallery_names);
		true_scores = [true_scores; distances(xor(s1,s2),i)];
		false_scores = [false_scores; distances(~s1,i)];
	end
	clearvars s1 s2 temp1 i;

	% disp(size(true_scores));
	% disp(max(true_scores));
	% disp(min(true_scores));
	% disp(size(false_scores));
	% disp(max(false_scores));
	% disp(min(false_scores));
	
	%%% compute CMC %%%
	fprintf('1. Computing CMC\n');
	t2 = tic;
	rec_rates = compute_cmc(distances, gallery_names, probe_names, variables.nameprefix);
	toc(t2)
	% set CMC outputs
	results.CMC_rec_rates = rec_rates;

	%%% compute ROC %%%
	fprintf('2. Computing ROC\n');
	t2 = tic;
	[ver_rate, miss_rate, rates_and_threshs] = compute_roc(true_scores, false_scores, output.roc_resolution);
	toc(t2)
	% set ROC outputs
	results.ROC_ver_rate = ver_rate;
	results.ROC_miss_rate = miss_rate;
	results.ROC_rates_and_threshs = rates_and_threshs;
	
	%%% compute DET %%%
	fprintf('3. Computing DET\n');
	t2 = tic;
	[pmiss, pfa] = compute_det(-true_scores, -false_scores, output.det_resolution);
	toc(t2)
	% set DET outputs
	results.DET_frr_rate = pmiss;
	results.DET_far_rate = pfa;
	clearvars true_scores false_scores;
	
	%%% print points on each curve to a file %%%
	if output.make_files == 1
		fprintf('4. Print points to file\n');
		t2 = tic;
		file = fopen([output.resultsdir id '_det.txt'], 'w');
		for i = 1:length(pfa)
		   fprintf(file, '%f,%f\n', pfa(i), pmiss(i));
		end
		fclose(file);

		file = fopen([output.resultsdir id '_roc.txt'], 'w');
		for i = 1:length(ver_rate)
		   fprintf(file, '%f,%f\n', miss_rate(i), ver_rate(i));
		end
		fclose(file);

		file = fopen([output.resultsdir id '_cmc.txt'], 'w');
		for i = 1:length(rec_rates)
		   fprintf(file, '%f\n', rec_rates(i));
		end
		fclose(file);
		toc(t2)	
	end

end


