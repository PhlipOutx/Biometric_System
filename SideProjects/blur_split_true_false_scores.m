function [true_scores false_scores top] = split_true_false_scores(id, distances, variables, input, output, iter)

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
	
	fid = fopen(input.blurvalues);
	blur_values = textscan(fid, '%s %f');
	fclose(fid);
	blur_names = blur_values{1};
	blur_values = blur_values{2};
	bnum = size(blur_names,1);

    %%% check that imagelists match distance matrix %%%
    [gnum1 pnum1] = size(distances);
    assert(pnum == pnum1, 'Distance matrix does not match imagelists: Probe count mismatch');
    assert(gnum == gnum1, 'Distance matrix does not match imagelists: Gallery count mismatch');
    clearvars pnum1 gnum1;
	
	%%% find id for each image %%%
	if isempty(input.probe_id) && isempty(input.gallery_id)
		psub = cellfun(@(x) x(1:variables.nameprefix), probe_names, 'UniformOutput', false);
		gsub = cellfun(@(x) x(1:variables.nameprefix), gallery_names, 'UniformOutput', false);

        asub=[psub;gsub];
        [avalues, aloc, aids] = unique(asub);
        
        [trash, psubloc] = ismember(psub,asub);
        pids = aids(psubloc); 
        
        [trash, gsubloc] = ismember(gsub,asub);
        gids = aids(gsubloc);
        
        clearvars avalues aloc trash psubloc gsubloc aids gsub psub;
	else
		fid = fopen(input.probe_id);
		pids = textscan(fid, '%d');
		fclose(fid);
		pids = pids{1};
		assert(pnum == size(pids,1), 'IDs does not match imagelists: Probe count mismatch');
		
		fid = fopen(input.gallery_id);
		gids = textscan(fid, '%d');
		fclose(fid);
		gids = gids{1};
		assert(gnum == size(gids,1), 'IDs does not match imagelists: Gallery count mismatch');
	end
	
	% compute absolute blur difference
	[btf, bindex] = ismember(probe_names, blur_names);
	t1 = blur_values(bindex(1:pnum));
	[btf, bindex] = ismember(gallery_names, blur_names);
	t2 = blur_values(bindex(1:gnum));
	abd = abs(bsxfun(@minus,t1,t2'));
	clearvars t1 t2 btf bindex blur_names blur_values;
	
	min_abd = min(min(abd))-0.001;
	max_abd = max(max(abd))+0.001;
	step = abs(max_abd-min_abd)/10.0;
	range = min_abd+step:step:max_abd;
	top = range(iter);
	bottom = range(iter)-step;

	%%% create list of tm and fm scores %%%
	nts = 0;cts = 1;nfs = 0;cfs = 1;
    for i = 1:pnum
		%s1 = strncmp(probe_names(i), gallery_names, variables.nameprefix);
		s1 = pids(i) == gids;
		s2 = strcmp(probe_names(i), gallery_names);
		nts = nts + sum(((abd(:,i) <= top) & xor(s1,s2)) & ((abd(:,i) > bottom) & xor(s1,s2)));
		%nts = nts + sum(xor(s1,s2));
		nfs = nfs + sum(((abd(:,i) <= top) & ~s1) & ((abd(:,i) > bottom) & ~s1));
		%nfs = nfs + sum(~s1);
	end
    true_scores = zeros(nts,1);
    false_scores = zeros(nfs,1);
	for i = 1:pnum
		% fill true_scores and false_scores
		%s1 = strncmp(probe_names(i), gallery_names, variables.nameprefix);
		s1 = pids(i) == gids;
		s2 = strcmp(probe_names(i), gallery_names);
		ts_l = ((abd(:,i) <= top) & xor(s1,s2)) & ((abd(:,i) > bottom) & xor(s1,s2));
		fs_l = ((abd(:,i) <= top) & ~s1) & ((abd(:,i) > bottom) & ~s1);
		true_scores(cts:cts+sum(ts_l)-1) = distances(ts_l,i);
		false_scores(cfs:cfs+sum(fs_l)-1) = distances(fs_l,i);
		cts = cts + sum(ts_l);
		cfs = cfs + sum(fs_l);
	end
	clearvars s1 s2 temp1 i;

end


