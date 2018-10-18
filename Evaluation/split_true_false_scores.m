function [true_scores false_scores] = split_true_false_scores(id, distances, variables, input, output)

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

	%%% create list of tm and fm scores %%%
	nts = 0;cts = 1;nfs = 0;cfs = 1;
    for i = 1:pnum
		%s1 = strncmp(probe_names(i), gallery_names, variables.nameprefix);
		s1 = pids(i) == gids;
		s2 = strcmp(probe_names(i), gallery_names);
		nts = nts + sum(xor(s1,s2));
		nfs = nfs + sum(~s1);
		printing(variables.printing, 'calculate number of true and false scores', i, pnum);
	end
	fprintf('\n');
    true_scores = zeros(nts,1);
    false_scores = zeros(nfs,1);
	for i = 1:pnum
		% fill true_scores and false_scores
		%s1 = strncmp(probe_names(i), gallery_names, variables.nameprefix);
		s1 = pids(i) == gids;
		s2 = strcmp(probe_names(i), gallery_names);
		true_scores(cts:cts+sum(xor(s1,s2))-1) = distances(xor(s1,s2),i);
		false_scores(cfs:cfs+sum(~s1)-1) = distances(~s1,i);
		cts = cts + sum(xor(s1,s2));
		cfs = cfs + sum(~s1);
		printing(variables.printing, 'generate true and false scores lists', i, pnum);
	end
	fprintf('\n');
	clearvars s1 s2 i;

end


