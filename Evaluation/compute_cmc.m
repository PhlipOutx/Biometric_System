function result = compute_cmc(id, distances, variables, input, output)
% INPUTS:
% distances             - a matrix of gallery x probe distances
% gallery_names         - a cell array of gallery names in the format
%                         <subject_id>_<recording_id>.<file_format>
% probe_names           - a cell array of probe names in the format
%                         <subject_id>_<recording_id>.<file_format>
% nameprefix            - length of the subject_id
%
% OUTPUTS:
% result                - recognition rates for each rank
%
% NOTES / COMMENTS
%
% The function was written with Matlab ver. 7.12.0.635 (R2011a).
%

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
    assert(pnum == pnum1, sprintf('Distance matrix does not match imagelists: Probe count mismatch [%d probes][%d distances]', pnum, pnum1));
    assert(gnum == gnum1, sprintf('Distance matrix does not match imagelists: Gallery count mismatch [%d gallery][%d distances]', gnum, gnum1));
    clearvars pnum1 gnum1;

    %%% Initilize variables %%%
    pnum = size(probe_names,1);
    gnum = size(gallery_names,1);

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
    
    %%% find Rank for each probe %%%
    missed = zeros(pnum,1);
    for i = 1:pnum
		%s1 = strncmp(probe_names(i), gallery_names, variables.nameprefix);
		s1 = pids(i) == gids;
		s2 = strcmp(probe_names(i), gallery_names);
        min_probe = min(distances(xor(s1,s2),i));
		if isempty(min_probe), min_probe=bitmax; end
        miss_match = length(unique(gids(distances(~s1,i) <= min_probe)));
        missed(i) = miss_match+1;
    end
    
    sorted = sort(missed);
    rank_count = zeros(1, sorted(pnum));
    for i = 1:pnum
        rank_count(sorted(i)) = rank_count(sorted(i)) + 1;
    end
    cumulative = cumsum(rank_count);
    cumulative = cumulative/pnum;

    result = cumulative;
end
