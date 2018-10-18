function results = PAEP(input)

	if isstruct(input)
		experiment = input;
	else
		experiment = xml_read(input);
	end

	if exist([experiment.output.resultsdir experiment.id experiment.output.paep_mat], 'file') == 0
		%load ([directories.outputdir filenames.dynamic], 'matrix');
		matrix = ones(experiment.variables.data.max_y, experiment.variables.data.max_x)*-1;
	else
		load ([experiment.output.resultsdir experiment.id experiment.output.paep_mat], 'matrix');
    end	

	image_data = paep_preproc(experiment.input.probe, experiment.input, experiment.variables);
	image_data2 = paep_preproc(experiment.input.gallery, experiment.input, experiment.variables);
	
	%%% Scan in image filenames %%%
	% fid = fopen(experiment.input.probe);
	% probe_names = textscan(fid, '%s');
	% fclose(fid);
	% probe_names = probe_names{1};
	% pnum = size(probe_names,1);
	
	% fid = fopen(experiment.input.gallery);
	% gallery_names = textscan(fid, '%s');
	% fclose(fid);
	% gallery_names = gallery_names{1};
	% gnum = size(gallery_names,1);
	
	%%% find id for each image %%%
	% if isempty(experiment.input.training_id)
		% psub = cellfun(@(x) x(1:experiment.variables.nameprefix), probe_names, 'UniformOutput', false);
		% gsub = cellfun(@(x) x(1:experiment.variables.nameprefix), gallery_names, 'UniformOutput', false);

		% asub=[psub;gsub];
		% [avalues, aloc, aids] = unique(asub);
		
		% [trash, psubloc] = ismember(psub,asub);
		% pids = aids(psubloc); 
		
		% [trash, gsubloc] = ismember(gsub,asub);
		% gids = aids(gsubloc);
		
		% clearvars avalues aloc trash psubloc gsubloc aids gsub psub;
	% else
		% fid = fopen(experiment.input.training_id);
		% pids = textscan(fid, '%d');
		% fclose(fid);
		% pids = pids{1};
		% gids = pids;
		% assert(pnum == size(pids,1), 'IDs does not match imagelists: Probe count mismatch');
	% end
	
	% temp = zeros(gnum,pnum);
	% for i = 1:pnum
		% s1 = pids(i) == gids;
		% s2 = strcmp(probe_names(i), gallery_names);
		% temp(:,i) = xor(s1,s2);
	% end
	% ts_mat=find(temp==1);
	% for i = 1:pnum
		% s1 = pids(i) == gids;
		% s2 = strcmp(probe_names(i), gallery_names);
		% temp(:,i) = ~s1;
	% end
	% fs_mat=find(temp==1);
	% clearvars s1 s2 i pids gids gallery_names probe_names pnum gnum temp; 
	
	time_total=0;
	completed=0;
	for i = 1:experiment.variables.data.max_y
		for j = 1:experiment.variables.data.max_x
			if matrix(i,j) == -1
				experiment.variables.cp = [i j];
				t1 = tic;
				probe_features = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.probe, image_data);
				gallery_features = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.gallery, image_data2);
				distances = slmetric_pw(gallery_features, probe_features, experiment.variables.distance);
				%true_scores = distances(ts_mat);
				%false_scores = distances(fs_mat);
				%clearvars distances;
				%score = paep_roc(true_scores, false_scores, experiment.output.roc_resolution)*100;
				%[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
				%score = ROC_rates_and_threshs.EER_er*100;
				CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
				score = CMC_rec_rates(1)*100;
				matrix(i,j) = score;
				time = toc(t1);
				time_total = time_total + time;
				completed = completed + 1;
				
				stemp = sprintf('[%6.2f/%6.2f seconds   %5.2f score   %d/%d j] running iteration', time, time_total/completed, score, j, experiment.variables.data.max_x);
				printing(experiment.variables.printing, stemp, i, experiment.variables.data.max_y);
				save ([experiment.output.resultsdir experiment.id experiment.output.paep_mat], 'matrix');
			end
		end
	end
	
end

function image_data = paep_preproc(imagelist,input,variables)
	%%% Scan in image filenames %%%
	fid = fopen(imagelist);
	images = textscan(fid, '%s');
	fclose(fid);
	images = images{1};
	inum = size(images,1);
	
	% preprocess image
	image_data = cell(inum,1);
	for i = 1:inum
		temp = char(images(i));
		im = imread([input.datadir temp]);
		image_data{i} = ipreproc(im, variables);
		printing(variables.printing, 'preprocessing image', i, inum);
	end
	fprintf('\n');
end


