function PAEPselect(input)

	if isstruct(input)
		experiment = input;
	else
		experiment = xml_read(input);
	end

	image_data = paep_preproc(experiment.input.probe, experiment.input, experiment.variables);
	image_data2 = paep_preproc(experiment.input.gallery, experiment.input, experiment.variables);
	    
    % get paep variables
	paepmaxima = load ([experiment.output.resultsdir experiment.id 'paepmaxima.mat'], 'paepmaxima');
	paepmaxima = paepmaxima.paepmaxima;
	
	PAEPselectwork(paepmaxima.w25, 'w25', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w23, 'w23', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w21, 'w21', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w19, 'w19', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w17, 'w17', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w15, 'w15', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w13, 'w13', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w11, 'w11', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w9, 'w9', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w7, 'w7', experiment, image_data, image_data2);
	PAEPselectwork(paepmaxima.w5, 'w5', experiment, image_data, image_data2);
	%PAEPselectwork(paepmaxima.w3, 'w3', experiment, image_data, image_data2);
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

function PAEPselectwork(data, window, experiment, image_data, image_data2)
	t1 = tic; 
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
	
    idx = find(data~=0); 
    pixel_perf = data(idx);
    [~,IX] = sort(pixel_perf,'descend');
    pixel_perf = pixel_perf(IX);
    idx = idx(IX);
	[i,j]=ind2sub([100,100],idx);
    cutoff1 = sum(pixel_perf > 1.0);
	cutoff2 = sum(pixel_perf > 0.5);
    
	% setup feature vectors
    if strcmp(experiment.variables.feature.name,'LBP')
        mapping = getmapping(experiment.variables.feature.samples, experiment.variables.feature.type);
        flength = mapping.num;
    elseif strcmp(experiment.variables.feature.name,'HOG')
        flength = experiment.variables.feature.bins;
    elseif strcmp(experiment.variables.feature.name,'LPQ')
        flength = 256;
    elseif strcmp(experiment.variables.feature.name,'WLD')
        flength = experiment.variables.feature.num_orientation*experiment.variables.feature.num_excitation;
    end
    probe_features = zeros(length(idx)*flength,pnum);
    gallery_features = zeros(length(idx)*flength,gnum);
	
    % extract features
    for k = 1:cutoff1
        experiment.variables.cp = [i(k) j(k)];
        p1 = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.probe, image_data);
    	g1 = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.gallery, image_data2);
    	probe_features((k-1)*flength+1:k*flength,:) = p1;
        gallery_features((k-1)*flength+1:k*flength,:) = g1;
		printing(experiment.variables.printing, 'features from patch', k, cutoff1);
    end
    distances = slmetric_pw(gallery_features, probe_features, experiment.variables.distance);
	CMC_rec_rates1 = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
    %fprintf('Rank-1 of top %d pixels: %f\n',cutoff1,CMC_rec_rates1(1)*100);
    
    % extract features
    for k = cutoff1+1:cutoff2
        experiment.variables.cp = [i(k) j(k)];
        p1 = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.probe, image_data);
    	g1 = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.gallery, image_data2);
    	probe_features((k-1)*flength+1:k*flength,:) = p1;
        gallery_features((k-1)*flength+1:k*flength,:) = g1;
		printing(experiment.variables.printing, 'features from patch', k, cutoff2);
    end
	s=whos('probe_features');
    distances = slmetric_pw(gallery_features, probe_features, experiment.variables.distance);
	CMC_rec_rates2 = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
    %fprintf('Rank-1 of all %d pixels: %f\n',cutoff2,CMC_rec_rates2(1)*100);
	
	% extract features
    for k = cutoff2+1:length(idx)
        experiment.variables.cp = [i(k) j(k)];
        p1 = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.probe, image_data);
    	g1 = paep_extract_features(experiment.variables, experiment.input, experiment.output, experiment.input.gallery, image_data2);
    	probe_features((k-1)*flength+1:k*flength,:) = p1;
        gallery_features((k-1)*flength+1:k*flength,:) = g1;
		printing(experiment.variables.printing, 'features from patch', k, length(idx));
    end
	s=whos('probe_features');
    distances = slmetric_pw(gallery_features, probe_features, experiment.variables.distance);
	CMC_rec_rates3 = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
    %fprintf('Rank-1 of all %d pixels: %f\n',length(idx),CMC_rec_rates2(1)*100);
	time = toc(t1);
	
	fprintf('\n%s \t %d \t %f \t %d \t %f \t %d \t %f \t %10.2f s \t %6.2f MB \n', window, cutoff1, CMC_rec_rates1(1)*100, cutoff2, CMC_rec_rates2(1)*100, length(idx), CMC_rec_rates3(1)*100, time, (s.bytes/1000000)*2);
end

