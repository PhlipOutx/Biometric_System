function features = extract_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir)

	if strcmp(variables.feature.name,'LBP') || strcmp(variables.feature.name,'HOG') || strcmp(variables.feature.name,'LPQ') || strcmp(variables.feature.name,'WLD')  || strcmp(variables.feature.name,'LCH')
		features = extract_local_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir);
	elseif strcmp(variables.feature.name,'PCA') || strcmp(variables.feature.name,'LDA') || strcmp(variables.feature.name,'KFA') || strcmp(variables.feature.name,'KPCA') || strcmp(variables.feature.name,'ICA')
		features = extract_subspace_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir);
	elseif strcmp(variables.feature.name,'SIFT') || strcmp(variables.feature.name,'SURF')
		features = extract_keypoint_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir);
	elseif strcmp(variables.feature.name,'GABOR')
		features = extract_gabor_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir);
	end
	
end


function features = extract_subspace_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir)

	% check for already created subspace training matrix
    if ~isempty(output.training_mat) && exist([outputdir id output.training_mat],'file')
        model = load([outputdir id output.training_mat], 'model');
        model = model.model;
    else
		% Scan in image filenames
        fid = fopen(input.training);
        images = textscan(fid, '%s');
        fclose(fid);
        images = images{1};
        inum = size(images,1);
		
		% allocate space
		temp = char(images(1));
		I = imread([datadir temp]);
		I = ipreproc(I, variables);
        train_data = zeros(numel(I), inum);
		
		% fill matrix
		fprintf('1. Loading training data [%d images]\n', inum);
		t2 = tic;
		for i = 1:inum
			% preprocess image
            temp = char(images(i));
            I = imread([datadir temp]);
            I = ipreproc(I, variables);
			% fill matrix
			train_data(:,i) = I(:);
			%printing(variables.printing, 'load training data', i, inum);
		end
		s=whos('train_data');
		fprintf('Training matrix is %d features x %d images [%6.2f MB]\n', size(train_data,1), size(train_data,2), s.bytes/1000000);
		clearvars I temp i fid inum s;
		toc(t2)
		
		%%% subspace training %%%
		fprintf('2. Computing %s subspace model\n', variables.feature.name);
		t2 = tic;
		if strcmp(variables.feature.name,'PCA')
			% build pca subspace
			model = pca_subspace(train_data, variables.feature.vectors);
			clearvars train_data;
		elseif strcmp(variables.feature.name,'LDA')
			% find subject id for each image
			if isempty(input.training_id)
				isub = cellfun(@(x) x(1:variables.nameprefix), images, 'UniformOutput', false);
				[ivalues, iloc, iids] = unique(isub);
				clearvars ivalues iloc images isub;
			else
				fid = fopen(input.training_id);
				iids = textscan(fid, '%d');
				fclose(fid);
				iids = iids{1};
			end
			% build lda subspace
			model = lda_subspace(train_data, iids, variables.feature.vectors);
			clearvars train_data iids;
		elseif strcmp(variables.feature.name,'KFA')
			% find subject id for each image
			if isempty(input.training_id)
				isub = cellfun(@(x) x(1:variables.nameprefix), images, 'UniformOutput', false);
				[ivalues, iloc, iids] = unique(isub);
				clearvars ivalues iloc images isub;
			else
				fid = fopen(input.training_id);
				iids = textscan(fid, '%d');
				fclose(fid);
				iids = iids{1};
			end
			% build kfa subspace
			model = kfa_subspace(train_data, iids, 'fpp', [0 0.7], variables.feature.vectors);
			clearvars train_data iids;
		elseif strcmp(variables.feature.name,'KPCA')
			model = kpca_subspace(train_data, 'poly', [0 2], variables.feature.vectors);
			clearvars train_data;
		elseif strcmp(variables.feature.name,'ICA')
			% build ica subspace
			[model, A, W] = fastica(train_data', 'verbose', 'off', 'lastEig', variables.feature.vectors, 'numOfIC', variables.feature.vectors);
			clearvars train_data;
		end
		s=whos('model');
		fprintf('Model is %d features x %d subspace faces [%6.2f MB] and captures %4.2f%% of variance\n', size(model.W,1), size(model.W,2), s.bytes/1000000, model.varcap*100);
		
		% save model
        if ~isempty(output.training_mat)
            if(~exist(outputdir, 'dir'))
               mkdir(outputdir);
            end
            save([outputdir id output.training_mat], 'model');
        end
		toc(t2)
	end

	% check for already created subspace features
    if ~isempty(savefeatures) && exist([outputdir id savefeatures],'file')
        features = load([outputdir id savefeatures], 'features');
        features = features.features;
    else
		% Scan in image filenames
        fid = fopen(imagelist);
        images = textscan(fid, '%s');
        fclose(fid);
        images = images{1};
        inum = size(images,1);
		
		% allocate space
		temp = char(images(1));
		I = imread([datadir temp]);
		I = ipreproc(I, variables);
        test_data = zeros(numel(I), inum);
		
		% fill matrix
		fprintf('1. Loading test data [%d images]\n', inum);
		t2 = tic;
		for i = 1:inum
			% preprocess image
            temp = char(images(i));
            I = imread([datadir temp]);
            I = ipreproc(I, variables);
			% fill matrix
			test_data(:,i) = I(:);
			%printing(variables.printing, 'load test data', i, inum);
		end
		s=whos('test_data');
		fprintf('Test matrix is %d features x %d images [%6.2f MB]\n', size(test_data,1), size(test_data,2), s.bytes/1000000);
		clearvars I temp i fid images inum s;
		toc(t2)
		
		%%% subspace projection %%%
		if strcmp(variables.feature.name,'PCA') || strcmp(variables.feature.name,'LDA') || strcmp(variables.feature.name,'ICA')
			features = linear_subspace_projection(test_data, model, 1);
			clearvars test_data;
		elseif strcmp(variables.feature.name,'KFA') || strcmp(variables.feature.name,'KPCA')
			features = nonlinear_subspace_projection(test_data, model);
			clearvars test_data;
		end
		
		% save features
        if ~isempty(savefeatures)
            if(~exist(outputdir, 'dir'))
               mkdir(outputdir);
            end
            save([outputdir id savefeatures], 'features', '-v7.3');
        end
		
	end

end



function features = extract_local_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir)

    % check for already created matrix
    if ~isempty(savefeatures) && exist([outputdir id savefeatures],'file')
        features = load([outputdir id savefeatures], 'features');
        features = features.features;
    else        
        %%% Scan in image filenames %%%
        fid = fopen(imagelist);
        images = textscan(fid, '%s');
        fclose(fid);
        images = images{1};
        inum = size(images,1);

        %%% setup variables %%%
		if strcmp(variables.feature.name,'LBP') &&  strcmp(variables.feature.type,'all')
			flength = 256;
        elseif strcmp(variables.feature.name,'LBP')
            mapping = getmapping(variables.feature.samples, variables.feature.type);
            flength = mapping.num;
		elseif strcmp(variables.feature.name,'HOG')
			flength = variables.feature.bins;
		elseif strcmp(variables.feature.name,'LPQ')
			flength = 256;
		elseif strcmp(variables.feature.name,'WLD')
			flength = variables.feature.num_orientation*variables.feature.num_excitation;
		elseif strcmp(variables.feature.name,'LCH')
			flength = 16;
        end
        patches = getpatches(variables);

        % allocate space
        features = zeros(patches.num*flength,inum);
		fprintf('#Extracting %d features from %d patches\n', flength, patches.num);

        %%% compute features %%%
        fprintf('1. Computing %s features on %d images\n', variables.feature.name, inum);
		if strcmp(variables.feature.name,'LBP') &&  strcmp(variables.feature.type,'all')
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% split into blocks
				if strcmp(variables.patches.type,'park')
					variables.data.max_x = size(I,2);
					variables.data.max_y = size(I,1);
					patches = getpatches(variables);
				end
				for j = 1:patches.num
					H1 = lbp(I(patches.y(j):patches.dy(j),patches.x(j):patches.dx(j)))';
					features((j-1)*flength+1:j*flength,i) = (H1/sum(H1));
				end
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		elseif strcmp(variables.feature.name,'LBP')  
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% split into blocks
				if strcmp(variables.patches.type,'park')
					variables.data.max_x = size(I,2);
					variables.data.max_y = size(I,1);
					patches = getpatches(variables);
				end
				for j = 1:patches.num
					H1 = lbp(I(patches.y(j):patches.dy(j),patches.x(j):patches.dx(j)), variables.feature.radius, mapping.samples, mapping, 'h')';
					features((j-1)*flength+1:j*flength,i) = (H1/sum(H1));
				end
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		elseif strcmp(variables.feature.name,'HOG')
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% split into blocks
				if strcmp(variables.patches.type,'park')
					variables.data.max_x = size(I,2);
					variables.data.max_y = size(I,1);
					patches = getpatches(variables);
				end
				for j = 1:patches.num
					H1 = hog_dissexp10(I(patches.y(j):patches.dy(j),patches.x(j):patches.dx(j)), variables.feature.bins, variables.feature.radius);
                    H1 = (H1/sum(H1));
					% NaN fix for GOH
					H1(isnan(H1)) = 0;
					features((j-1)*flength+1:j*flength,i) = H1;
				end
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		elseif strcmp(variables.feature.name,'LPQ')
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% split into blocks
				if strcmp(variables.patches.type,'park')
					variables.data.max_x = size(I,2);
					variables.data.max_y = size(I,1);
					patches = getpatches(variables);
				end
				for j = 1:patches.num
					H1 = lpq(I(patches.y(j):patches.dy(j),patches.x(j):patches.dx(j)), variables.feature.winSize, variables.feature.decorr, variables.feature.freqestim, variables.feature.mode);
					features((j-1)*flength+1:j*flength,i) = (H1/sum(H1));
				end
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		elseif strcmp(variables.feature.name,'WLD')
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% split into blocks
				if strcmp(variables.patches.type,'park')
					variables.data.max_x = size(I,2);
					variables.data.max_y = size(I,1);
					patches = getpatches(variables);
				end
				for j = 1:patches.num
					H1 = wld(I(patches.y(j):patches.dy(j),patches.x(j):patches.dx(j)), variables.feature.num_orientation, variables.feature.num_excitation);
					features((j-1)*flength+1:j*flength,i) = (H1/sum(H1));
				end
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		elseif strcmp(variables.feature.name,'LCH')
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% split into blocks
				if strcmp(variables.patches.type,'park')
					variables.data.max_x = size(I,2);
					variables.data.max_y = size(I,1);
					patches = getpatches(variables);
				end
				for j = 1:patches.num
					H1 = lch(I(patches.y(j):patches.dy(j),patches.x(j):patches.dx(j),:), [1 2], [0 1; 0 1]);
					features((j-1)*flength+1:j*flength,i) = (H1/sum(H1));
				end
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		end

        % save matrix
        if ~isempty(savefeatures)
            if(~exist(outputdir, 'dir'))
               mkdir(outputdir);
            end
            save([outputdir id savefeatures], 'features', '-v7.3');
        end
    end

end


function features = extract_keypoint_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir)

    % check for already created matrix
    if ~isempty(savefeatures) && exist([outputdir id savefeatures],'file')
        features = load([outputdir id savefeatures], 'features');
        features = features.features;
    else        
        %%% Scan in image filenames %%%
        fid = fopen(imagelist);
        images = textscan(fid, '%s');
        fclose(fid);
        images = images{1};
        inum = size(images,1);

        %%% setup variables %%%
        if strcmp(variables.feature.name,'SIFT')
            flength = 1;
		elseif strcmp(variables.feature.name,'SURF')
			flength = 1;
			Options.verbose=false;
        end

        % allocate space
        features = cell(1,inum);

        %%% compute features %%%
        fprintf('1. Computing %s features on %d images\n', variables.feature.name, inum);
		if strcmp(variables.feature.name,'SIFT')  
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% SIFT calculations
				IS = single(I);
				[F1,D1] = vl_sift(IS);
				features{i} = D1;
				
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		elseif strcmp(variables.feature.name,'SURF')
			for i = 1:inum
				% preprocess image
				temp = char(images(i));
				I = imread([datadir temp]);
				I = ipreproc(I, variables);

				% SURF calculations
				vec = OpenSurf(I, Options);
				if isempty(fieldnames(vec))
					features{i} = zeros(64,1);
				else
					features{i} = reshape([vec.descriptor],64,[]);
				end
		
				printing(variables.printing, variables.feature.name, i, inum);
			end
			fprintf('\n');
		end

        % save matrix
        if ~isempty(savefeatures)
            if(~exist(outputdir, 'dir'))
               mkdir(outputdir);
            end
            save([outputdir id savefeatures], 'features', '-v7.3');
        end
    end

end

function features = extract_gabor_features_dissexp10(id, variables, input, output, imagelist, datadir, savefeatures, outputdir)

    % check for already created matrix
    if ~isempty(savefeatures) && exist([outputdir id savefeatures],'file')
        features = load([outputdir id savefeatures], 'features');
        features = features.features;
    else        
        %%% Scan in image filenames %%%
        fid = fopen(imagelist);
        images = textscan(fid, '%s');
        fclose(fid);
        images = images{1};
        inum = size(images,1);

        %%% setup variables %%%
		filter_bank = construct_Gabor_filters_PhD(variables.feature.gabor_o, variables.feature.gabor_s, [variables.data.max_y variables.data.max_x]);

        % allocate space
        %features = cell(1,inum);

        %%% compute features %%%
        fprintf('1. Computing %s features on %d images\n', variables.feature.name, inum);

		
		for i = 1:inum
			% preprocess image
			temp = char(images(i));
			I = imread([datadir temp]);
			I = ipreproc(I, variables);

			% Gabor calculations
			I = double(I);
			H1 = filter_image_with_Gabor_bank_PhD(I,filter_bank,64);
			features(:,i) = H1;
			
			printing(variables.printing, variables.feature.name, i, inum);
		end
		fprintf('\n');


        % save matrix
        if ~isempty(savefeatures)
            if(~exist(outputdir, 'dir'))
               mkdir(outputdir);
            end
            save([outputdir id savefeatures], 'features', '-v7.3');
        end
    end

end




