function features = paep_extract_features(variables, input, output, imagelist, image_data)
	%%% Scan in image filenames %%%
	fid = fopen(imagelist);
	images = textscan(fid, '%s');
	fclose(fid);
	images = images{1};
	inum = size(images,1);
	
	% setup variables
	if strcmp(variables.feature.name,'LBP')
		mapping = getmapping(variables.feature.samples, variables.feature.type);
		flength = mapping.num;
	elseif strcmp(variables.feature.name,'HOG')
		flength = variables.feature.bins;
	elseif strcmp(variables.feature.name,'LPQ')
		flength = 256;
	elseif strcmp(variables.feature.name,'WLD')
		flength = variables.feature.num_orientation*variables.feature.num_excitation;
	end
	%patches = getpatches(variables.patches.type, variables.patches.size_x, variables.patches.size_y, variables.data.max_x, variables.data.max_y);
	
	% allocate space
	features = zeros(flength,inum);
	
	%%% compute features %%%
	if strcmp(variables.feature.name,'LBP')  
		for i = 1:inum
			% split into blocks
			b1 = variables.cp(1) - floor(variables.patches.size_y/2);
			b2 = variables.cp(1) + floor(variables.patches.size_y/2);
			b3 = variables.cp(2) - floor(variables.patches.size_x/2);
			b4 = variables.cp(2) + floor(variables.patches.size_x/2);
			if b1 < 1, b1 = 1; end
			if b2 > variables.data.max_y, b2 = variables.data.max_y; end
			if b3 < 1, b3 = 1; end
			if b4 > variables.data.max_x, b4 = variables.data.max_x; end
			H1 = lbp(image_data{i}(b1:b2,b3:b4), variables.feature.radius, mapping.samples, mapping, 'h')';
			features(:,i) = H1;
			%printing(variables.printing, variables.feature.name, i, inum);
		end
	elseif strcmp(variables.feature.name,'HOG')
		for i = 1:inum
			% split into blocks
			b1 = variables.cp(1) - floor(variables.patches.size_y/2);
			b2 = variables.cp(1) + floor(variables.patches.size_y/2);
			b3 = variables.cp(2) - floor(variables.patches.size_x/2);
			b4 = variables.cp(2) + floor(variables.patches.size_x/2);
			if b1 < 1, b1 = 1; end
			if b2 > variables.data.max_y, b2 = variables.data.max_y; end
			if b3 < 1, b3 = 1; end
			if b4 > variables.data.max_x, b4 = variables.data.max_x; end
			H1 = hog(image_data{i}(b1:b2,b3:b4), variables.feature.bins);
			% NaN fix for GOH
			H1(isnan(H1)) = 0;
			features(:,i) = H1;
			%printing(variables.printing, variables.feature.name, i, inum);
		end
	elseif strcmp(variables.feature.name,'LPQ')
		for i = 1:inum
			% split into blocks
			b1 = variables.cp(1) - floor(variables.patches.size_y/2);
			b2 = variables.cp(1) + floor(variables.patches.size_y/2);
			b3 = variables.cp(2) - floor(variables.patches.size_x/2);
			b4 = variables.cp(2) + floor(variables.patches.size_x/2);
			if b1 < 1, b1 = 1; end
			if b2 > variables.data.max_y, b2 = variables.data.max_y; end
			if b3 < 1, b3 = 1; end
			if b4 > variables.data.max_x, b4 = variables.data.max_x; end
			H1 = lpq(image_data{i}(b1:b2,b3:b4), variables.feature.winSize, variables.feature.decorr, variables.feature.freqestim, variables.feature.mode);
			features(:,i) = H1;
			%printing(variables.printing, variables.feature.name, i, inum);
		end
	elseif strcmp(variables.feature.name,'WLD')
		for i = 1:inum
			% split into blocks
			b1 = variables.cp(1) - floor(variables.patches.size_y/2);
			b2 = variables.cp(1) + floor(variables.patches.size_y/2);
			b3 = variables.cp(2) - floor(variables.patches.size_x/2);
			b4 = variables.cp(2) + floor(variables.patches.size_x/2);
			if b1 < 1, b1 = 1; end
			if b2 > variables.data.max_y, b2 = variables.data.max_y; end
			if b3 < 1, b3 = 1; end
			if b4 > variables.data.max_x, b4 = variables.data.max_x; end
			H1 = wld(image_data{i}(b1:b2,b3:b4), variables.feature.num_orientation, variables.feature.num_excitation);
			features(:,i) = H1;
			%printing(variables.printing, variables.feature.name, i, inum);
		end
	end
end

