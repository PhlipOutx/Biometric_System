function results = Keypoints_image(input)
% input(struct):    Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%                   list of required fields.
%					
% input(string):    Filepath to xml file that will be parsed into a struct.
%
% results:          See Section 3.0 for list of outputs.

    start_t = tic;	
    if isstruct(input)
        experiment = input;
    else
        % parse xml experiment file
        experiment = xml_read(input);
    end
    
    fprintf('###Beginning %s\n', experiment.id);
    
    % extract gallery features
    fprintf('##Extracting gallery features\n');
    fprintf('#From images in %s\n', experiment.input.gallery);

    %%% Scan in image filenames %%%
    fid = fopen(experiment.input.gallery);
    galleries = textscan(fid, '%s');
    fclose(fid);
    galleries = galleries{1};
    gnum = size(galleries,1);
    
    fid = fopen(experiment.input.probe);
    probes = textscan(fid, '%s');
    fclose(fid);
    probes = probes{1};
    pnum = size(probes,1);
    
    %%% Create empty matrix the size of an image %%%
    temp = char(galleries(1));
    I = imread([experiment.input.datadir temp]);
    I = ipreproc(I, experiment.variables);
    map = zeros(size(I,1), size(I,2));
    num_keypoints = 0;

    if strcmp(experiment.variables.feature.name,'SIFT')  	
        fprintf('1. Computing keypoint locations from %d images\n', gnum);
        for i = 1:gnum
            % preprocess image
            temp = char(galleries(i));
            I = imread([experiment.input.datadir temp]);
            I = ipreproc(I, experiment.variables);

            % SIFT calculations
            IS = single(I);
            [F1,~] = vl_sift(IS);
            for j = 1:size(F1,2);
                map(round(F1(2,j)), round(F1(1,j))) = map(round(F1(2,j)), round(F1(1,j))) + 1;
            end
            num_keypoints = num_keypoints + size(F1,2);
            
            printing(experiment.variables.printing, experiment.variables.feature.name, i, gnum);
        end
        fprintf('\n');
        fprintf('2. Computing keypoint locations from %d images\n', pnum);
        for i = 1:pnum
            % preprocess image
            temp = char(probes(i));
            I = imread([experiment.input.datadir temp]);
            I = ipreproc(I, experiment.variables);

            % SIFT calculations
            IS = single(I);
            [F1,~] = vl_sift(IS);
            for j = 1:size(F1,2);
                map(round(F1(2,j)), round(F1(1,j))) = map(round(F1(2,j)), round(F1(1,j))) + 1;
            end
            num_keypoints = num_keypoints + size(F1,2);
            
            printing(experiment.variables.printing, experiment.variables.feature.name, i, pnum);
        end

        % save results matrix
        map2 = double(map) / double(max(max(map)));
        map2 = uint8(round(map2 * 255.0));
        imwrite(map2,[experiment.id '_keypoints.ppm'],'ppm');
        
        set(gcf, 'Visible', 'off')
        h = imagesc(map);    
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        colorbar;
        saveas(h,[experiment.id '_heatmap.ppm'],'ppm')
        
        % print out results
        fprintf('\nSIFT_Image Results: %s\n', experiment.id);
        fprintf('=============================================================\n');
        fprintf('Found %d keypoints\n', num_keypoints);
        fprintf('From %d images\n', pnum+gnum);
        fprintf('Largest number of keypoints in one pixel: %d\n', max(max(map)));
        fprintf('Least number of keypoints in one pixel: %d\n', min(min(map)));
        fprintf('Number of pixels with zero keypoints: %d\n', sum(map(:)==0));
        end_t = toc(start_t);
        fprintf('=============================================================\n');
        fprintf('Elapsed time: %f seconds\n', end_t);
    elseif strcmp(experiment.variables.feature.name,'SURF')
        Options.verbose=false;
        fprintf('1. Computing keypoint locations from %d images\n', gnum);
        for i = 1:gnum
            % preprocess image
            temp = char(galleries(i));
            I = imread([experiment.input.datadir temp]);
            I = ipreproc(I, experiment.variables);

            % SURF calculations
            vec = OpenSurf(I, Options);
            for j = 1:size(vec,2);
                map(round(vec(j).y), round(vec(j).x)) = map(round(vec(j).y), round(vec(j).x)) + 1;
            end
            num_keypoints = num_keypoints + size(vec,2);
            
            printing(experiment.variables.printing, experiment.variables.feature.name, i, gnum);
        end
        fprintf('\n');
        fprintf('2. Computing keypoint locations from %d images\n', pnum);
        for i = 1:pnum
            % preprocess image
            temp = char(probes(i));
            I = imread([experiment.input.datadir temp]);
            I = ipreproc(I, experiment.variables);

            % SURF calculations
            vec = OpenSurf(I, Options);
            for j = 1:size(vec,2);
                map(round(vec(j).y), round(vec(j).x)) = map(round(vec(j).y), round(vec(j).x)) + 1;
            end
            num_keypoints = num_keypoints + size(vec,2);
            
            printing(experiment.variables.printing, experiment.variables.feature.name, i, pnum);
        end

        % save results matrix
        map2 = double(map) / double(max(max(map)));
        map2 = uint8(round(map2 * 255.0));
        imwrite(map2,[experiment.id '_keypoints.ppm'],'ppm');
        
        set(gcf, 'Visible', 'off')
        h = imagesc(map);    
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        colorbar;
        saveas(h,[experiment.id '_heatmap.ppm'],'ppm')
        
        % print out results
        fprintf('\nSIFT_Image Results: %s\n', experiment.id);
        fprintf('=============================================================\n');
        fprintf('Found %d keypoints\n', num_keypoints);
        fprintf('From %d images\n', pnum+gnum);
        fprintf('Largest number of keypoints in one pixel: %d\n', max(max(map)));
        fprintf('Least number of keypoints in one pixel: %d\n', min(min(map)));
        fprintf('Number of pixels with zero keypoints: %d\n', sum(map(:)==0));
        end_t = toc(start_t);
        fprintf('=============================================================\n');
        fprintf('Elapsed time: %f seconds\n', end_t);
    end
end
