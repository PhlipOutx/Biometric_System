function results = SURF_image(input)
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
    images = textscan(fid, '%s');
    fclose(fid);
    images = images{1};
    inum = size(images,1);
    
    %%% Create empty matrix the size of an image %%%
    image_data = cell(inum,1);
    
    fprintf('1. Computing mean image from %d images\n', inum);
    for i = 1:inum
        % preprocess image
        temp = char(images(i));
        I = imread([experiment.input.datadir temp]);
        image_data{i} = ipreproc(I, experiment.variables);
        
        printing(experiment.variables.printing, experiment.variables.feature.name, i, inum);
    end
    fprintf('\n');
    
    mean_image=zeros(size(image_data{1}));
    for k = 1:length(image_data)
       mean_image = mean_image+double(image_data{k});
    end
    mean_image = mean_image/double(length(image_data));
    mean_image=uint8(round(mean_image));
    imwrite(mean_image,[experiment.id 'mean.pgm'],'pgm');
    
    end_t = toc(start_t);
    fprintf('=============================================================\n');
    fprintf('Elapsed time: %f seconds\n', end_t);
end
