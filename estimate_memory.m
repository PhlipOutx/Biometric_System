function estimate_memory(input)
%input
% input(struct):    Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%                   list of required fields.
%					
% input(string):    Filepath to xml file that will be parsed into a struct.
%
%output
% Prints estimated memory to be used during each part of the experiments.
% Ensure that there is enough memory on the machine running the experiment to
% accomidate the experiment. 

    if isstruct(input)
        experiment = input;
    else
        % parse xml experiment file
        experiment = xml_read(input);
    end
    
    % get available memory
    %[user,sys] = memory;
    
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
    
    if exist(experiment.input.training,'file')
        fid = fopen(experiment.input.training);
        training_names = textscan(fid, '%s');
        fclose(fid);
        training_names = training_names{1};
        tnum = size(training_names,1);
    end
    
    % other vars
    pixels = experiment.variables.data.max_x * experiment.variables.data.max_y;
    patches = (experiment.variables.data.max_x / experiment.variables.patches.size_x) * (experiment.variables.data.max_y / experiment.variables.patches.size_y);
    
    % calculate
    if strcmp(experiment.variables.feature.name,'SIFT')
        gdur  = gnum*128*100*8;
        gsize = gdur;
        pdur  = gsize + pnum*128*100*8;
        psize = pnum*128*100*8;
    elseif strcmp(experiment.variables.feature.name,'SURF')
        gdur  = gnum*64*100*8;
        gsize = gdur;
        pdur  = gsize + pnum*64*100*8;
        psize = pnum*64*100*8;
    elseif strcmp(experiment.variables.feature.name,'LBP')
        mapping = getmapping(experiment.variables.feature.samples, experiment.variables.feature.type);
        gdur  = gnum*patches*mapping.num*8;
        gsize = gdur;
        pdur  = gsize + pnum*patches*mapping.num*8;
        psize = pnum*patches*mapping.num*8;
    elseif strcmp(experiment.variables.feature.name,'HOG')
        gdur  = gnum*patches*experiment.variables.feature.bins*8;
        gsize = gdur;
        pdur  = gsize + pnum*patches*experiment.variables.feature.bins*8;
        psize = pnum*patches*experiment.variables.feature.bins*8;
    elseif strcmp(experiment.variables.feature.name,'LPQ')
        gdur  = gnum*patches*256*8;
        gsize = gdur;
        pdur  = gsize + pnum*patches*256*8;
        psize = pnum*patches*256*8;
    elseif strcmp(experiment.variables.feature.name,'WLD')
        gdur  = gnum*patches*experiment.variables.feature.num_orientation*experiment.variables.feature.num_excitation*8;
        gsize = gdur;
        pdur  = gsize + pnum*patches*experiment.variables.feature.num_orientation*experiment.variables.feature.num_excitation*8;
        psize = pnum*patches*experiment.variables.feature.num_orientation*experiment.variables.feature.num_excitation*8;
    elseif strcmp(experiment.variables.feature.name,'PCA')
        gdur1 = pixels*tnum*8 + 2*tnum*tnum*5*8 + pixels*experiment.variables.feature.vectors;
        gdur2 = pixels*gnum*8 + pixels*experiment.variables.feature.vectors*8 + tnum*experiment.variables.feature.vectors*8 + gnum*experiment.variables.feature.vectors*8;
        gdur  = max([gdur1 gdur2]);
        gsize = gnum*experiment.variables.feature.vectors*8;
        pdur  = gsize + pixels*pnum*8 + pixels*experiment.variables.feature.vectors*8 + tnum*experiment.variables.feature.vectors*8 + pnum*experiment.variables.feature.vectors*8;
        psize = pnum*experiment.variables.feature.vectors*8;
    elseif strcmp(experiment.variables.feature.name,'LDA')
        gdur1 = pixels*tnum*8 + tnum*tnum*5*8;
        gdur2 = pixels*tnum*5*8 + tnum*tnum*2*8;
        gdur3 = pixels*tnum*3*8 + tnum*tnum*7*8;
        gdur4 = pixels*tnum*4*8 + tnum*tnum*5*8;
        gdur5 = pixels*gnum*8 + pixels*tnum*8 + tnum*experiment.variables.feature.vectors*8 + gnum*experiment.variables.feature.vectors*8;
        gdur  = max([gdur1 gdur2 gdur3 gdur4 gdur5]);
        gsize = gnum*experiment.variables.feature.vectors*8;
        pdur  = gsize + pixels*pnum*8 + pixels*experiment.variables.feature.vectors*8 + tnum*experiment.variables.feature.vectors*8 + pnum*experiment.variables.feature.vectors*8;
        psize = pnum*experiment.variables.feature.vectors*8;
    elseif strcmp(experiment.variables.feature.name,'KFA')
        gdur1 = pixels*tnum*2*8 + tnum*tnum*10*8;
        gdur2 = pixels*tnum*2*8 + tnum*tnum*6*8 + tnum*experiment.variables.feature.vectors*3*8 + experiment.variables.feature.vectors*experiment.variables.feature.vectors*8;
        gdur3 = pixels*tnum*8 + tnum*tnum*2*8 + tnum*experiment.variables.feature.vectors*2*8 + pixels*gnum*8 + gnum*tnum*8 + gnum*experiment.variables.feature.vectors*8;
        gdur  = max([gdur1 gdur2]);
        gsize = gnum*experiment.variables.feature.vectors*8;
        pdur  = gsize + pixels*tnum*8 + tnum*tnum*2*8 + tnum*experiment.variables.feature.vectors*2*8 + pixels*pnum*8 + pnum*tnum*8 + pnum*experiment.variables.feature.vectors*8; 
        psize = pnum*experiment.variables.feature.vectors*8;
    elseif strcmp(experiment.variables.feature.name,'KPCA')
        gdur1 = pixels*tnum*2*8 + tnum*tnum*6*8 + tnum*experiment.variables.feature.vectors*2*8;
        gdur2 = pixels*tnum*8 + tnum*tnum*2*8 + tnum*experiment.variables.feature.vectors*2*8 + pixels*gnum*8 + gnum*tnum*8 + gnum*experiment.variables.feature.vectors*8;
        gdur  = max([gdur1 gdur2]);
        gsize = gnum*experiment.variables.feature.vectors*8;
        pdur  = gsize + pixels*tnum*8 + tnum*tnum*2*8 + tnum*experiment.variables.feature.vectors*2*8 + pixels*pnum*8 + pnum*tnum*8 + pnum*experiment.variables.feature.vectors*8;
        psize = pnum*experiment.variables.feature.vectors*8;
    end
    dist  = gnum*pnum*8;
    eval  = gnum*pnum*8*2;
    
        
    
    % print
    %fprintf('You have [%9.2f MB] of available memeory\n', (sys.PhysicalMemory.Available-200000000)/1000000);
    fprintf('There are [%d] gallery images and [%d] probe images\n', gnum, pnum);
    fprintf('[%u] similarity scores will be calculated\n', uint64(gnum*pnum));
    fprintf('You need an estimated [%9.2f MB] of memory during gallery features\n', (gdur)/1000000);
    fprintf('You need an estimated [%9.2f MB] of memory for gallery features\n',    (gsize)/1000000);
    fprintf('You need an estimated [%9.2f MB] of memory during probe features\n',   (pdur)/1000000);
    fprintf('You need an estimated [%9.2f MB] of memory for probe features\n',      (psize)/1000000);
    fprintf('You need an estimated [%9.2f MB] of memory during distance matrix\n',  (gsize+psize+dist)/1000000);
    fprintf('You need an estimated [%9.2f MB] of memory for distance matrix\n',     (dist)/1000000);
    fprintf('You need an estimated [%9.2f MB] of memory during evaluation\n',       (eval)/1000000);
    
end