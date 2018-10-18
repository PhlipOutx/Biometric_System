function variable_selection(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if isstruct(input)
        experiment = input;
    else
        % parse xml experiment file
        experiment = xml_read(input);
    end
    
    if strcmp(experiment.variables.feature.name,'LBP') || strcmp(experiment.variables.feature.name,'HOG') || strcmp(experiment.variables.feature.name,'LPQ') || strcmp(experiment.variables.feature.name,'WLD')
        variable_selection_local(experiment);
    elseif strcmp(experiment.variables.feature.name,'PCA') || strcmp(experiment.variables.feature.name,'LDA') || strcmp(experiment.variables.feature.name,'KFA') || strcmp(experiment.variables.feature.name,'KPCA') || strcmp(experiment.variables.feature.name,'ICA')
        variable_selection_global(experiment);
    elseif strcmp(experiment.variables.feature.name,'SIFT') || strcmp(experiment.variables.feature.name,'SURF')
        variable_selection_keypoint(experiment);
    end
    
end
    
function variable_selection_local(experiment)

    % allocate
    patch_size = [10, 15, 20, 25, 30, 35, 40, 45, 50];
    if strcmp(experiment.variables.feature.name, 'LPQ')==1
        patch_size = [15, 20, 25, 30, 35, 40, 45, 50];
    end
    data_size = [2;3;4;5;6;7;8;9;10]*patch_size;
    rank1s = zeros(size(data_size));
    eers = ones(size(data_size))*100;
    vrs = zeros(size(data_size));	
    
    for i = 1:size(rank1s,1)
        for j = 1:size(rank1s,2)
            fprintf('\n\n BEGINNING i=%d j=%d patch=%d data=%d \n\n', i, j, patch_size(j), data_size(i,j));
            experiment.variables.patches.size_x = patch_size(j);
            experiment.variables.patches.size_y = patch_size(j);
            experiment.variables.data.max_x = data_size(i,j);
            experiment.variables.data.max_y = data_size(i,j);
        
            results = biometric_system(experiment);
            rank1s(i,j) = results.CMC_rec_rates(1)*100;
            eers(i,j) = results.ROC_rates_and_threshs.EER_er*100;
            vrs(i,j) = results.ROC_rates_and_threshs.VER_01FAR_ver*100;
        end
        if experiment.variables.select_type == 1 && i > 1
            [C1 I1] = max(rank1s,[],2);
            [C2 I2] = max(C1);
            if C2 ~= C1(i) && C2 ~= C1(i-1), break, end
        elseif experiment.variables.select_type == 2 && i > 1
            [C1 I1] = min(eers,[],2);
            [C2 I2] = min(C1);
            if C2 ~= C1(i) && C2 ~= C1(i-1), break, end
        elseif experiment.variables.select_type == 3 && i > 1
            [C1 I1] = max(vrs,[],2);
            [C2 I2] = max(C1);
            if C2 ~= C1(i) && C2 ~= C1(i-1), break, end
        end
    end
    
    fprintf('\n###RESULTS###\n');
    fprintf('\n###Rank-1s###\n');
    disp(rank1s);
    fprintf('\n###EERs###\n');
    disp(eers);
    fprintf('\n###VR at 0.1%% FAR###\n');
    disp(vrs);

    if experiment.variables.select_type == 1
        [C1 I1] = max(rank1s,[],2);
        [C2 I2] = max(C1);
        fprintf('%d block size\n', patch_size(I1(I2)));
        fprintf('%d x %d image size\n', data_size(I2,I1(I2)), data_size(I2,I1(I2)));
        fprintf('Rank-1: %f\n', rank1s(I2,I1(I2)));
    elseif experiment.variables.select_type == 2
        [C1 I1] = min(eers,[],2);
        [C2 I2] = min(C1);
        fprintf('%d block size\n', patch_size(I1(I2)));
        fprintf('%d x %d image size\n', data_size(I2,I1(I2)), data_size(I2,I1(I2)));
        fprintf('EER: %f\n', eers(I2,I1(I2)));
    elseif experiment.variables.select_type == 3
        [C1 I1] = max(vrs,[],2);
        [C2 I2] = max(C1);
        fprintf('%d block size\n', patch_size(I1(I2)));
        fprintf('%d x %d image size\n', data_size(I2,I1(I2)), data_size(I2,I1(I2)));
        fprintf('VR at 0.1%% FAR: %f\n', vrs(I2,I1(I2)));
    end
end
    
    
function variable_selection_global(experiment)

    max_vec = experiment.variables.feature.vectors;
    use_x = experiment.variables.data.max_x;
    use_y = experiment.variables.data.max_y;
    
    % allocate
    ret_vecs = [5,10,25,50,75,100,125,150,175,200,250,300,350,400,500,600,700,800,1000,1500,2000,2500,3000,4000,5000,7500,10000];
    ret_vecs = ret_vecs(ret_vecs < max_vec);
    ret_vecs = [ret_vecs,max_vec-1];
    rank1s = zeros(size(ret_vecs));
    eers = ones(size(ret_vecs))*100;
    vrs = zeros(size(ret_vecs));
    
    for i = 1:numel(rank1s)
        fprintf('\n\n BEGINNING i=%d ret_vecs=%d\n\n', i, ret_vecs(i));
        experiment.variables.feature.vectors = ret_vecs(i);
        
        results = biometric_system(experiment);
        rank1s(i) = results.CMC_rec_rates(1)*100;
        eers(i) = results.ROC_rates_and_threshs.EER_er*100;
        vrs(i) = results.ROC_rates_and_threshs.VER_01FAR_ver*100;
        if experiment.variables.select_type == 1
            [C1 I1] = max(rank1s);
            if I1 < i-5, break, end
        elseif experiment.variables.select_type == 2
            [C1 I1] = min(eers);
            if I1 < i-5, break, end
        elseif experiment.variables.select_type == 3
            [C1 I1] = max(vrs);
            if I1 < i-5, break, end
        end
    end
    
    fprintf('\n###RESULTS###\n');
    fprintf('\n###Rank-1s###\n');
    disp(rank1s);
    fprintf('\n###EERs###\n');
    disp(eers);
    fprintf('\n###VR at 0.1%% FAR###\n');
    disp(vrs);

    if experiment.variables.select_type == 1
        [C1 I1] = max(rank1s);
        fprintf('%d x %d image size\n', use_x, use_y);
        fprintf('Vectors to retain: %d\n', ret_vecs(I1));
        fprintf('Rank-1: %f\n', rank1s(I1));
    elseif experiment.variables.select_type == 2
        [C1 I1] = min(eers);
        fprintf('%d x %d image size\n', use_x, use_y);
        fprintf('Vectors to retain: %d\n', ret_vecs(I1));
        fprintf('EER: %f\n', eers(I1));
    elseif experiment.variables.select_type == 3
        [C1 I1] = max(vrs);
        fprintf('%d x %d image size\n', use_x, use_y);
        fprintf('Vectors to retain: %d\n', ret_vecs(I1));
        fprintf('VR at 0.1%% FAR: %f\n', vrs(I1));
    end
end
    
function variable_selection_keypoint(experiment)

    % allocate
    data_size = 60:10:800;
    rank1s = zeros(size(data_size));	
    eers = ones(size(data_size))*100;	
    vrs = zeros(size(data_size));

    for i = 1:numel(rank1s)
        fprintf('\n\n BEGINNING i=%d data_size=%d \n\n', i, data_size(i));
        experiment.variables.data.max_x = data_size(i);
        experiment.variables.data.max_y = data_size(i);
        
        results = biometric_system(experiment);
        rank1s(i) = results.CMC_rec_rates(1)*100;
        eers(i) = results.ROC_rates_and_threshs.EER_er*100;
        vrs(i) = results.ROC_rates_and_threshs.VER_01FAR_ver*100;
        if experiment.variables.select_type == 1
            [C1 I1] = max(rank1s);
            if I1 < i-5, break, end
        elseif experiment.variables.select_type == 2
            [C1 I1] = min(eers);
            if I1 < i-5, break, end
        elseif experiment.variables.select_type == 3
            [C1 I1] = max(vrs);
            if I1 < i-5, break, end
        end
    end
    
    fprintf('\n###RESULTS###\n');
    fprintf('\n###Rank-1s###\n');
    disp(rank1s);
    fprintf('\n###EERs###\n');
    disp(eers);
    fprintf('\n###VR at 0.1%% FAR###\n');
    disp(vrs);

    if experiment.variables.select_type == 1
        [C1 I1] = max(rank1s);
        fprintf('%d x %d image size\n', data_size(I1), data_size(I1));
        fprintf('Rank-1: %f\n', rank1s(I1));
    elseif experiment.variables.select_type == 2
        [C1 I1] = min(eers);
        fprintf('%d x %d image size\n', data_size(I1), data_size(I1));
        fprintf('EER: %f\n', eers(I1));
    elseif experiment.variables.select_type == 3
        [C1 I1] = max(vrs);
        fprintf('%d x %d image size\n', data_size(I1), data_size(I1));
        fprintf('VR at 0.1%% FAR: %f\n', vrs(I1));
    end
end
    