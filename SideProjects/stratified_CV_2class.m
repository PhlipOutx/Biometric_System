%stratified k-fold cross validation (2 class)
%%%%%%%%%INPUTS
%features, features to be used for classification
%NEED to have these in the XML FILE
%variables.cross_validation.training_class, path to file with class per image
%                           .classes, how many classes present/training on
%                           .rounds, number of rounds for cross validation
%                           .type, what type of classifier 
%**(only 'svm' supported at the moment)
%                           .svm_command, command to pass to libsvm
%input.training, file path ot the imagelist
%     .training_id, file path to the subject id file
%%%%%%%%OUTPUTS
%k_test_lists, imagelists for each round
%k_performances, accuracy for each round
%predicted_class, predicted class for each round

function [k_test_lists, k_performances, predicted_class] = ...
    stratified_CV_2class(id, features, variables, input, output)

    %%% Scan in image filenames %%%
    fid = fopen(input.training);
    image_names = textscan(fid, '%s');
    fclose(fid);
    image_names = image_names{1};
    inum = size(image_names,1);

    %%% Scan in class ids %%%
    fid = fopen(variables.cross_validation.training_class);
    image_class = textscan(fid, '%s');
    fclose(fid);
    image_class = image_class{1};
    assert(inum == size(image_class, 1), 'Class ID size does not match imagelist size')
    
    %%% check number of unique classes %%%
    [image_class_val, ~, image_class_ids] = unique(image_class);%image_class_loc
    cnum=size(image_class_val,1);
    assert(cnum == variables.cross_validation.classes, ['Number of classes' ...
        'does not match the number present in file']);
    assert(cnum==2,'Not set up to handle more than 2 classes yet')
    %if cnum>2
    %    display(['WARNING: not sure on C-Class SVM yet, I take no ' ...
    %    'responsibility for what follows']);
    %end
    
    %%% find subject id for each image %%%
    if isempty(input.training_id)
        iids = cellfun(@(x) x(1:variables.nameprefix), image_names, 'UniformOutput', false);
    else
        fid = fopen(input.training_id);
        iids = textscan(fid, '%d');
        fclose(fid);
        iids = iids{1};
        assert(inum == size(iids,1), 'Subject IDs size does not match imagelist size')
    end
    
    %%% get unique subjects
    [~, subject_loc, subject_ids] = unique(iids);%subject_vals
    snum = size(subject_ids,1);
    
    %%% get class of unique subjects
    subject_class = image_class_ids(subject_loc);
    
    statistics=zeros(cnum,1);           %split of classes (totals)
    class_ndxs=cell(cnum,1);            %where people in each class are
    subjects=cell(cnum,1);              %order to use the subjects in each class
    for i = 1:cnum
        loc=find(subject_class==i);     %find everybody in class i
        tf=subject_class==i;
        num=sum(tf);                    %how many subjects are in class i
        statistics(i) = num;
        class_ndxs{i}=loc;              %store everbody that is in class i 
        subjects{i}=class_ndxs{i}(randperm(num));%random permutation of 
                                                 %everybody by class
    end
    
    clear tf loc i;
    
    %%% Start cross validation
    rounds=variables.cross_validation.rounds;
    k_test_lists  = cell(rounds,1);     %image test list for each round
    k_performances = zeros(rounds,1);   %accuracy for each round
    predicted_class = cell(rounds,1);   %predicted class each round
    test_size=floor(1/rounds*snum);     %number of subjects per round
    test_size_stats=round(statistics*test_size/snum);%by subject
    
    for k = 1:rounds
        %%% get indices to start and stop test set
        %%% to pull subjects into test set
        start_ndx=test_size_stats*(k-1)+1;
        stop_ndx=test_size_stats*k;
        if k==rounds
            stop_ndx=statistics;
        end
        
        %%% build test set
        test_ndx=[];
        for c=1:cnum
            %find all the image subject ids that match the test subjects
            %[tf loc]=ismember(subjects{c}(start_ndx(c):stop_ndx(c)),iids);
            %DON'T THINK THIS WILL WORK WITH MULTIPLE IMAGES per subject
            %try this instead
            for ndx=start_ndx(c):stop_ndx(c)
                loc=find(iids==subjects{c}(ndx));
                test_ndx=vertcat(test_ndx,loc);
            end
        end
        
        %%% get subset of features and class ids
        test_features=features(:,test_ndx);
        test_class=image_class_ids(test_ndx);
        
        %%% build training set from everybody not in the test set
        tf = ismember(1:inum, test_ndx);
        train_features=features(:,~tf);
        train_class=image_class_ids(~tf);
        
        %%% train and test according to classifier type
        switch lower(variables.cross_validation.type)
            case 'svm'
                 model = svmtrain(train_class, train_features', ...
                     variables.cross_validation.svm_command);
                [predicted_label, accuracy, ~] = svmpredict(test_class, test_features', model);
                predicted_class{k}=predicted_label;
                k_performances(k)=accuracy(1); 
            case 'nn'
                disp('ERROR: classifier not supported yet');
                break;
            otherwise
                disp('ERROR: what classifier did you want to use?');
                break;
        end;
        %get accuracy
        
        k_test_lists=image_names{test_ndx};
        
        
    end
    
