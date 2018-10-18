function results = cross_validation_experiment(input)
% input(struct):    Matlab structure containing fields that hold variables
%                   used in a cross validation experiment. See Section 2.0 for
%                   list of required fields. 
%**
%(Needs input.training and input.training_ids, but no probe or gallery files)
%Also needs
%(
%variables.cross_validation.training_class, path to file with class per image
%                          .classes, how many classes present/training on
%                          .rounds, number of rounds for cross validation
%                          .type, what type of classifier 
%**(only 'svm' supported at the moment)
%                          .svm_command, command to pass to libsvm
%)in addition to features inputs
%					
% input(string):    Filepath to xml file that will be parsed into a struct.
%
% results:          See Section 3.0 for list of outputs. (only performance)

    start_t = tic;	
    if isstruct(input)
        experiment = input;
    else
        % parse xml experiment file
        experiment = xml_read(input);
    end
    
    %%% extract training features %%%
    t1 = tic;
    fprintf('##Extracting features\n');
    fprintf('#From images in %s\n', experiment.input.training);
    fprintf('#With training set in %s\n', experiment.input.training);
    training_features = extract_features(experiment.id, experiment.variables, ...
                               experiment.input, experiment.output, experiment.input.training, ...
                               experiment.input.datadir, experiment.output.training_mat, ...
                               experiment.output.resultsdir);
    s=whos('training_features');
    fprintf('Produced training feature matrix [%d features x %d images] [%6.2f MB]\n', ...
    size(training_features,1), size(training_features,2), s.bytes/1000000);
    toc(t1)
    fprintf('\n');
    
    %%% go into cross validation experiment %%%
    [test_lists, accuracy, predicted_class]=stratified_CV_2class(experiment.id, training_features, ...
                                                        experiment.variables, experiment.input, ...
                                                        experiment.output);
                                                    
    %%% print results %%%
    end_t = toc(start_t);
    fprintf('Cross validation accuracy:\t%5.2f\n', sum(accuracy)/experiment.variables.cross_validation.rounds);
    fprintf('Elapsed time: %f seconds\n', end_t);
    
    %%% set results %%%
    results.cv_accuracy=sum(accuracy)/experiment.variables.cross_validation.rounds;
    results.accuracy = accuracy;
    results.predicted_class = predicted_class;
end
