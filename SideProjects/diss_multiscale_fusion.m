function results = diss_multiscale_fusion(input1, input2, input3, input4, input5, input6)
% input(struct):	Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%					list of required fields.
%					
% input(string):	Filepath to xml file that will be parsed into a struct.
%
% results:			See Section 3.0 for list of outputs.

	start_t = tic;	
	if isstruct(input1)
		experiment1 = input1;
	elseif isempty(input1) == 0
		experiment1 = xml_read(input1);
    end
	
	if isstruct(input2)
		experiment2 = input2;
	elseif isempty(input2) == 0
		experiment2 = xml_read(input2);
    end
	
	if isstruct(input3)
		experiment3 = input3;
	elseif isempty(input3) == 0
		experiment3 = xml_read(input3);
    end
	
	if isstruct(input4)
		experiment4 = input4;
	elseif isempty(input4) == 0
		experiment4 = xml_read(input4);
    end
	
	if isstruct(input5)
		experiment5 = input5;
	elseif isempty(input5) == 0
		experiment5 = xml_read(input5);
    end
	
	if isstruct(input6)
		experiment6 = input6;
	elseif isempty(input6) == 0
		experiment6 = xml_read(input6);
    end
   
    
	fprintf('###Beginning %s\n', experiment1.id);
	fprintf('Image size: %d X %d\n', experiment1.variables.data.max_y, experiment1.variables.data.max_x);
	fprintf('Patch size: %d X %d\n', experiment1.variables.patches.size_y, experiment1.variables.patches.size_x);
	fprintf('Patch type: %s\n', experiment1.variables.patches.type);
	fprintf('Distance:   %s\n', experiment1.variables.distance);
	fprintf('Feature:    %s\n', experiment1.variables.feature.name);
	fprintf('Distance:   %s\n', experiment1.output.distance_mat);
	fprintf('Results:    %s\n', experiment1.output.results_mat);
		
	
	% load distance matricies
	% results = load([experiment1.output.resultsdir w1id 'distance.mat'], 'distances'); distances = results.distances; clearvars results;
	% results = load([experiment1.output.resultsdir w3id 'distance.mat'], 'distances'); distances = distances + results.distances; clearvars results;
	% results = load([experiment1.output.resultsdir w5id 'distance.mat'], 'distances'); distances = distances + results.distances; clearvars results;
	% results = load([experiment1.output.resultsdir w7id 'distance.mat'], 'distances'); distances = distances + results.distances; clearvars results;
	% results = load([experiment1.output.resultsdir w9id 'distance.mat'], 'distances'); distances = distances + results.distances; clearvars results;
	
	
	weights = [1.0 1.0 1.0 1.0 1.0 1.0];
	subregion = 1;
	for k = 1:2
		for j=1:length(weights)
			t2 = tic;
			Rank1s=zeros(21,1);
			for i=0:20
				results = load([experiment1.output.resultsdir experiment1.id 'distance.mat'], 'distances'); 
				if j == 1
					distances = (i*.25)*results.distances{subregion};
				else
					distances = (weights(1)*results.distances{subregion});
				end
				results = load([experiment2.output.resultsdir experiment2.id 'distance.mat'], 'distances');
				if j == 2
					distances = distances + (i*.25)*results.distances{subregion};
				else
					distances = distances + (weights(2)*results.distances{subregion});
				end
				results = load([experiment3.output.resultsdir experiment3.id 'distance.mat'], 'distances');
				if j == 3
					distances = distances + (i*.25)*results.distances{subregion};
				else
					distances = distances + (weights(3)*results.distances{subregion});
				end
				results = load([experiment4.output.resultsdir experiment4.id 'distance.mat'], 'distances');
				if j == 4
					distances = distances + (i*.25)*results.distances{subregion};
				else
					distances = distances + (weights(4)*results.distances{subregion});
				end
				
				if isempty(input5) == 0
					results = load([experiment5.output.resultsdir experiment5.id 'distance.mat'], 'distances'); 
					if j == 5
						distances = distances + (i*.25)*results.distances{subregion};
					else
						distances = distances + (weights(5)*results.distances{subregion});
					end
				end
				
				if isempty(input6) == 0
					results = load([experiment6.output.resultsdir experiment6.id 'distance.mat'], 'distances'); 
					if j == 6
						distances = distances + (i*.25)*results.distances{subregion};
					else
						distances = distances + (weights(5)*results.distances{subregion});
					end
					clearvars results;
				end
		
				CMC_rec_rates = compute_cmc(experiment1.id, distances, experiment1.variables, experiment1.input, experiment1.output);
				Rank1s(i+1) = CMC_rec_rates(1)*100;
				fprintf('[%3.2f w][%3.2f%%]', (i*.25), Rank1s(i+1));
			end
			fprintf('\n');
			[C,I] = max(Rank1s);
			weights(j) = (I-1)*0.25;
			time = toc(t2);
			fprintf('[%6.2f sec][%3.2f%% score][k=%d][j=%d] %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n\n',time,C,k,j,weights)
		end
	end
	
	% compute distance matrix
    t1 = tic;
    fprintf('##Computing distance matrix with %s\n', experiment1.variables.distance);
	results = load([experiment1.output.resultsdir experiment1.id 'distance.mat'], 'distances'); 
		distances = (weights(1)*results.distances{subregion});
	results = load([experiment2.output.resultsdir experiment2.id 'distance.mat'], 'distances');
		distances = distances + (weights(2)*results.distances{subregion});
	results = load([experiment3.output.resultsdir experiment3.id 'distance.mat'], 'distances');
		distances = distances + (weights(3)*results.distances{subregion});
	results = load([experiment4.output.resultsdir experiment4.id 'distance.mat'], 'distances');
		distances = distances + (weights(4)*results.distances{subregion});
	if isempty(input5) == 0
	results = load([experiment5.output.resultsdir experiment5.id 'distance.mat'], 'distances'); 
		distances = distances + (weights(5)*results.distances{subregion});
	end
	if isempty(input6) == 0
	results = load([experiment6.output.resultsdir experiment6.id 'distance.mat'], 'distances'); 
		distances = distances + (weights(6)*results.distances{subregion});
	end
	clearvars results;
	s=whos('distances');
	fprintf('Produced distance matrix [%d gallery images x %d probe images] [%6.2f MB]\n', size(distances,1), size(distances,2), s.bytes/1000000);
	%if ~isempty(experiment1.output.distance_mat)
	%	save([experiment1.output.resultsdir experiment1.id experiment1.output.distance_mat], 'distances', '-v7.3');
	%end
    toc(t1)
	fprintf('\n');
	
	% evaluate the experiment1
	t1 = tic;
	fprintf('##Evaluating the experiment1\n');
	% compute CMC
	fprintf('1. Computing CMC\n');
	t2 = tic;
	results.CMC_rec_rates = compute_cmc(experiment1.id, distances, experiment1.variables, experiment1.input, experiment1.output);
	toc(t2)
	
	fprintf('2. Split true and false scores\n');
	t2 = tic;
	[true_scores false_scores] = split_true_false_scores(experiment1.id, distances, experiment1.variables, experiment1.input, experiment1.output);
	clearvars distances;
	toc(t2)
	
	% compute ROC
	fprintf('3. Computing ROC\n');
	t2 = tic;
	[results.ROC_ver_rate, results.ROC_miss_rate, results.ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment1.output.roc_resolution);
	toc(t2)
	
	toc(t1)
	fprintf('\n');
	
	% print points on each curve to a file
	if experiment1.output.make_files == 1
		fprintf('5. Print points to file\n');
		t2 = tic;
		file = fopen([experiment1.output.resultsdir experiment1.id '_det.txt'], 'w');
		for i = 1:length(results.ROC_ver_rate)
		   fprintf(file, '%f,%f\n', results.ROC_miss_rate(i), 1-results.ROC_ver_rate(i));
		end
		fclose(file);

		file = fopen([experiment1.output.resultsdir experiment1.id '_roc.txt'], 'w');
		for i = 1:length(results.ROC_ver_rate)
		   fprintf(file, '%f,%f\n', results.ROC_miss_rate(i), results.ROC_ver_rate(i));
		end
		fclose(file);

		file = fopen([experiment1.output.resultsdir experiment1.id '_cmc.txt'], 'w');
		for i = 1:length(results.CMC_rec_rates)
		   fprintf(file, '%f\n', results.CMC_rec_rates(i));
		end
		fclose(file);
		toc(t2)	
	end
	
	% plotting results
	if experiment1.output.make_figures == 1
		t1 = tic;
		fprintf('##Plotting the results\n');
		
		hdet = plot_det(1-results.ROC_ver_rate, results.ROC_miss_rate, results.ROC_rates_and_threshs.EER_er, experiment1.variables.feature.name, experiment1.output.plot_color, experiment1.output.plot_linewidth, 1);
		saveas(hdet,[experiment1.output.resultsdir experiment1.id '_det.fig']);
		close(hdet);
		
		hcmc = plot_cmc(results.CMC_rec_rates, experiment1.variables.feature.name, experiment1.output.plot_color, experiment1.output.plot_linewidth, 2);
		saveas(hcmc,[experiment1.output.resultsdir experiment1.id '_cmc.fig']);
		close(hcmc);
		
		hroc = plot_roc(results.ROC_ver_rate, results.ROC_miss_rate, results.ROC_rates_and_threshs.EER_er, experiment1.variables.feature.name, experiment1.output.plot_color, experiment1.output.plot_linewidth, 3);
		saveas(hroc,[experiment1.output.resultsdir experiment1.id '_roc.fig']);
		close(hroc);
		
		toc(t1)
		fprintf('\n');
	end
    end_t = toc(start_t);
	
	% save results matrix
	if ~isempty(experiment1.output.results_mat)
		if(~exist(experiment1.output.resultsdir, 'dir'))
		   mkdir(experiment1.output.resultsdir);
		end
		save([experiment1.output.resultsdir '_fusion_' experiment1.id experiment1.output.results_mat], 'results');
	end
    
    % print out results
	fprintf('\nexperiment1al Results: %s\n', experiment1.id);
    fprintf('=============================================================\n');
    fprintf('Identification experiment1s:\n');
    fprintf('Rank-1   recognition rate: %7.4f%%\n', results.CMC_rec_rates(1)*100);
	if length(results.CMC_rec_rates) > 10
		fprintf('Rank-10  recognition rate: %7.4f%%\n', results.CMC_rec_rates(10)*100);
	else
		fprintf('Rank-%d recognition rate: %7.4f%%\n', length(results.CMC_rec_rates), results.CMC_rec_rates(length(results.CMC_rec_rates))*100);
	end
	if length(results.CMC_rec_rates) > 100
		fprintf('Rank-100 recognition rate: %7.4f%%\n', results.CMC_rec_rates(100)*100);
	else
		fprintf('Rank-%d recognition rate: %7.4f%%\n', length(results.CMC_rec_rates), results.CMC_rec_rates(length(results.CMC_rec_rates))*100);
	end
    fprintf('Verification/Authentication experiment1s:\n');
    fprintf('EER: %7.4f%%\n', results.ROC_rates_and_threshs.EER_er*100);
	fprintf('Dprime: %7.4f\n', results.ROC_rates_and_threshs.dprime);
    fprintf('VR at 1%% FAR:    %7.4f%%\n', results.ROC_rates_and_threshs.VER_1FAR_ver*100);
    fprintf('VR at 0.1%% FAR:  %7.4f%%\n', results.ROC_rates_and_threshs.VER_01FAR_ver*100);
    fprintf('VR at 0.01%% FAR: %7.4f%%\n', results.ROC_rates_and_threshs.VER_001FAR_ver*100);
	fprintf('Recap:\n');
	%fprintf('%7.4f\t%7.4f\t%7.4f\n',results.CMC_rec_rates(1)*100,results.ROC_rates_and_threshs.EER_er*100,results.ROC_rates_and_threshs.VER_01FAR_ver*100);
	fprintf('%7.4f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n',results.CMC_rec_rates(1)*100,weights);
    fprintf('=============================================================\n');
    fprintf('Elapsed time: %f seconds\n', end_t);
end

function ret = checkSameList(gallery, probe)
	%%% Scan in image filenames %%%
	fid = fopen(gallery);
	glist = textscan(fid, '%s');
	fclose(fid);
	glist = glist{1};
	gnum = size(glist,1);
	
	%%% Scan in image filenames %%%
	fid = fopen(probe);
	plist = textscan(fid, '%s');
	fclose(fid);
	plist = plist{1};
	pnum = size(plist,1);
	
	if pnum ~= gnum
		ret=0;
	elseif (sum(strcmp(plist, glist)) == gnum)
		ret=1;
	else
		ret=0;
	end
end
