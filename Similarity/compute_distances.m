function distances = compute_distances(gallery_features, probe_features, id, variables, input, output)

	% check for already created subspace training matrix
    if ~isempty(output.distance_mat) && exist([output.resultsdir id output.distance_mat],'file')
        distances = load([output.resultsdir id output.distance_mat], 'distances');
        distances = distances.distances;
		
		% if strcmp(variables.feature.name,'SIFT')
			% fprintf('Continuing SIFT distance\n');
			% gnum = size(gallery_features,2);
			% pnum = size(probe_features,2);
			% for i = 1:pnum
				% if sum(distances(:,i)==0) > 1
					% ti = probe_features{i};
					% for j = 1:gnum
						% tj = gallery_features{j};
						% [matches,scores] = vl_ubcmatch(ti,tj);
						% if length(scores) == 0
							% tempd = NaN;
						% else
							% tempd = (sum(scores)/size(scores,2))/size(scores,2);
						% end
						% distances(j,i) = tempd;
					% end
					% printing(variables.printing, 'SIFT distance', i, pnum);
					% if mod(uint32(pnum/50),i) == 0
						% save([output.resultsdir id output.distance_mat], 'distances', '-v7.3');
					% end
				% end
			% end
			% distances(isnan(distances)) = max(max(distances))+1;
		% elseif strcmp(variables.feature.name,'SURF')
			% fprintf('Continuing SURF distance\n');
			% gnum = size(gallery_features,2);
			% pnum = size(probe_features,2);
			% for i = 1:pnum
				% if sum(distances(:,i)==0) > 1
					% ti = probe_features{i};
					% for j = 1:gnum
						% tj = gallery_features{j};
						% [matches,scores] = vl_ubcmatch(ti,tj);
						% if length(scores) == 0
							% tempd = NaN;
						% else
							% tempd = (sum(scores)/size(scores,2))/size(scores,2);
						% end
						% distances(j,i) = tempd;
					% end
					% printing(variables.printing, 'SURF distance', i, pnum);
					% if mod(uint32(pnum/50),i) == 0
						% save([output.resultsdir id output.distance_mat], 'distances', '-v7.3');
					% end
				% end
			% end
			% distances(isnan(distances)) = max(max(distances))+1;
		% end
    else
		if strcmp(variables.feature.name,'SIFT')
			fprintf('Computing SIFT distance\n');
			gnum = size(gallery_features,2);
			pnum = size(probe_features,2);
			istart = 1;
			time10 = 0;
			timeNow = 0;
			if ~isempty(output.distance_mat) && exist([output.resultsdir 'temp' id output.distance_mat],'file')
				distances = load([output.resultsdir 'temp' id output.distance_mat], 'distances', 'i', 'time10');
				istart = distances.i;
				time10 = distances.time10;
				distances = distances.distances;
			else
				distances = zeros(gnum,pnum);
			end
			for i = istart:pnum
				s_t = tic;
				ti = probe_features(:,i);
				for j = 1:gnum
					tj = gallery_features(:,j);
                    tempd = 0; nothing = 1;
                    for k = 1:size(tj,1)
                        [~,scores] = vl_ubcmatch(ti{k},tj{k});
                        if ~isempty(scores)
                            tempd = tempd+(sum(scores)/size(scores,2))/size(scores,2);
                            nothing = 0;
                        end
                    end
                    if nothing == 1, tempd = NaN; end
					distances(j,i) = tempd;
				end
				e_t = toc(s_t);
				timeNow = e_t;
				if i == 10
					time10 = e_t;
				end
				printing(variables.printing, sprintf('SIFT distance (10th: %4.0fs Last: %4.0fs)', time10, timeNow), i, pnum);
				if mod(i,50) == 0 && ~isempty(output.distance_mat)
					save([output.resultsdir 'temp' id output.distance_mat], 'distances', 'i', 'time10', '-v7.3');
				end
			end
			distances(isnan(distances)) = max(max(distances))+1;
		elseif strcmp(variables.feature.name,'SURF')
			fprintf('Computing SURF distance\n');
			gnum = size(gallery_features,2);
			pnum = size(probe_features,2);
			istart = 1;
			time10 = 0;
			timeNow = 0;
			if ~isempty(output.distance_mat) && exist([output.resultsdir 'temp' id output.distance_mat],'file')
				distances = load([output.resultsdir 'temp' id output.distance_mat], 'distances', 'i', 'time10');
				istart = distances.i;
				time10 = distances.time10;
				distances = distances.distances;
			else
				distances = zeros(gnum,pnum);
			end
			for i = istart:pnum
				s_t = tic;
				ti = probe_features{i};
				for j = 1:gnum
					tj = gallery_features{j};
					[~,scores] = vl_ubcmatch(ti,tj);
					if length(scores) == 0
						tempd = NaN;
					else
						tempd = (sum(scores)/size(scores,2))/size(scores,2);
					end
					distances(j,i) = tempd;
				end
				e_t = toc(s_t);
				timeNow = e_t;
				if i == 10
					time10 = e_t;
				end
				printing(variables.printing, sprintf('SURF distance (10th: %4.0fs Last: %4.0fs)', time10, timeNow), i, pnum);
				if mod(i,50) == 0 && ~isempty(output.distance_mat)
					save([output.resultsdir 'temp' id output.distance_mat], 'distances', 'i', 'time10', '-v7.3');
				end
			end
			distances(isnan(distances)) = max(max(distances))+1;
		elseif strcmp(variables.distance,'libsvm')
			gids = variables.gids;
			pids = variables.pids;
			pnum = size(pids,1);
			gnum = size(gids,1);
			num_per_class = histc(gids,unique(gids));
			
			fprintf('1. Train SVM Classifier with %d images\n',gnum);
			t2 = tic;
			% load or create the SVM model
			if ~isempty(output.distance_mat) && exist([output.resultsdir id 'model.mat'],'file')
				model = load([output.resultsdir id 'model.mat'], 'model');
				model = model.model;
			else
				model = svmtrain(gids, gallery_features', '-q -s 0 -t 0');
				if ~isempty(output.distance_mat)
					if(~exist(output.resultsdir, 'dir'))
					   mkdir(output.resultsdir);
					end
					save([output.resultsdir id 'model.mat'], 'model', '-v7.3');
				end
			end
			clearvars gallery_features;
			s=whos('model');
			fprintf('Produced SVM model [%6.2f MB]\n', s.bytes/1000000);
			toc(t2)
			
			fprintf('2. Predict Classification with SVM\n');
			t2 = tic;
			distances = zeros(pnum,1);
			for k=1:pnum
				[predict_label, ~, ~] = svmpredict(pids(k), probe_features(:,k)', model);
				if predict_label == pids(k)
					distances(k) = 1;
				else
					distances(k) = 0;
				end
				printing(variables.printing, sprintf('predict label [%07d] actual [%07d]', predict_label, pids(k)), k, pnum);
			end
			clearvars probe_features model;
			fprintf('\n');
			toc(t2)
			
			% fprintf('2. Predict Classification with SVM\n');
			% t2 = tic;
            % predict_label=zeros(pnum,1);
            % dec_values=zeros(pnum,model.nr_class*(model.nr_class-1)/2);
			% for i=1:pnum
				% [predict_label(i,1), accuracy, dec_values(i,:)] = svmpredict(pids(i), probe_features(:,i)', model);
				% printing(variables.printing, 'predict label', i, pnum);
			% end
			% clearvars probe_features;
			% s=whos('dec_values');
			% fprintf('Produced decision values matrix [%d probes x %d classes] [%6.2f MB]\n', size(dec_values,1), size(dec_values,2), s.bytes/1000000);
			% toc(t2)
			
			% fprintf('3. build votes array\n');
			% t2 = tic;
			% distances = zeros(pnum,model.nr_class);
			% %wvotes = distances;
			% index = 1;
			% for k=1:pnum
			    % index=1;
			    % for i=1:model.nr_class
			        % for j=i+1:model.nr_class
			            % if dec_values(k,index) > 0
			                % distances(k,i)=distances(k,i)+1;
			                % %wvotes(k,i)=wvotes(k,i)+abs(dec_values(k,index));
			            % else
			                % distances(k,j)=distances(k,j)+1;
			                % %wvotes(k,j)=wvotes(k,j)+abs(dec_values(k,index));
			            % end
			            % index=index+1;
			        % end 
                % end
                % printing(variables.printing, 'votes', k, pnum);
			% end
			% clearvars dec_values;
			% fprintf('\n');
			% distances1=abs(distances-max(max(distances))+1)';
			% distances = zeros(gnum,pnum);
			% index=1;
			% for i=1:model.nr_class
				% for j=1:num_per_class(i)
					% distances(index,:) = distances1(i,:);
					% index=index+1;
				% end
			% end
			% clearvars model;
			% toc(t2)
			
		else
			distances = slmetric_pw(gallery_features, probe_features, variables.distance);
		end
		% save distances matrix
        if ~isempty(output.distance_mat)
            if(~exist(output.resultsdir, 'dir'))
               mkdir(output.resultsdir);
            end
            save([output.resultsdir id output.distance_mat], 'distances', '-v7.3');
        end
	end
end