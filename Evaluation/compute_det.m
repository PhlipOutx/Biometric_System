function [pmiss, pfa] = compute_det(true_scores, false_scores, resolution)
% INPUTS:
% true_scores           - a vector of true scores (distances)
% false_scores          - a vector of false scores (distances)
%
% OUTPUTS:
% pmiss                 - a vector of size 1 x length(true_scores;false_scores)+1 
%                         containing the false reject probability evaluated at 
%                         the corresponding treshold values
% pfa                   - a vector of size 1 x length(true_scores;false_scores)+1 
%                         containing the false accept probability evaluated at 
%                         the corresponding treshold values
%
% NOTES / COMMENTS
%
% The function was written with Matlab ver. 7.12.0.635 (R2011a).
%
% adapted from the DET-Curve Plotting software for use with MATLAB
% National Institute of Standards and Technology, Information Technology Laboratory
%

	%%% Initilize variables %%%
	num_ts = numel(true_scores);
	num_fs = numel(false_scores);
	num_total = num_ts + num_fs;

	%%% Initilize return variables %%%
	pmiss = zeros(num_ts+num_fs+1, 1);
	pfa = zeros(num_ts+num_fs+1, 1);

	scores(1:num_fs,1) = false_scores;
	scores(1:num_fs,2) = 0;
	scores(num_fs+1:num_total,1) = true_scores;
	scores(num_fs+1:num_total,2) = 1;
	scores = DETsort(scores,1:size(scores,2));

	sumtrue = cumsum(scores(:,2),1);
	sumfalse = num_fs - ([1:num_total]'-sumtrue);

	pmiss(1) = 0;
	pfa(1) = 1.0;
	pmiss(2:num_total+1) = sumtrue ./ num_ts;
	pfa(2:num_total+1) = sumfalse ./ num_fs;
    
    delta = round(length(pmiss)/resolution);
    pmiss = pmiss(1:delta:length(pmiss));
    pfa = pfa(1:delta:length(pfa));

end


function [y,ndx] = DETsort(x,col)
% DETsort Sort rows, the first in ascending, the remaining in decending
% thereby postponing the false alarms on like scores.
% based on SORTROWS

	ndx = (1:size(x,1))';

	% sort 2nd column ascending
	[v,ind] = sort(x(ndx,2));
	ndx = ndx(ind);

	% reverse to decending order
	ndx(1:size(x,1)) = ndx(size(x,1):-1:1);

	% now sort first column ascending
	[v,ind] = sort(x(ndx,1));
	ndx = ndx(ind);
	y = x(ndx,:);

end
