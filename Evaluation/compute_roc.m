function [ver_rate, miss_rate, rates_and_threshs] = compute_roc(true_scores, false_scores, resolution)
% INPUTS:
% true_scores           - a vector of true scores (distances)
% false_scores          - a vector of false scores (distances)
% resolution            - the number of points in the ROC curve
%
% OUTPUTS:
% ver_rate              - a vector of size 1 x (resolution+1) containing the 
%                         verification rates of the ROC curve evaluated at 
%                         the corresponding treshold values; the values in
%                         the vector stand for 1-FRR; where FRR is the
%                         false rejection rate
% miss_rate             - a vector of size 1 x (resolution+1) containing the miss 
%                         verification rates of the ROC curve evaluated at 
%                         the corresponding treshold values; the values in
%                         the vector stand for 1-FAR; where FAR is the
%                         false acceptance rate
% rates_and_threshs     - a structure with data computed at characteristic
%                         operating points on the ROC curve; the structure 
%                         contains the following fields: 
%                         
%         .EER_er       - equal error rate - the value where FAR=FRR
%         .EER_tr       - the threshold needed to obtain the EER 
%         .EER_frr      - the value of the FRR at the EER 
%         .EER_ver      - the value of the verification rate at the EER                      
%         .EER_far      - the value of the FAR at the EER 
% 
% 
%         .FRR_01FAR_er  - the value of the half total error rate at the 
%                          ROC operating point that ensures that the false 
%                          acceptance rate is 10 times higher than the 
%                          false rejection rate; i.e., FAR=0.1FRR
%         .FRR_01FAR_tr  - the threshold needed to obtain FAR=0.1FRR 
%         .FRR_01FAR_frr - the value of the FRR at FAR=0.1FRR 
%         .FRR_01FAR_ver - the value of the verification rate at FAR=0.1FRR                      
%         .FRR_01FAR_far - the value of the FAR at FAR=0.1FRR 
% 
% 
%         .FRR_10FAR_er  - the value of the half total error rate at the 
%                          ROC operating point that ensures that the false 
%                          rejection rate is 10 times higher than the 
%                          false acceptance rate; i.e., FAR=10FRR
%         .FRR_10FAR_tr  - the threshold needed to obtain FAR=10FRR 
%         .FRR_10FAR_frr - the value of the FRR at FAR=10FRR 
%         .FRR_10FAR_ver - the value of the verification rate at FAR=10FRR                      
%         .FRR_10FAR_far - the value of the FAR at FAR=10FRR
% 
% 
%         .VER_001FAR_er  - the value of the half total error rate at the 
%                           ROC operating point where the FAR equals 0.01% 
%         .VER_001FAR_tr  - the threshold needed to a FAR of 0.01%
%         .VER_001FAR_frr - the value of the FRR a FAR of 0.01% 
%         .VER_001FAR_ver - the value of the verification rate at a FAR of 
%                           0.01%                      
%         .VER_001FAR_far - the value of the FAR at a FAR of 0.01% (note 
%                           that this is the actual value of the FAR, since 
%                           there might not be anough data to obtain an FAR 
%                           of exactly 0.01%)
% 
% 
%         .VER_01FAR_er  - the value of the half total error rate at the 
%                          ROC operating point where the FAR equals 0.1% 
%         .VER_01FAR_tr  - the threshold needed to a FAR of 0.1%
%         .VER_01FAR_frr - the value of the FRR a FAR of 0.1% 
%         .VER_01FAR_ver - the value of the verification rate at a FAR of 
%                          0.1%                      
%         .VER_01FAR_far - the value of the FAR at a FAR of 0.1% (note 
%                          that this is the actual value of the FAR, since 
%                          there might not be anough data to obtain an FAR 
%                          of exactly 0.1%)
% 
% 
%         .VER_1FAR_er  - the value of the half total error rate at the 
%                          ROC operating point where the FAR equals 1% 
%         .VER_1FAR_tr  - the threshold needed to a FAR of 1%
%         .VER_1FAR_frr - the value of the FRR a FAR of 1% 
%         .VER_1FAR_ver - the value of the verification rate at a FAR of 
%                          1%                      
%         .VER_1FAR_far - the value of the FAR at a FAR of 1% (note 
%                          that this is the actual value of the FAR, since 
%                          there might not be anough data to obtain an FAR 
%                          of exactly 1%)
%                         
%
% NOTES / COMMENTS
%
% The function was written with Matlab ver. 7.12.0.635 (R2011a).
%
% adapted from the PhD face recognition toolbox
% Copyright (c) 2011 Vitomir Štruc
%

	%%% Initilize return variables %%%
	ver_rate = [];
	miss_rate = [];
	rates_and_threshs = [];
	
	if isempty(true_scores)
		true_scores = 0;
	end
	if isempty(false_scores)
		false_scores = 0;
	end

	%%% Initilize variables %%%
	num_ts = numel(true_scores);
	num_fs = numel(false_scores);
	min_ts = min(true_scores);
	max_fs = max(false_scores);
	
	if (min_ts == max_fs)
		rates_and_threshs.dprime = 0;
		rates_and_threshs.EER_er  = 1;
		rates_and_threshs.EER_tr  = 0;
		rates_and_threshs.EER_frr = 1;
		rates_and_threshs.EER_ver = 0;
		rates_and_threshs.EER_far = 1;
		rates_and_threshs.FRR_01FAR_er  = 1;
		rates_and_threshs.FRR_01FAR_tr  = 0;
		rates_and_threshs.FRR_01FAR_frr = 1;
		rates_and_threshs.FRR_01FAR_ver = 0;
		rates_and_threshs.FRR_01FAR_far = 1;
		rates_and_threshs.FRR_10FAR_er  = 1;
		rates_and_threshs.FRR_10FAR_tr  = 0;
		rates_and_threshs.FRR_10FAR_frr = 1;
		rates_and_threshs.FRR_10FAR_ver = 0;
		rates_and_threshs.FRR_10FAR_far = 1;
		rates_and_threshs.VER_001FAR_er  = 1;
		rates_and_threshs.VER_001FAR_tr  = 0;
		rates_and_threshs.VER_001FAR_frr = 1;
		rates_and_threshs.VER_001FAR_ver = 0;
		rates_and_threshs.VER_001FAR_far = 1;
		rates_and_threshs.VER_01FAR_er  = 1;
		rates_and_threshs.VER_01FAR_tr  = 0;
		rates_and_threshs.VER_01FAR_frr = 1;
		rates_and_threshs.VER_01FAR_ver = 0;
		rates_and_threshs.VER_01FAR_far = 1;
		rates_and_threshs.VER_1FAR_er  = 1;
		rates_and_threshs.VER_1FAR_tr  = 0;
		rates_and_threshs.VER_1FAR_frr = 1;
		rates_and_threshs.VER_1FAR_ver = 0;
		rates_and_threshs.VER_1FAR_far = 1;
		ver_rate = 0;
		miss_rate = 1;
		return
	end

	% computing false reject errors
	delta = (max_fs-min_ts)/resolution;
	counter = 1;
	templen = length(min_ts:delta:max_fs);
	fre = zeros(1,resolution);
	for threshold = min_ts:delta:max_fs
		errors = sum(true_scores < threshold);
		fre(1,counter) = 1-(errors/num_ts);
		counter = counter+1;
		printing(2, 'computing false reject errors', counter, templen);
	end
	fprintf('\n');

	% computing false accept errors
	counter = 1;
	fae = zeros(1,resolution);
	for threshold = min_ts:delta:max_fs
		errors = sum(false_scores < threshold);
		fae(1,counter) = (errors/num_fs);
		counter = counter+1;
		printing(2, 'computing false accept errors', counter, templen);
	end
	fprintf('\n');


	% find thresholds for EER, FRR = 0.1FAR, FRR = 10FAR, VR @ 0.01%FAR, VR @ 0.1%FAR, VR @ 1%FAR
	maxi1 = Inf;	maxi2 = Inf;	maxi3 = Inf;	maxi4 = Inf;	maxi5 = Inf;	maxi6 = Inf;
	for i = 1:resolution+1
		% EER
		if abs(fae(i)-fre(i)) < maxi1
		   index1 = i;
		   maxi1 = abs(fae(i)-fre(i));
		end
		
		% FRR = 0.1FAR
		if abs(0.1*fae(i)-fre(i)) < maxi2
		   index2 = i;
		   maxi2 = abs(0.1*fae(i)-fre(i));
		end
		
		% FRR = 10FAR
		if abs(10*fae(i)-fre(i)) < maxi3
		   index3 = i;
		   maxi3 = abs(10*fae(i)-fre(i));
		end
		
		% VR @ 0.01%FAR
		if abs(fae(i)-0.01/100) < maxi4
		   index4 = i;
		   maxi4 = abs(fae(i)-0.01/100);
		end
		
		% VR @ 0.1%FAR
		if abs(fae(i)-0.1/100) < maxi5
		   index5 = i;
		   maxi5 = abs(fae(i)-0.1/100);
		end
		
		% VR @ 1%FAR
		if abs(fae(i)-1/100) < maxi6
		   index6 = i;
		   maxi6 = abs(fae(i)-1/100);
		end
		
		printing(2, 'compute performance', i, resolution+1);
	end
	fprintf('\n');
	
	% D'
	rates_and_threshs.dprime = abs(mean(true_scores) - mean(false_scores))/sqrt((std(true_scores)^2+std(false_scores)^2)/2);

	% EER
	C=fae+fre;
	rates_and_threshs.EER_er  = C(index1)/2;
	rates_and_threshs.EER_tr  = min_ts+(index1-1)*delta;
	rates_and_threshs.EER_frr = sum(true_scores>(min_ts+(index1-1)*delta))/num_ts;
	rates_and_threshs.EER_ver = 1-rates_and_threshs.EER_frr;
	rates_and_threshs.EER_far = sum(false_scores<(min_ts+(index1-1)*delta))/num_fs;

	% FRR = 0.1FAR
	rates_and_threshs.FRR_01FAR_er  = C(index2)/2;
	rates_and_threshs.FRR_01FAR_tr  = min_ts+(index2-1)*delta;
	rates_and_threshs.FRR_01FAR_frr = sum(true_scores>(min_ts+(index2-1)*delta))/num_ts;
	rates_and_threshs.FRR_01FAR_ver = 1-rates_and_threshs.FRR_01FAR_frr;
	rates_and_threshs.FRR_01FAR_far = sum(false_scores<(min_ts+(index2-1)*delta))/num_fs;

	% FRR = 10FAR
	rates_and_threshs.FRR_10FAR_er  = C(index3)/2;
	rates_and_threshs.FRR_10FAR_tr  = min_ts+(index3-1)*delta;
	rates_and_threshs.FRR_10FAR_frr = sum(true_scores>(min_ts+(index3-1)*delta))/num_ts;
	rates_and_threshs.FRR_10FAR_ver = 1-rates_and_threshs.FRR_10FAR_frr;
	rates_and_threshs.FRR_10FAR_far = sum(false_scores<(min_ts+(index3-1)*delta))/num_fs;

	% VR @ 0.01%FAR
	rates_and_threshs.VER_001FAR_er  = C(index4)/2;
	rates_and_threshs.VER_001FAR_tr  = min_ts+(index4-1)*delta;
	rates_and_threshs.VER_001FAR_frr = sum(true_scores>(min_ts+(index4-1)*delta))/num_ts;
	rates_and_threshs.VER_001FAR_ver = 1-rates_and_threshs.VER_001FAR_frr;
	rates_and_threshs.VER_001FAR_far = sum(false_scores<(min_ts+(index4-1)*delta))/num_fs;

	% VR @ 0.1%FAR
	rates_and_threshs.VER_01FAR_er  = C(index5)/2;
	rates_and_threshs.VER_01FAR_tr  = min_ts+(index5-1)*delta;
	rates_and_threshs.VER_01FAR_frr = sum(true_scores>(min_ts+(index5-1)*delta))/num_ts;
	rates_and_threshs.VER_01FAR_ver = 1-rates_and_threshs.VER_01FAR_frr;
	rates_and_threshs.VER_01FAR_far = sum(false_scores<(min_ts+(index5-1)*delta))/num_fs;

	% VR @ 1%FAR
	rates_and_threshs.VER_1FAR_er  = C(index6)/2;
	rates_and_threshs.VER_1FAR_tr  = min_ts+(index6-1)*delta;
	rates_and_threshs.VER_1FAR_frr = sum(true_scores>(min_ts+(index6-1)*delta))/num_ts;
	rates_and_threshs.VER_1FAR_ver = 1-rates_and_threshs.VER_1FAR_frr;
	rates_and_threshs.VER_1FAR_far = sum(false_scores<(min_ts+(index6-1)*delta))/num_fs;

	% set outputs
	ver_rate = 1- fre;
	miss_rate = fae;
end