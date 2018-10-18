function results = diss_aging2_exp(input)
% input(struct):	Matlab structure containing fields that hold variables
%                   used in a biometric experiment. See Section 2.0 for
%					list of required fields.
%					
% input(string):	Filepath to xml file that will be parsed into a struct.
%
% results:			See Section 3.0 for list of outputs.

	
	if isstruct(input)
		experiment = input;
	else
		% parse xml experiment file
		experiment = xml_read(input);
	end
    
	fprintf('###Beginning %s\n', experiment.id);
	
	
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG.srt', 'gfLBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year0.srt', 'y0LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year1.srt', 'y1LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year2.srt', 'y2LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year3.srt', 'y3LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year4.srt', 'y4LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year5.srt', 'y5LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year6.srt', 'y6LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year7.srt', 'y7LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year8.srt', 'y8LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year9.srt', 'y9LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year10.srt', 'y10LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year11.srt', 'y11LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year12.srt', 'y12LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year13.srt', 'y13LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year14.srt', 'y14LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year15.srt', 'y15LBP');
	
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG.srt', 'gfHOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year0.srt', 'y0HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year1.srt', 'y1HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year2.srt', 'y2HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year3.srt', 'y3HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year4.srt', 'y4HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year5.srt', 'y5HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year6.srt', 'y6HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year7.srt', 'y7HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year8.srt', 'y8HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year9.srt', 'y9HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year10.srt', 'y10HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year11.srt', 'y11HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year12.srt', 'y12HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year13.srt', 'y13HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year14.srt', 'y14HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year15.srt', 'y15HOG');
	
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG.srt', 'gfLPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year0.srt', 'y0LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year1.srt', 'y1LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year2.srt', 'y2LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year3.srt', 'y3LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year4.srt', 'y4LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year5.srt', 'y5LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year6.srt', 'y6LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year7.srt', 'y7LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year8.srt', 'y8LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year9.srt', 'y9LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year10.srt', 'y10LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year11.srt', 'y11LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year12.srt', 'y12LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year13.srt', 'y13LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year14.srt', 'y14LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year15.srt', 'y15LPQ');
	
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG.srt', 'gfSIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year0.srt', 'y0SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year1.srt', 'y1SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year2.srt', 'y2SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year3.srt', 'y3SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year4.srt', 'y4SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year5.srt', 'y5SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year6.srt', 'y6SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year7.srt', 'y7SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year8.srt', 'y8SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year9.srt', 'y9SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year10.srt', 'y10SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year11.srt', 'y11SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year12.srt', 'y12SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year13.srt', 'y13SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year14.srt', 'y14SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year15.srt', 'y15SIFT');
	
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG.srt', 'gfGABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year0.srt', 'y0GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year1.srt', 'y1GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year2.srt', 'y2GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year3.srt', 'y3GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year4.srt', 'y4GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year5.srt', 'y5GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year6.srt', 'y6GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year7.srt', 'y7GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year8.srt', 'y8GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year9.srt', 'y9GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year10.srt', 'y10GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year11.srt', 'y11GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year12.srt', 'y12GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year13.srt', 'y13GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year14.srt', 'y14GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year15.srt', 'y15GABOR');
	
	
	
	
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG0.srt', 'g0LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG1.srt', 'g1LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG2.srt', 'g2LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG3.srt', 'g3LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG4.srt', 'g4LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG5.srt', 'g5LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG6.srt', 'g6LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG7.srt', 'g7LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG8.srt', 'g8LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG9.srt', 'g9LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG10.srt', 'g10LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG11.srt', 'g11LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG12.srt', 'g12LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG13.srt', 'g13LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG14.srt', 'g14LBP');
	extractTheseFeatures(experiment, 'LBP', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG15.srt', 'g15LBP');
	
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG0.srt', 'g0HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG1.srt', 'g1HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG2.srt', 'g2HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG3.srt', 'g3HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG4.srt', 'g4HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG5.srt', 'g5HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG6.srt', 'g6HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG7.srt', 'g7HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG8.srt', 'g8HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG9.srt', 'g9HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG10.srt', 'g10HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG11.srt', 'g11HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG12.srt', 'g12HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG13.srt', 'g13HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG14.srt', 'g14HOG');
	extractTheseFeatures(experiment, 'HOG', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG15.srt', 'g15HOG');
	
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG0.srt', 'g0LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG1.srt', 'g1LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG2.srt', 'g2LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG3.srt', 'g3LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG4.srt', 'g4LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG5.srt', 'g5LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG6.srt', 'g6LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG7.srt', 'g7LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG8.srt', 'g8LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG9.srt', 'g9LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG10.srt', 'g10LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG11.srt', 'g11LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG12.srt', 'g12LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG13.srt', 'g13LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG14.srt', 'g14LPQ');
	extractTheseFeatures(experiment, 'LPQ', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG15.srt', 'g15LPQ');
	
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG0.srt', 'g0SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG1.srt', 'g1SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG2.srt', 'g2SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG3.srt', 'g3SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG4.srt', 'g4SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG5.srt', 'g5SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG6.srt', 'g6SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG7.srt', 'g7SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG8.srt', 'g8SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG9.srt', 'g9SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG10.srt', 'g10SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG11.srt', 'g11SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG12.srt', 'g12SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG13.srt', 'g13SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG14.srt', 'g14SIFT');
	extractTheseFeatures(experiment, 'SIFT', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG15.srt', 'g15SIFT');
	
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG0.srt', 'g0GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG1.srt', 'g1GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG2.srt', 'g2GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG3.srt', 'g3GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG4.srt', 'g4GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG5.srt', 'g5GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG6.srt', 'g6GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG7.srt', 'g7GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG8.srt', 'g8GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG9.srt', 'g9GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG10.srt', 'g10GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG11.srt', 'g11GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG12.srt', 'g12GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG13.srt', 'g13GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG14.srt', 'g14GABOR');
	extractTheseFeatures(experiment, 'GABOR', '\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG15.srt', 'g15GABOR');
	
	
	
	% LBP
	% LBPrank1 = zeros(1,16);
	% LBPeer = zeros(1,16);
	% for j = 0:15
		% features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'LBP.mat'], 'features'); gallery_features = features.features; clearvars features;
		% features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'LBP.mat'], 'features'); probe_features = features.features; clearvars features;
		
		% experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		% experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		% experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
		% experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
		% distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		% clearvars gallery_features probe_features;
		% CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% LBPrank1(j+1) = CMC_rec_rates(1)*100;
		% [true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% clearvars distances;
		% [~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		% LBPeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		% printing(experiment.variables.printing, 'LBP experiments', j+1, 16);
	% end
    % fprintf('\nLBP Results: Rank-1\n');
	% for j = 1:16 
		% fprintf('%f,',LBPrank1(1,j));
	% end
	% fprintf('\n');
	% fprintf('\nLBP Results: EER\n');
	% for j = 1:16 
		% fprintf('%f,',LBPeer(1,j));
	% end
	% fprintf('\n\n');
	
	
	
	% HOG
	% HOGrank1 = zeros(1,16);
	% HOGeer = zeros(1,16);
	% for j = 0:15
		% features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'HOG.mat'], 'features'); gallery_features = features.features; clearvars features;
		% features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'HOG.mat'], 'features'); probe_features = features.features; clearvars features;
		
		% experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		% experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		% experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
		% experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
		% distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		% clearvars gallery_features probe_features;
		% CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% HOGrank1(j+1) = CMC_rec_rates(1)*100;
		% [true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
		% clearvars distances;
		% [~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		% HOGeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		% printing(experiment.variables.printing, 'HOG experiments', j+1, 16);
	% end
    % fprintf('\nHOG Results: Rank-1\n');
	% for j = 1:16 
		% fprintf('%f,',HOGrank1(1,j));
	% end
	% fprintf('\n');
	% fprintf('\nHOG Results: EER\n');
	% for j = 1:16 
		% fprintf('%f,',HOGeer(1,j));
	% end
	% fprintf('\n\n');
	
	
	% LPQ
	% LPQrank1 = zeros(1,16);
	% LPQeer = zeros(1,16);
	% for j = 0:15
		% features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'LPQ.mat'], 'features'); gallery_features = features.features; clearvars features;
		% features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'LPQ.mat'], 'features'); probe_features = features.features; clearvars features;
		
		% experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		% experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		% experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
		% experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
		% distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		% clearvars gallery_features probe_features;
		% CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% LPQrank1(j+1) = CMC_rec_rates(1)*100;
		% [true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
		% clearvars distances;
		% [~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		% LPQeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		% printing(experiment.variables.printing, 'LPQ experiments', j+1, 16);
	% end
    % fprintf('\nLPQ Results: Rank-1\n');
	% for j = 1:16 
		% fprintf('%f,',LPQrank1(1,j));
	% end
	% fprintf('\n');
	% fprintf('\nLPQ Results: EER\n');
	% for j = 1:16 
		% fprintf('%f,',LPQeer(1,j));
	% end
	% fprintf('\n\n');
	
	% SIFT
	% experiment.variables.feature.name = 'SIFT';
	% SIFTrank1 = zeros(1,16);
	% SIFTeer = zeros(1,16);
	% for j = 0:15
		% features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'SIFT.mat'], 'features'); gallery_features = features.features; clearvars features;
		% features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'SIFT.mat'], 'features'); probe_features = features.features; clearvars features;
		
		% experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		% experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		% experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
		% experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
		% distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		% clearvars gallery_features probe_features;
		% CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		% SIFTrank1(j+1) = CMC_rec_rates(1)*100;
		% [true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
		% clearvars distances;
		% [~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		% SIFTeer(j+1) = ROC_rates_and_threshs.EER_er*100;
		% printing(experiment.variables.printing, 'SIFT experiments', j+1, 16);
	% end
    % fprintf('\nSIFT Results: Rank-1\n');
	% for j = 1:16 
		% fprintf('%f,',SIFTrank1(1,j));
	% end
	% fprintf('\n');
	% fprintf('\nSIFT Results: EER\n');
	% for j = 1:16 
		% fprintf('%f,',SIFTeer(1,j));
	% end
	% fprintf('\n\n');
	
	% GABOR
	experiment.variables.feature.name = 'GABOR';
	GABORrank1 = zeros(1,16);
	GABOReer = zeros(1,16);
	for j = 0:15
		features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'GABOR.mat'], 'features'); gallery_features = features.features; clearvars features;
		features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'GABOR.mat'], 'features'); probe_features = features.features; clearvars features;
		
		experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
		experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
		experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
		experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
		distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
		clearvars gallery_features probe_features;
		CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
		GABORrank1(j+1) = CMC_rec_rates(1)*100;
		[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
		clearvars distances;
		[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
		GABOReer(j+1) = ROC_rates_and_threshs.EER_er*100;
		printing(experiment.variables.printing, 'GABOR experiments', j+1, 16);
	end
    fprintf('\nGABOR Results: Rank-1\n');
	for j = 1:16 
		fprintf('%f,',GABORrank1(1,j));
	end
	fprintf('\n');
	fprintf('\nGABOR Results: EER\n');
	for j = 1:16 
		fprintf('%f,',GABOReer(1,j));
	end
	fprintf('\n\n');
 
end

function extractTheseFeatures(experiment, feat1, list1, name1)

	% feature extraction
	if(~exist([experiment.output.resultsdir experiment.id name1 '.mat'], 'file'))
		experiment.variables.feature.name = feat1;
		
		% extract gallery features
		t1 = tic;
		features = extract_features(experiment.id, experiment.variables, experiment.input, experiment.output, list1, experiment.input.datadir, '', '');
		s=whos('features');
		fprintf('Produced feature matrix [%d features x %d images] [%6.2f MB]\n', size(features,1), size(features,2), s.bytes/1000000);
		toc(t1)
		fprintf('\n');
		
		save([experiment.output.resultsdir experiment.id name1 '.mat'], 'features', '-v7.3'); clearvars features;
	end

end
