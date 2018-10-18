function results = diss_aging_exp(input)
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
	
	
	flength = [59 12 256];
	if strcmp(experiment.variables.patches.type,'miller');
		loc = [1,11;12,22;23,26;27,30;31,48;49,63;64,84];
	elseif strcmp(experiment.variables.patches.type,'miller2');
		loc = [1,7;8,12;13,16;17,20;21,28;29,34;35,74];
	end
	
	
	
	
	
	% LBP
	LBPrank1 = zeros(8,16);
	LBPeer = zeros(8,16);
	flength = 59;
	for i = 1:8
		for j = 0:15
			features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'LBP.mat'], 'features'); gfLBP = features.features; clearvars features;
			features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'LBP.mat'], 'features'); pfLBP = features.features; clearvars features;
			if i==8
				gallery_features = gfLBP;
				probe_features =  pfLBP;
			else
				gallery_features =  gfLBP(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
				probe_features =  pfLBP(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
			end
			clearvars pfLBP gfLBP;
			
			experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
			experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
			experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
			experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
			distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
			clearvars gallery_features probe_features;
			CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
			LBPrank1(i,j+1) = CMC_rec_rates(1)*100;
			[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
			clearvars distances;
			[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
			LBPeer(i,j+1) = ROC_rates_and_threshs.EER_er*100;
			printing(experiment.variables.printing, 'LBP experiments', (i-1)*16+j+1, 128);
		end
	end
	clearvars gfLBP;
    fprintf('\nLBP Results: Rank-1\n');
	for i = 1:8
		for j = 1:16 
			fprintf('%f,',LBPrank1(i,j));
		end
		fprintf('\n');
	end
	fprintf('\nLBP Results: EER\n');
	for i = 1:8
		for j = 1:16 
			fprintf('%f,',LBPeer(i,j));
		end
		fprintf('\n');
	end
	fprintf('\n');
	
	
	
	% HOG
	HOGrank1 = zeros(8,16);
	HOGeer = zeros(8,16);
	flength = 12;
	for i = 1:8
		for j = 0:15
			features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'HOG.mat'], 'features'); gfHOG = features.features; clearvars features;
			features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'HOG.mat'], 'features'); pfHOG = features.features; clearvars features;
			if i==8
				gallery_features = gfHOG;
				probe_features =  pfHOG;
			else
				gallery_features =  gfHOG(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
				probe_features =  pfHOG(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
			end
			clearvars pfHOG gfHOG;
			
			experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
			experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
			experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
			experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
			distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
			clearvars gallery_features probe_features;
			CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
			HOGrank1(i,j+1) = CMC_rec_rates(1)*100;
			[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
			clearvars distances;
			[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
			HOGeer(i,j+1) = ROC_rates_and_threshs.EER_er*100;
			printing(experiment.variables.printing, 'HOG experiments', (i-1)*16+j+1, 128);
		end
	end
	clearvars gfHOG;
    fprintf('\nHOG Results: Rank-1\n');
	for i = 1:8
		for j = 1:16 
			fprintf('%f,',HOGrank1(i,j));
		end
		fprintf('\n');
	end
	fprintf('\nHOG Results: EER\n');
	for i = 1:8
		for j = 1:16 
			fprintf('%f,',HOGeer(i,j));
		end
		fprintf('\n');
	end
	fprintf('\n');
	
	
	% LPQ
	LPQrank1 = zeros(8,16);
	LPQeer = zeros(8,16);
	flength = 256;
	for i = 1:8
		for j = 0:15
			features = load([experiment.output.resultsdir experiment.id 'g' num2str(j) 'LPQ.mat'], 'features'); gfLPQ = features.features; clearvars features;
			features = load([experiment.output.resultsdir experiment.id 'y' num2str(j) 'LPQ.mat'], 'features'); pfLPQ = features.features; clearvars features;
			if i==8
				gallery_features = gfLPQ;
				probe_features =  pfLPQ;
			else
				gallery_features =  gfLPQ(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
				probe_features =  pfLPQ(((loc(i,1)-1)*flength)+1:loc(i,2)*flength,:);
			end
			clearvars pfLPQ gfLPQ;
			
			experiment.input.gallery_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_yearG' num2str(j) '.srt'];
			experiment.input.probe_id = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_ids_year' num2str(j) '.srt'];
			experiment.input.gallery = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_yearG' num2str(j) '.srt'];
			experiment.input.probe = ['\\skynet\users\pemille\imagelists\2012_11_Dissertation2\pinellas_lefteye_year' num2str(j) '.srt'];
			distances = compute_distances(gallery_features, probe_features, experiment.id, experiment.variables, experiment.input, experiment.output);
			clearvars gallery_features probe_features;
			CMC_rec_rates = compute_cmc(experiment.id, distances, experiment.variables, experiment.input, experiment.output);
			LPQrank1(i,j+1) = CMC_rec_rates(1)*100;
			[true_scores false_scores] = split_true_false_scores(experiment.id, distances, experiment.variables, experiment.input, experiment.output); 
			clearvars distances;
			[~, ~, ROC_rates_and_threshs] = compute_roc(true_scores, false_scores, experiment.output.roc_resolution);
			LPQeer(i,j+1) = ROC_rates_and_threshs.EER_er*100;
			printing(experiment.variables.printing, 'LPQ experiments', (i-1)*16+j+1, 128);
		end
	end
	clearvars pfLPQ;
    fprintf('\nLPQ Results: Rank-1\n');
	for i = 1:8
		for j = 1:16 
			fprintf('%f,',LPQrank1(i,j));
		end
		fprintf('\n');
	end
	fprintf('\nLPQ Results: EER\n');
	for i = 1:8
		for j = 1:16 
			fprintf('%f,',LPQeer(i,j));
		end
		fprintf('\n');
	end
	fprintf('\n');
 
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
