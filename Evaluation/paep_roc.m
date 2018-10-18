function score = paep_roc(true_scores, false_scores, resolution)

	% VR at 0.1% FAR
	false_scores = sort(false_scores);
    threshold = false_scores(round(numel(false_scores)*0.001));
    score = sum(true_scores < threshold)/numel(true_scores);

end