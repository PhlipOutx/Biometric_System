function h = plot_det (pmiss, pfa, eer, name, color, linewidth, flag)
%
%  Plot_DET plots detection performance tradeoff on a DET plot
%  and returns the handle for plotting.
%
%  pmiss and pfa are the vectors of miss and corresponding false
%  alarm probabilities to be plotted.
%
%  The usage of Plot_DET is analogous to the standard matlab
%  plot function.
%
% See DET_usage for an example of how to use Plot_DET.
%
% linewidth : controls the line thickness. The default thickness
% is 0.5. A value between 2 and 5 will give a nice thick line.
%

	len = max(size(pmiss));
	assert(len == max(size(pfa)), 'vector size of pmiss and pfa not equal in call to plot_det');


    %%% create figure if desired %%%
	DET_limits = [0.05 50 ];
	pticks = [ 0.1 0.2 0.5 1 2 5 10 20 40 ];
	ticklabels = [ '0.1'; '0.2'; '0.5'; ' 1 '; ' 2 '; ' 5 '; '10 '; '20 '; '40 ' ];

	h = figure(flag);
	hold on;

	axis('square');

	set(gca, 'xlim', ppndf(DET_limits / 100));
	set(gca, 'xtick', ppndf(pticks / 100));
	set(gca, 'XTickLabel', ticklabels);
	set(gca, 'xgrid', 'on');

	set(gca, 'ylim', ppndf(DET_limits / 100));
	set(gca, 'ytick', ppndf(pticks / 100));
	set(gca, 'YTickLabel', ticklabels);
	set(gca, 'ygrid', 'on');

	set(gca, 'LineWidth', 2);
	set(gca, 'FontSize', 14);

	pfa2 = pfa;
	pmiss2 = pmiss;
	pfa2(len) = pfa(len - 1);
	pmiss2(len) = 90;

	p = plot(ppndf(pfa2), ppndf(pmiss2), color, 'LineWidth', linewidth);
	scatter(ppndf(eer), ppndf(eer), color, 'LineWidth', linewidth, 'Marker', 'd')
	xlabel('False Accept Rate (%)', 'FontSize', 18);
	ylabel('False Reject Rate (%)', 'FontSize', 18);
	legend(p, sprintf('%s - EER: %4.2f%%',name,eer*100), 'FontSize', 14, 'Location', 'Best');
	grid on;
	
end

function norm_dev = ppndf(cum_prob)
    norm_dev = sqrt(2) * erfinv(2 * cum_prob - 1);
end