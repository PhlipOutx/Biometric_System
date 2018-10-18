function h = plot_cmc (rec_rates, name, color, linewidth, flag)
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

	len = max(size(rec_rates));

    %%% create figure if desired %%%
	xlimits = [0.05 101];
	xticks = [ 1 10 100 1000 ];
	xticklabels = [ ' 1  '; ' 10 '; '100 '; '1000' ];
	ylimits = [0.05 101];
	yticks = [ 10 20 30 40 50 60 70 80 90 100 ];
	yticklabels = [ ' 10'; ' 20'; ' 30'; ' 40'; ' 50'; ' 60'; ' 70'; ' 80'; ' 90'; '100' ];

	h = figure(flag);
	p = semilogx(1:len, rec_rates, color, 'LineWidth', linewidth);
	%hold on;
	
	set(gca, 'xlim', xlimits);
	set(gca, 'xtick', xticks);
	set(gca, 'XTickLabel', xticklabels);
	set(gca, 'xgrid', 'on');

	set(gca, 'ylim', ylimits / 100);
	set(gca, 'ytick', yticks / 100);
	set(gca, 'YTickLabel', yticklabels);
	set(gca, 'ygrid', 'on');

	set(gca, 'LineWidth', 2);
	set(gca, 'FontSize', 14);

	axis([1 1001 0.001 1.001]);

	xlabel('Rank', 'FontSize', 18);
	ylabel('Recognition Rate (%)', 'FontSize', 18);
	legend(p, sprintf('%s - Rank1: %4.2f%%',name,rec_rates(1)*100), 'FontSize', 14, 'Location', 'Best');
	grid on;
	
end
