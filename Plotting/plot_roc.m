function h = plot_roc (ver_rate, miss_rate, eer, name, color, linewidth, flag)
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

    %%% create figure if desired %%%
	xlimits = [0.005 100];
	xticks = [ 0.01 0.1 1.0 10 100 ];
	xticklabels = [ '0.01'; '0.1 '; '1.0 '; ' 10 '; '100 ' ];
	ylimits = [0.005 1];
	yticks = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 ];
	yticklabels = [ ' 10'; ' 20'; ' 30'; ' 40'; ' 50'; ' 60'; ' 70'; ' 80'; ' 90'; '100' ];

	h = figure(flag);
	p = semilogx(miss_rate*100, ver_rate, color, 'LineWidth', linewidth);
	%hold on;
	
	set(gca, 'xlim', xlimits);
	set(gca, 'xtick', xticks);
	set(gca, 'XTickLabel', xticklabels);
	set(gca, 'xgrid', 'on');

	set(gca, 'ylim', ylimits);
	set(gca, 'ytick', yticks);
	set(gca, 'YTickLabel', yticklabels);
	set(gca, 'ygrid', 'on');

	set(gca, 'LineWidth', 2);
	set(gca, 'FontSize', 14);

	axis([0.001 100.001 0.001 1.001]);

	xlabel('False Accept Rate (%)', 'FontSize', 18);
	ylabel('Verification Rate (%)', 'FontSize', 18);
	legend(p, sprintf('%s - EER: %4.2f%%',name,eer*100), 'FontSize', 14, 'Location', 'Best');
	grid on;
	
end
