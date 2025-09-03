% Plot metrics
metric_var_d = Vary_d_n_l;
metric_var_freq = Vary_freq;

plot_stiff_metrics(metric_var_d,metric_var_freq)
% print('Effect-of-params-AD-stiff','-dpdf','-r0')

plot_perm_metrics(metric_var_d,metric_var_freq)
% print('Effect-of-params-AL-perm','-dpdf','-r0')   

plot_freq_metrics_extra(metric_var_freq) % not in paper
% print('Effect-of-freq&l-Extra','-dpdf','-r0')

plot_stiff_metrics_AS_appendix(metric_var_d,metric_var_freq)
% print('Effect-of-params-AS-stiff-appendix','-dpdf','-r0')

plot_stiff_metrics_AD_appendix(metric_var_d,metric_var_freq)
% print('Effect-of-params-AD-stiff-appendix','-dpdf','-r0') 
