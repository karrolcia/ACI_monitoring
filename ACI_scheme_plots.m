%% ACI Scheme - plots only
%% Create scatter plots of the selected data
% Create scatter plots based on the LWP interval
s = 20.0 ;
% Select only the data points that where retrieved with the microphysical
% retrieval
is_data = find(nansum(re_hm,2) ~= 0);

data_lidar = (int_beta_atten(is_data,:));
data_radar = (int_Z(is_data,:));
data_doppler = (doppler_v(is_data,:));
data_reff = (re_hm(is_data,:));
data_ext_cloud = ext_hm(is_data,:);
data_lwc_cloud = lwc_hm(is_data,:);
data_Nd = N_hm(is_data);
data_tau_cloud = tau_hm(is_data);

%% Define level to compare
aerosols = I_max - 9;
aerosols(aerosols <= 0) = 1;
    
cloud = I_max+2;
aero_level = aerosols(is_data) ; % Attenuated Backscatter
cloud_level = cloud(is_data)   ; % Radar Reflectivity Factor

for i = 1:length(is_data)
data_lidar_level(i) = data_lidar(i,aero_level(i));
data_radar_level(i) = data_radar(i,cloud_level(i));
data_doppler_level(i) = data_doppler(i,cloud_level(i));
data_reff_level(i) = data_reff(i,cloud_level(i));
end

%% Define different LWP intervals
lwp_cloud = lwp(is_data);
if min(lwp) < 30
    min_lwp = 30;
else min_lwp = min(lwp) ;
end
max_lwp = 90;
% max_lwp = max(lwp_cloud);
id_cloud = (lwp_cloud > min_lwp & lwp_cloud < max_lwp );
% data_Nd = data_Nd(id_cloud);
% data_tau_cloud = data_tau_cloud(id_cloud);

%% Define subplot net
sx = 3; % number of rows
sy = 2; %number of columns
n_bins = sx*sy; % must be divided by 2
bin_size =  round((max_lwp - min_lwp)/n_bins);
c_min = min(lwp_cloud);
c_max = max(lwp_cloud);

for i = 1:n_bins
  
 id_lwp = find(lwp_cloud >= min_lwp + bin_size*(i-1) & ...
         lwp_cloud <= min_lwp + bin_size*i);
 lwp_id{i} = id_lwp;
 % Divide data into bins based on the LWP
    for j = 1:length(lwp_id{i})
      corr_lidar(i,j) = (data_lidar(lwp_id{i}(j),aero_level(lwp_id{i}(j))));
      corr_radar(i,j) = (data_radar(lwp_id{i}(j),cloud_level(lwp_id{i}(j))));
      corr_reff(i,j) = (data_reff(lwp_id{i}(j),cloud_level(lwp_id{i}(j))));
%       corr_ext(i,j) = (data_ext_cloud(lwp_id{i}(j),cloud_level(lwp_id{i}(j))));
%       corr_lwc(i,j) = (data_lwc_cloud(lwp_id{i}(j),cloud_level(lwp_id{i}(j))));
    end
     
 % Calculate correlation coefficient for each scatter plot
  data_lidar_corr{i} = double(corr_lidar(i,:));
  data_radar_corr{i} = double(corr_radar(i,:));
  data_reff_corr{i} = double(corr_reff(i,:));
%   data_ext_corr{i} = double(corr_ext(i,:));
%   data_lwc_corr{i} = double(corr_lwc(i,:));
%   data_Nd_corr{i} = double(data_Nd(lwp_id{i})) ;
%   data_tau_corr{i} = double(data_tau_cloud(lwp_id{i})) ;
  data_lwp_corr{i} = double(lwp_cloud(lwp_id{i})) ;
  data_time{i} = double(time(lwp_id{i}));
  
  % Calculate correlations
  [correlations_Z_beta{i}, pvalue_Z_beta{i}] = ...
      corrcoef(log(data_lidar_corr{i}),log(data_radar_corr{i}));
  [correlations_reff_beta{i}, pvalue_reff_beta{i}] = ...
      corrcoef(log(data_lidar_corr{i}),log(data_reff_corr{i}));
%   [correlations_ext_beta{i}, pvalue_ext_beta{i}] = ...
%       corrcoef(data_lidar_corr{i},data_ext_corr{i});
%   [correlations_lwc_beta{i}, pvalue_lwc_beta{i}] = ...
%       corrcoef(data_lidar_corr{i},data_lwc_corr{i});
%   [correlations_Nd_beta{i}, pvalue_Nd_beta{i} ]= ...
%       corrcoef(log(data_lidar_corr{i}),log(data_Nd_corr{i}));
%   [correlations_tau_beta{i}, pvalue_tau_beta{i}] = ...
%       corrcoef(data_lidar_corr{i},data_tau_corr{i});
  
  % Calculate slope and R^2
  [poly_Z_beta{i}, stuct_Z_beta{i}, mu_Z_beta{i}] = ...
      polyfit(data_lidar_corr{i}, data_radar_corr{i},1);
  [poly_Z_aci{i}, stuct_Z_aci{i}, mu_Z_aci{i}] = ...
      polyfit(log(data_lidar_corr{i}), log(data_radar_corr{i}),1);
  [poly_reff_beta{i}, stuct_reff_beta{i}, mu_reff_beta{i}] = ...
      polyfit(data_lidar_corr{i}, data_reff_corr{i},1);
  [poly_reff_aci{i}, stuct_aci_beta{i}] = ...
      polyfit(log(data_lidar_corr{i}), log(data_reff_corr{i}),1);
%   [poly_ext_beta{i}, stuct_ext_beta{i}, mu_ext_beta{i}] = ...
%       polyfit(data_lidar_corr{i}, data_ext_corr{i},1);
%   [poly_lwc_beta{i}, stuct_lwc_beta{i}, mu_lwc_beta{i}] = ...
%       polyfit(data_lidar_corr{i}, data_lwc_corr{i},1);
%   [poly_Nd_beta{i}, stuct_Nd_beta{i}, mu_Nd_beta{i}] = ...
%       polyfit(data_lidar_corr{i}, data_Nd_corr{i}',1);
%   [poly_Nd_aci{i}, stuct_Nd_aci{i}] = ...
%       polyfit(log(data_lidar_corr{i}), log(data_Nd_corr{i}'),1);
%   [poly_tau_beta{i}, stuct_tau_beta{i}, mu_tau_beta{i}] = ...
%       polyfit(data_lidar_corr{i}, data_tau_corr{i}',1);
    
%   [yfit_Z_beta{i}, delta_Z_beta{i}] = ...
%       polyval(poly_Z_beta{i},data_lidar_corr{i},stuct_Z_beta{i});
%   [yfit_reff_beta{i}, delta_reff_beta{i}] = ...
%       polyval(poly_reff_beta{i},data_lidar_corr{i},stuct_reff_beta{i});
%   [yfit_Nd_beta{i}, delta_Nd_beta{i}] = ...
%       polyval(poly_Nd_beta{i},data_lidar_corr{i},stuct_Nd_beta{i});
%     
%   yresid_Z_beta {i} =  data_radar_corr{i} - yfit_Z_beta{i};
%   yresid_reff_beta {i} =  data_reff_corr{i} - yfit_reff_beta{i};
%   yresid_Nd_beta {i} =  data_Nd_corr{i} - yfit_Nd_beta{i}';
%   SSresid_Z_beta{i} = sum(yresid_Z_beta{i}.^2);
%   SSresid_reff_beta{i} = sum(yresid_reff_beta{i}.^2);
%   SSresid_Nd_beta{i} = sum(yresid_Nd_beta{i}.^2);
%   SStotal_Z_beta{i} = (length(data_radar_corr{i})-1) * var(data_radar_corr{i});
%   SStotal_reff_beta{i} = (length(data_reff_corr{i})-1) * var(data_reff_corr{i});
%   SStotal_Nd_beta{i} = (length(data_Nd_corr{i})-1) * var(data_Nd_corr{i});
%   rsq_Z_beta{i} = 1 - SSresid_Z_beta{i}./SStotal_Z_beta{i};
%   rsq_reff_beta{i} = 1 - SSresid_reff_beta{i}./SStotal_reff_beta{i};
%   rsq_Nd_beta{i} = 1 - SSresid_Nd_beta{i}./SStotal_Nd_beta{i};

clear id_lwp corr_*
end

%% Plot time series of measurements
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) ' data time series'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'units','centimeters','Position',[2 50 8.3 15]);
% subplot(3,1,1)
% p11 = pcolor(time,height,Z');
% set(p11, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_height])
% ylabel('Height [km]')
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
% set(gca, 'FontSize',8)
% c=colorbar('north');
% c.Label.String = 'Z [dBZ]';
% c.Label.FontSize = 8;
% % xlabel('Time [UTC]')
% t = title('Radar Reflectivity Factor', 'FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.045;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t
% 
% subplot(3,1,2)
% p12 = pcolor(time,height,beta_atten');
% set(p12, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_height])
% % xlabel('Time [UTC]')
% ylabel('Height [km]')
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [1 1.5 2.5]); 
% set(gca, 'FontSize',8)
% t = title('Attenuated Backscatter Coefficient', 'FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% caxis([0 5*1.e-5])
% c=colorbar('north');
% c.Label.String = '\beta [m^{-1} sr^{-1}]';
% set(c,'YTick',1*1.e-5:1*1.e-5:4*1.e-5)
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%2.e')) )
% c.Label.FontSize = 8;
% grid on
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.045;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t
% 
% subplot(3,1,3)
% plot(time,lwp,'.-')
% hold on
% plot([time(1) time(end)],[min_lwp min_lwp],'r-')
% hold on
% plot([time(1) time(end)],[0 0],'k-')
% axis tight
% hold on
% plot([time(1) time(end)],[max_lwp max_lwp],'r-')
% ylabel('LWP [g m^{-2}]')
% xlabel('Time [UTC]')
% set(gca, 'FontSize',8)
% set(gca, 'Color', 'none')
% t = title('Liquid Water Path', 'FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% clear t
% 
% 
% fig_name = ([num2str(year) num2str(month) num2str(day) '_data_timeseries']);
% export_fig(sprintf(fig_name), '-pdf', '-eps', '-transparent')
% savefig(fig_name)
% clear fig_name ax cpos 

%% Plot time series of retrievals
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) ' retrieval time series'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'units','centimeters','Position',[2 50 8.3 10]);
% subplot(2,1,1)
% p11 = pcolor(time,height,re_hm');
% set(p11, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_height])
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
% set(gca,'XTick',time_start+2:2:time_end-2)
% ylabel('Height [km]')
% set(gca, 'FontSize',8)
% c=colorbar('north');
% set(c,'YTick',4:2:12)
% c.Label.String = 'r_e [\mum]';
% c.Label.FontSize = 8;
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.06;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c 
% t = title('Cloud Droplets Effective Radius','FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% clear t
% 
% subplot(2,1,2)
% plot(time,N_hm,'*','markers',1.5)
% ylabel('N_d [cm^{-3}] * 10^3')
% xlabel('Time [UTC]')
% set(gca,'xlim',[time_start time_end])
% set(gca,'YTick',500:1000:2500)
% set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
% set(gca,'XTick',time_start+2:2:time_end-2)
% set(gca, 'FontSize',8)
% t = title('Cloud Droplets Number Concentration','FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% clear t
% 
% 
% fig_name = ([num2str(year) num2str(month) num2str(day) '_wcp_timeseries']);
% export_fig(sprintf(fig_name), '-pdf', '-eps', '-transparent')
% savefig(fig_name)
% clear fig_name

%% Plot scatter plot of ATB and Z
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) '-' 'correlations beta vs Z'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[2 30 15 15]);
% for i = 1:n_bins
%  % Plot the scatter plot
%  cx(i) = subplot(sx,sy,i) ; 
%  scatter(data_lidar_corr{i}, data_radar_corr{i}, s,  data_lwp_corr{i} , 'filled')
%  hold on
%  corr_lidar_radar = fit((data_lidar_corr{i}'),(data_radar_corr{i}'),'poly1');
%  mdl_lidar_radar{i} = fitlm(log(data_lidar_corr{i}'),log(data_radar_corr{i}'));
%  plot(corr_lidar_radar)
%  legend('off')
%  hold off
%  set(gca, 'YScale', 'log')
%  set(gca, 'XScale', 'log')
%  set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
%  set(gca, 'FontSize',8)
%  caxis([c_min c_max])
%  xlim([0.5*1.e-3 3*1.e-3])
%  ylim([1*1.e-20 1*1.e-18])
%  a=get(gca,'XLim');
%  x=max(a)-(max(a)-1.1*min(a));
%  b=get(gca,'YLim');
%  y=max(b)-(max(b)-min(b))/5;
%  text(x,y,[' \itr = ',num2str(correlations_Z_beta{i}(2),'%.2f'),...
%   '\newline','\itr^2 = ',num2str(mdl_lidar_radar{i}.Rsquared.Ordinary,'%.2f')])%, ... 
%   ...  %'\newline','\itm = ',num2str(poly_Z_beta{i}(1),'%.2e'), ...
% %     '\newline','\itACI_Z = ',num2str((poly_Z_aci{i}(1))*-1,'%.2f')],'FontSize',10);
%  title([num2str(round(min_lwp + (bin_size*(i-1)))) ' < LWP < ',  ... 
%    num2str(round(min_lwp + (bin_size*i)))],'FontSize',10,'FontWeight','normal')
% end
% H2=labelEdgeSubPlots('\beta [sr^{-1}* 10^{-3}]','Z [m^6 m^{-3} * m]');
% h2=colorbar;
% set(h2, 'Position', [.925 .11 .0181 .81])
% h2.Label.String = 'LWP \newline [g m^{-2}]';
% h2.Label.FontSize = 8;
% h2.Label.Rotation = 0; 
% h2.Label.Position = [2 15 0];

% fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_Z_' ...
%     num2str(round(bin_size))]);
% export_fig(sprintf(fig_name), '-pdf', '-eps', '-transparent')
% savefig(fig_name)
% clear TitleFigure xlim ylim x y a b fig_name h2

%% Plot scatter plot of ATB and Reff
TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'correlations beta vs Reff'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units','centimeters','Position',[10 30 15 15]);
for i = 1:n_bins
 cx(i) = subplot(sx,sy,i) ; 
 scat_fig_R(i) = scatter((data_lidar_corr{i}),(data_reff_corr{i}), s,  data_lwp_corr{i} , 'fill');
 hold on
 [corr_lidar_reff, gof_reff] = fit((data_lidar_corr{i}'),(data_reff_corr{i}'),'poly1');
 mdl_lidar_reff{i} = fitlm(log(data_lidar_corr{i}'),log(data_reff_corr{i}'));
 plot(corr_lidar_reff)
 legend('off')
 set(gca, 'YScale', 'log')
 set(gca, 'XScale', 'log')
 set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
 set(gca, 'FontSize',8)
 caxis([c_min c_max])
 xlim([0.5*1.e-3 3*1.e-3])
 ylim([1 100])
 a=get(gca,'XLim');
 x=max(a)-(max(a)-3.25*min(a));
 b=get(gca,'YLim');
 y=max(b)-(max(b)-min(b))/1.5;
 text(x,y,['\it r = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
    '\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f')],...
    'FontSize',10);
    ...%'\newline','\itm = ',num2str(poly_reff_beta{i}(1),'%.2f') ...
%     '\newline','\itACI_{r} = ',num2str((poly_reff_aci{i}(1))*-1,'%.2f')]
 title([num2str(round(min_lwp + (bin_size*(i-1)))) ' < LWP < ',  ... 
   num2str(round(min_lwp + (bin_size*i)))],'FontSize',10,'FontWeight','normal')
end
H=labelEdgeSubPlots('\beta [sr^{-1} * 10^{-3}]','r_e [\mum]');
h2=colorbar;
set(h2, 'Position', [.925 .11 .0181 .81])
h2.Label.String = 'LWP \newline [g m^{-2}]';
h2.Label.FontSize = 8;
h2.Label.Rotation = 0; 
h2.Label.Position = [2 15 0];
% 
fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_r_' ...
    num2str(round(bin_size))]);
export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
clear TitleFigure xlim ylim x y a b fig_name fig h2
% 
%% Plot the scatter plot of ATB and N_d - ACI_N is not dependent on LWP, just one plot
  
[correlations_Nd_beta, pvalue_Nd_beta]= ...
      corrcoef(log(data_lidar_level(id_cloud)),log(data_Nd(id_cloud)));
[poly_Nd_beta, stuct_Nd_beta, mu_Nd_beta] = ...
      polyfit(data_lidar_level(id_cloud), data_Nd(id_cloud)',1);
[poly_Nd_aci, stuct_Nd_aci] = ...
      polyfit(log(data_lidar_level(id_cloud)), log(data_Nd(id_cloud)'),1);

TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'correlations beta vs N_d'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 8.3 8.3]);
scatter(data_lidar_level(id_cloud), data_Nd(id_cloud), s, lwp_cloud(id_cloud) , 'fill')
hold on
[corr_lidar_Nd, gof_Nd] = fit(data_lidar_level(id_cloud)',data_Nd(id_cloud),'poly1');
mdl_lidar_Nd = fitlm(log(data_lidar_level(id_cloud)),log(data_Nd(id_cloud)));
plot(corr_lidar_Nd)
legend('off')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
set(gca,'FontSize',8)
xlim([0.5*1.e-3 3*1.e-3])
ylim([100 2500])
a=get(gca,'XLim');
x=max(a)-(max(a)-3.5*min(a));
b=get(gca,'YLim');
y=max(b)-(max(b)-min(b))/1.025;
text(x,y,[ '\it r = ',num2str(correlations_Nd_beta(2),'%.2f'), ...
    '\newline','\itr^2 = ',num2str(mdl_lidar_Nd.Rsquared.Ordinary,'%.2f')], ...
    'FontSize',10);
    ...%'\newline','\itm = ',num2str(poly_Nd_beta{i}(1),'%.2e') ...
%     '\newline','\itACI_{N} = ',num2str(poly_Nd_aci(1),'%.2f')],'FontSize',10);
caxis([c_min c_max])
title([num2str(min_lwp) ' < LWP < ', num2str(max_lwp)],'FontSize',10, ...
     'FontWeight','normal')
xlabel('\beta [sr^{-1} * 10^{-3}]')
ylabel('N_d [cm^{-3}] ')
h2=colorbar;
xlabel(h2,'LWP [g m^{-2}]')
set(h2,'Location','southoutside');
ax = gca;
axpos = ax.Position;
cpos = h2.Position;
cpos(4) = 0.5*cpos(4);
h2.Position = cpos;
ax.Position = axpos;
% 
fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_N']);
export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
clear TitleFigure xlim ylim x y a b fig_name h2

%% Statistics

% for i = 1:n_bins
% 
% stats.ACI_r(i) = struct('r', {(correlations_reff_beta{i}(2))}, ...
%     'rsq', {(rsq_reff_beta{i})}, ...
%     'm', {poly_reff_beta{i}(1)}, ...
%     'ACI_r', {poly_reff_aci{i}(1)*-1});
% stats.ACI_Z(i) = struct('r', {(correlations_Z_beta{i}(2))}, ...
%     'rsq', {(rsq_Z_beta{i})}, ...
%     'm', {poly_Z_beta{i}(1)}, ...
%     'ACI_Z', {poly_Z_aci{i}(1)*-1});
% end
% stats.ACI_N = struct('r', {(correlations_Nd_beta(2))}, ...
%     'rsq', {(rsq_Nd_beta)}, ...
%     'm', {poly_Nd_beta(1)}, ...
%     'ACI_N', {poly_Nd_aci(1)*-1});
% 
% %% Write statistics to a table
% ACI_r = struct2table(stats.ACI_r);
% writetable(ACI_r, [ num2str(year) num2str(month) num2str(day) ...
%     'ACI_r_statistics_bin' num2str(round(bin_size))]);
% ACI_N = struct2table(stats.ACI_N);
% writetable(ACI_N, [ num2str(year) num2str(month) num2str(day) ...
%     'ACI_N_statistics_bin' num2str(round(bin_size))]);
% ACI_Z = struct2table(stats.ACI_Z);
% writetable(ACI_Z, [ num2str(year) num2str(month) num2str(day) ...
%     'ACI_Z_statistics_bin' num2str(round(bin_size))]);

%% Plot the scatter plot of ATB and extinction
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) '-' 'correlations beta vs extinction'];
% figure('name', TitleFigure, 'NumberTitle','off');
% for i = 1:n_bins
%  % Plot the scatter plot
%  cx(i) = subplot(sx,sy,i) ; 
%  scatter(data_lidar_corr{i}, data_ext_corr{i}, s,  data_lwp_corr{i} , 'fill')
%  hold on
%  corr_lidar_radar = fit(data_lidar_corr{i}', data_radar_corr{i}','poly1');
%  plot(corr_lidar_radar)
%  legend('off')
%  hold off
%  hold on
%  lsline
%  set(gca, 'YScale', 'log')
%  set(gca, 'XScale', 'log')
%  caxis([c_min c_max])
%  xlim([0.5*1.e-3 3*1.e-3])
%  ylim([0.01 0.3])
%  a=get(gca,'XLim');
%  x=max(a)-(max(a)-min(a))/2;
%  b=get(gca,'YLim');
%  y=max(b)-(max(b)-min(b))/2;
%  text(x,y,['\rho = ',num2str(correlations_ext_beta{i}(2),'%.2f'),...
%     '\newline','slope = ',num2str(poly_ext_beta{i}(1),'%.2e')]);
%  title(['LWP ' num2str(round(min_lwp + (bin_size*(i-1)))) '-',  ... 
%    num2str(round(min_lwp + (bin_size*i)))  ' g/m^2'],'FontSize',16,'FontWeight','bold')
%  xlabel('\beta [sr^{-1}]')
%  ylabel('Cloud Extinction [1/m]')
%  set(gca, 'FontSize',8)
% end
% h2=colorbar('FontSize',8);
% set(h2, 'Position', [.93 .11 .0181 .8150])
% xlabel(h2,'g/m^2')
% clear TitleFigure xlim ylim x y a b

%% Plot the scatter plot of ATB and LWC
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) '-' 'correlations beta vs LWC'];
% figure('name', TitleFigure, 'NumberTitle','off');
% for i = 1:n_bins
%  % Plot the scatter plot
%  cx(i) = subplot(sx,sy,i) ; 
%  scatter(data_lidar_corr{i}, data_lwc_corr{i}, s,  data_lwp_corr{i} , 'fill')
%  hold on
%  corr_lidar_radar = fit(data_lidar_corr{i}', data_radar_corr{i}','poly1');
%  plot(corr_lidar_radar)
%  legend('off')
%  hold off
%  hold on
%  lsline
%  set(gca, 'YScale', 'log')
%  set(gca, 'XScale', 'log')
%  caxis([c_min c_max])
%  xlim([0.5*1.e-3 3*1.e-3])
%  ylim([0.01 0.8])
%  %  a=get(gca,'XLim');
%  x=max(a)-(max(a)-min(a))/2;
%  b=get(gca,'YLim');
%  y=max(b)-(max(b)-min(b))/2;
%  text(x,y,['\rho = ',num2str(correlations_lwc_beta{i}(2),'%.2f'),...
%     '\newline','slope = ',num2str(poly_lwc_beta{i}(1),'%.2e')]);
%  title(['LWP ' num2str(round(min_lwp + (bin_size*(i-1)))) '-',  ... 
%    num2str(round(min_lwp + (bin_size*i)))  ' g/m^2'],'FontSize',16,'FontWeight','bold')
%  xlabel('\beta [sr^{-1}]')
%  ylabel('LWC [g/m^3]')
%   set(gca, 'FontSize',8)
% end
% h2=colorbar('FontSize',8);
% set(h2, 'Position', [.93 .11 .0181 .8150])
% xlabel(h2,'g/m^2')
% clear TitleFigure xlim ylim x y a b

%% Scatter plot of ATB and tau
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) '-' 'correlations beta vs optical thickness'];
% figure('name', TitleFigure, 'NumberTitle','off');
% for i = 1:n_bins
%  % Plot the scatter plot
%  cx(i) = subplot(sx,sy,i) ; 
%  scatter(data_lidar_corr{i}, data_tau_corr{i}, s,  data_lwp_corr{i} , 'fill')
%  str1 = [' \rho = ' num2str(correlations_tau_beta{i}(2),'% 1.2f')];
%  text(0.00055, 50, str1);
%  str2 = [' slope = ' num2str(poly_tau_beta{i}(1),'% 4.3e')];
%  text(0.00055, 45, str2);
%  hold on
%  lsline
%  set(gca, 'YScale', 'log')
%  set(gca, 'XScale', 'log')
%  xlim([0.5*1.e-3 3*1.e-3])
%  ylim([0 50])
%  caxis([c_min c_max])
%  title(['LWP ' num2str(round(min_lwp + (bin_size*(i-1)))) '-',  ... 
%    num2str(round(min_lwp + (bin_size*i)))  ' g/m^2'],'FontSize',16,'FontWeight','bold')
%  xlabel('Attenuated Backscatter Coefficient [1/m/sr]')
%  ylabel('\tau')
%  set(gca, 'FontSize',8)
% end
% h2=colorbar('FontSize',16);
% set(h2, 'Position', [.93 .11 .0181 .8150])
% xlabel(h2,'g/m^2')
% clear TitleFigure xlim ylim

% % Plot time series of Z, Z_water, ATB and ATB selected
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) ' selected data'];
% figure('name', TitleFigure, 'NumberTitle','off');
% subplot(2,2,1)
% cmap=colormap;
% cmap(1,:)=1;
% p2 = pcolor(time,height,Z_water');
% set(p2, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_height])
% colormap(cmap);
% c=colorbar;
% ylabel(c,'Z [dBZ]')
% hold on
% p22 = plot(time,cb_layer,'k.');
% hold on
% p23 = plot(time,ct_layer,'.');
% hold on
% p24 = plot(time(is_data),height(aerosols(is_data)),'m.');
% hold on
% p25 = plot(time(is_data),height(cloud(is_data)),'g.');
% legend([p22 p23 p24 p25], {'cb layer','ct layer', 'aerosols level', ...
%     'cloud level'})
% ylabel('Height above MSL [m]')
% xlabel('Time [UTC]')
% title('Z single-layer, no drizzle cloud')
% grid on
% 
% subplot(2,2,2)
% cmap=colormap;
% cmap(1,:)=1;
% p5 = pcolor(time,height,re_hm');
% set(p5, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_height])
% hold on
% p51 = plot(time(is_data),height(cloud(is_data)),'g.');
% legend([p51], {'cloud level'})
% colormap(cmap);
% c=colorbar;
% ylabel(c,' R_e_f_f [\mum]')
% ylabel('Height above MSL [m]')
% xlabel('Time [UTC]')
% title('Cloud droplets Effective Radius (R_e_f_f)')
% grid on
% 
% subplot(2,2,3)
% p4 = pcolor(time,height,int_beta_atten');
% set(p4, 'EdgeColor', 'none');
% axis tight
% xlabel('Time [UTC]')
% ylabel('Height above MSL [m]')
% title('Integrated Attenuated Backscatter Coefficient (\beta)')
% grid on
% hold on
% p41 = plot(time(is_data),height(aerosols(is_data)),'m.');
% hold on
% p42 = plot(time(is_data),height(cloud(is_data)),'g.');
% legend([p41 p42], {'aerosols level','cloud level'})
% colormap(cmap);
% c=colorbar;
% ylabel(c,'\beta [sr^{-1}]')
% 
% subplot(2,2,4)
% plot(time,lwp,'.-')
% hold on
% plot([time(1) time(end)],[min_lwp min_lwp],'r-')
% hold on
% plot([time(1) time(end)],[0 0],'k-')
% axis tight
% hold on
% plot([time(1) time(end)],[max_lwp max_lwp],'r-')
% ylabel('LWP [g/m^2]')
% xlabel('Time [UTC]')
% title('LWP')
% grid on
% clear TitleFigure
% 

%Plot the cloudnet categorisation
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' num2str(year)];
% figure('name', TitleFigure, 'NumberTitle','off');
% p5 = pcolor(time,height,category_bits');
% set(p5, 'EdgeColor', 'none');
% xlabel('Time [UTC]')
% ylabel('Height above MSL [m]')
% title('CLOUDNET target categorisation')
% caxis([0 5])
% 
% % legend([p51 p52], {'aerosols level','cloud level'})
% % mymap = [0.5, 0.1, 0.5
% %     0.1, 0.5, 0.8
% %     0.2, 0.7, 0.6
% %     0.8, 0.7, 0.3
% %     0.9, 1, 0];
% % colormap(mymap);
% bar_labels = {'cloud','rain','drizzle','aerosols','ice'};
% lcolorbar(bar_labels);
% % ylabel(c,'category bits')
% clear TitleFigure

% % Plot time series of LWP, N_d and selected height for Z, ATB and R_eff
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) 'retrievals'];
% figure('name', TitleFigure, 'NumberTitle','off');
% subplot(5,1,1)
% plot(time(is_data),lwp(is_data),'.-')
% hold on
% plot([time(is_data(1)) time(is_data(end))],[min_lwp min_lwp],'r-')
% hold on
% plot([time(is_data(1)) time(is_data(end))],[0 0],'k-')
% axis tight
% hold on
% plot([time(is_data(1)) time(is_data(end))],[max_lwp max_lwp],'r-')
% ylabel('LWP [g/m^2]')
% xlabel('Time [UTC]')
% title('LWP')
% grid on
% 
% subplot(5,1,2)
% plot(time(is_data),data_radar_level)
% axis tight
% xlabel('Time [UTC]')
% title('Integrated Radar Reflectivity [m^6/m^3 * m])')
% grid on
% 
% subplot(5,1,3)
% plot(time(is_data),data_lidar_level)
% axis tight
% xlabel('Time [UTC]')
% title('Integrated Attenuated Backscatter [[sr^{-1}]]')
% grid on
% 
% subplot(5,1,4)
% plot(time(is_data),data_reff_level)
% axis tight
% xlabel('Time [UTC]')
% title('Cloud Droplet Effective Radius (r_e)')
% grid on
% 
% subplot(5,1,5)
% plot(time,N_hm)
% axis tight
% xlabel('Time [UTC]')
% title('Cloud Droplet Number Concentration N_d')
% grid on
% 
% clear TitleFigure


