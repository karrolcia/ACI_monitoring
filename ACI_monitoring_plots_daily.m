%% Make a loop over the fields in cloudnet.(date)
%% ACI Scheme - plots only
%ACI_scheme_plotsv2
max_height = 3000;
%% Create scatter plots of the selected data
% Create scatter plots based on the LWP interval
s = 20.0 ;
% Select only the data points that where retrieved with the microphysical
% retrieval 
    
%% Define level to compare
is_data = find(nansum(retrieval_knist.(date).re_hm,2) ~= 0); % find only timesteps with retrieved data
% define level of aerosols to compare
for it = 1:n_time
aerosol(it) = I_max(it)- 9; 
cloud(it) = I_max(it)+2;
end

% define level of cloud to compare

aero_level = aerosol(is_data) ; % Attenuated Backscatter
aero_level = aero_level(~isnan(aero_level));

cloud_level = cloud(is_data) ;  % Radar Reflectivity Factor
cloud_level = cloud_level(~isnan(cloud_level));
aero_level = aero_level(~isnan(cloud_level));
cloud_level = cloud_level(~isnan(aero_level));
is_data = is_data(~isnan(aero_level));
is_data = is_data(~isnan(cloud_level));


data_lidar = (int_beta_atten(is_data,:));
data_radar = (int_Z(is_data,:));
data_doppler = cloudnet.(date).v(is_data,:);
data_reff = (retrieval_knist.(date).re_hm(is_data,:));
data_ext_cloud = retrieval_knist.(date).ext_hm(is_data,:);
data_lwc_cloud = retrieval_knist.(date).lwc_hm(is_data,:);
data_Nd = retrieval_knist.(date).N_hm(is_data);
data_tau_cloud = retrieval_knist.(date).tau_hm(is_data); 

% Meteorological information
data_pressure = cloudnet.(date).pressure(is_data,:);
data_temperature = cloudnet.(date).temperature(is_data,:);
data_humidity = cloudnet.(date).specific_humidity(is_data,:);


for i = 1:length(is_data)
data_lidar_level(i) = data_lidar(i,aero_level(i));
data_radar_level(i) = data_radar(i,cloud_level(i));
data_doppler_level(i) = data_doppler(i,cloud_level(i));
data_reff_level(i) = data_reff(i,cloud_level(i));
data_ext_level(i) = data_ext_cloud(i,cloud_level(i));
data_lwc_level(i) = data_lwc_cloud(i,cloud_level(i));

data_pressure_interval(i) = data_pressure(i,aero_level(i));
data_temperature_interval(i) = data_temperature(i,aero_level(i));
data_humidity_interval(i) = data_humidity(i,aero_level(i)); 
end
data_reff_level(data_reff_level>=7)=NaN;
data_reff_level = data_reff_level(~isnan(data_lidar_level));
data_reff_level = data_reff_level(~isnan(data_reff_level));
data_lidar_level = data_lidar_level(~isnan(data_lidar_level));
data_lidar_level = data_lidar_level(~isnan(data_reff_level));

data_pressure_interval = data_pressure_interval(~isnan(data_reff_level));
data_pressure_interval = data_pressure_interval(~isnan(data_lidar_level));
data_temperature_interval = data_temperature_interval(~isnan(data_reff_level));
data_temperature_interval = data_temperature_interval(~isnan(data_lidar_level));
data_humidity_interval = data_humidity_interval(~isnan(data_reff_level));
data_humidity_interval = data_humidity_interval(~isnan(data_lidar_level));

%% Define different LWP intervals
lwp_cloud = cloudnet.(date).lwp(is_data);
lwp_cloud = lwp_cloud(~isnan(data_reff_level));
lwp_cloud = lwp_cloud(~isnan(data_lidar_level));

% lwp_cloud = lwp_cloud(~isnan(data_reff_level));
if min(lwp_cloud) < 30
    min_lwp =30;
else min_lwp = min(lwp_cloud) ;
end
if max(lwp_cloud) > 150
max_lwp =150;
else max_lwp = max(lwp_cloud);
end
max_lwp =90;
id_cloud = find(lwp_cloud > min_lwp & lwp_cloud < max_lwp );

%% Define subplot net
sx = 3; % number of rows
sy = 2; %number of columns
n_bins = sx*sy; % must be divided by 2
bin_size =  round((max_lwp - min_lwp)/n_bins);
c_min = min_lwp;
c_max = max_lwp;

for i = 1:n_bins
  
 id_lwp = find(lwp_cloud >= min_lwp + bin_size*(i-1) & ...
         lwp_cloud <= min_lwp + bin_size*i);
 lwp_id{i} = id_lwp;
 % Divide data into bins based on the LWP
    for j = 1:length(lwp_id{i})
      corr_lidar(i,j) = data_lidar_level(lwp_id{i}(j));
      corr_radar(i,j) = data_radar_level(lwp_id{i}(j));
      corr_reff(i,j) = data_reff_level(lwp_id{i}(j));
    end
    
 % Calculate correlation coefficient for each scatter plot
  data_lidar_corr{i} = corr_lidar(i,:);
  data_radar_corr{i} = corr_radar(i,:);
  
  data_reff_corr{i} = corr_reff(i,:);
  data_reff_corr{i} = data_reff_corr{i}(~isnan(data_reff_corr{i}));
  data_lidar_corr{i} = data_lidar_corr{i}(~isnan(data_reff_corr{i}));
  data_lidar_corr{i} = data_lidar_corr{i}(~isnan(data_lidar_corr{i}));
  data_reff_corr{i} = data_reff_corr{i}(~isnan(data_lidar_corr{i}));
  
  data_radar_corr{i} = data_radar_corr{i}(~isnan(data_reff_corr{i}));
  data_Nd_corr{i} = (data_Nd(lwp_id{i})) ;
  data_Nd_corr{i} = data_Nd_corr{i}(~isnan(data_reff_corr{i}));
  data_Nd_corr{i} = data_Nd_corr{i}(~isnan(data_lidar_corr{i}));
  data_lwp_corr{i} = (lwp_cloud(lwp_id{i})) ;
  data_lwp_corr{i} = data_lwp_corr{i}(~isnan(data_reff_corr{i}));
  data_lwp_corr{i} = data_lwp_corr{i}(~isnan(data_lidar_corr{i}));
  
  data_time{i} = (time(lwp_id{i}));
  data_time{i} = data_time{i}(~isnan(data_reff_corr{i}));
  % Calculate correlations

  [correlations_reff_beta{i}, pvalue_reff_beta{i}] = ...
      corrcoef(log(data_lidar_corr{i}),...
      log(data_reff_corr{i}));
  [poly_reff_aci{i}, stuct_aci_reff{i},mu_reff_aci{i}] = ...
      polyfit(log(data_lidar_corr{i}), log(data_reff_corr{i}),1);
 
  ACI_corr.(date).ACI_reff{i} = -1*poly_reff_aci{i}(1);
clear id_lwp corr_* 
end

%% Plot scatter plot of ATB and Reff
TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'correlations beta vs Reff'];

figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units','centimeters','Position',[10 30 20 15]);
    for i = 1:n_bins
     cx(i) = subplot(sx,sy,i) ; 
     scat_fig_R(i) = scatter((data_lidar_corr{i}),(data_reff_corr{i}), s, ...
         data_lwp_corr{i} , 'fill');
     hold on
     [corr_lidar_reff, gof_reff] = fit(log(data_lidar_corr{i}'),log(data_reff_corr{i}'),'poly1');
     mdl_lidar_reff{i} = fitlm(log(data_lidar_corr{i}'),log(data_reff_corr{i}'));
%      if mdl_lidar_reff{i}.Coefficients.pValue(2) <= 0.1
%         plot(corr_lidar_reff)
%      end   
     legend('off')
%      set(gca, 'YScale', 'log')
     set(gca, 'XScale', 'log')
%      set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
     set(gca, 'FontSize',8)
     caxis([c_min c_max])
%      xlim([-8 -6])
%      ylim([-3 3])
     xlim([1.5*1.e-4 2.5*1.e-3])
     ylim([1 10])
     a=get(gca,'XLim');
     x=max(a)-(max(a)-1.1*min(a));
     b=get(gca,'YLim');
     y=max(b)-(max(b)-min(b))/1.25;
     text(x,y,['\it R = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
      ...%'\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
      '\newline','\itACI_{r} = ',num2str(ACI_corr.(date).ACI_reff{i},'%.2f')],...
      'FontSize',10)
%       '\newline','\itp-value = ',num2str(pvalue_reff_beta{i}(2),'%.2e')], ...
%       'FontSize',10);
     ACI_corr.(date).reff_r{i} = correlations_reff_beta{i}(2);
     ACI_corr.(date).reff_r2{i} = mdl_lidar_reff{i}.Rsquared.Ordinary;
%
     title([num2str(round(min_lwp + (bin_size*(i-1)))) ' < LWP < ',  ... 
       num2str(round(min_lwp + (bin_size*i)))],'FontSize',10,'FontWeight','normal')
    end
    H=labelEdgeSubPlots('\beta [sr^{-1}]','r_e [\mum]');
    h2=colorbar;
    set(h2, 'Position', [.925 .11 .0181 .81])
    h2.Label.String = 'LWP \newline [g m^{-2}]';
    h2.Label.FontSize = 8;
    h2.Label.Rotation = 0; 
    h2.Label.Position = [2 15 0];
% % 
fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_r_' ...
     num2str(round(bin_size))]);
export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
clear TitleFigure xlim ylim x y a b fig_name fig h2
% date
for ii=1:n_bins
    fprintf('%8.0f < LWP < %8.0f & %8.2f & %8.2f & %8.2f & %8.0f \\\\ \n', ...
        round(min_lwp + (bin_size*(ii-1))), ...
        round(min_lwp + (bin_size*ii)), ...
        (poly_reff_aci{ii}(1)*-1), ...
        correlations_reff_beta{ii}(2), ...    
        mdl_lidar_reff{ii}.Rsquared.Ordinary, ...
        numel(data_lidar_corr{ii}))
         n_data(ii) =   numel(data_lidar_corr{ii});
end 
data_points = sum(n_data)
%% Plot the scatter plot of ATB and N_d - ACI_N is not dependent on LWP, 
%  just one plot
  
[correlations_Nd_beta, pvalue_Nd_beta]= ...
      corrcoef(log(data_lidar_level(id_cloud)),log(data_Nd(id_cloud)));
[poly_Nd_aci, stuct_Nd_aci] = ...
      polyfit(log(data_lidar_level(id_cloud)), log(data_Nd(id_cloud)),1);

TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'correlations beta vs N_d'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 8.3 8.3]);
    scatter(data_lidar_level(id_cloud), data_Nd(id_cloud), s, lwp_cloud(id_cloud) , 'fill')
    hold on
    % [corr_lidar_Nd, gof_Nd] = fit(data_lidar_level(id_cloud)',data_Nd(id_cloud),'poly1');
    mdl_lidar_Nd = fitlm(log(data_lidar_level(id_cloud)),log(data_Nd(id_cloud)));
    % plot(corr_lidar_Nd)
    legend('off')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
%     set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
    set(gca,'FontSize',8)
    xlim([1.5*1.e-4 2.5*1.e-3])
    ylim([10 2500])
     a=get(gca,'XLim');
     x=max(a)-(max(a)-1.1*min(a));
     b=get(gca,'YLim');
     y=max(b)-(max(b)-min(b))/1.01;
    text(x,y,[ '\itR = ',num2str(correlations_Nd_beta(2),'%.2f'), ...
        ...%'\newline','\itr^2 = ',num2str(mdl_lidar_Nd.Rsquared.Ordinary,'%.2f'), ...
        '\newline','\itACI_{N} = ',num2str(poly_Nd_aci(1),'%.2f')], ...
        'FontSize',10);
    caxis([c_min c_max])
    title([num2str(round(min_lwp)) ' < LWP < ', num2str(round(max_lwp))],'FontSize',10, ...
         'FontWeight','normal')
    xlabel('\beta [sr^{-1}]')
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
    ACI_corr.(date).ACI_N = poly_Nd_aci(1);
    ACI_corr.(date).N_r = correlations_Nd_beta(2);
    ACI_corr.(date).N_r2 = mdl_lidar_Nd.Rsquared.Ordinary;
% % 
fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_N']);
export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
% 
% ACI = cell2mat(ACI_corr.(date).ACI_reff);
% reff_r = cell2mat(ACI_corr.(date).reff_r);
% reff_r2 = cell2mat(ACI_corr.(date).reff_r2);
% lwp_s = min_lwp:bin_size:max_lwp;
% TitleFigure=('ACI vs. Correlation Coefficient');
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 50 15 15]);
% scatter(ACI, reff_r,50,lwp_s(1:end-1),'fill')
% caxis([min(lwp_s) max(lwp_s)])
% lsline
% title('ACI vs. Correlation Coefficient','FontSize',10, ...
%          'FontWeight','normal')
% xlabel('ACI')
% ylabel('r ')
% h2=colorbar;
% xlabel(h2,'LWP [g m^{-2}]')
%     set(h2,'Location','southoutside');
%     ax = gca;
%     axpos = ax.Position;
%     cpos = h2.Position;
%     cpos(4) = 0.5*cpos(4);
%     h2.Position = cpos;
%     ax.Position = axpos;


%% Statistics

% for i = 1:n_bins
% 
% stats.ACI_r(i) = struct('r', {(correlations_reff_beta{i}(2))}, ...
%     'rsq', {(rsq_reff_beta{i})}, ...
%     'm', {poly_reff_beta{i}(1)}, ...
%     'ACI_r', {poly_reff_aci{i}(1)*-1});
% % stats.ACI_Z(i) = struct('r', {(correlations_Z_beta{i}(2))}, ...
% %     'rsq', {(rsq_Z_beta{i})}, ...
% %     'm', {poly_Z_beta{i}(1)}, ...
% %     'ACI_Z', {poly_Z_aci{i}(1)*-1});
% end
% stats.ACI_N = struct('r', {(correlations_Nd_beta(2))}, ...
%     'rsq', {(rsq_Nd_beta)}, ...
%     'm', {poly_Nd_beta(1)}, ...
%     'ACI_N', {poly_Nd_aci(1)*-1});

%% Write statistics to a table
% ACI_r = struct2table(stats.ACI_r);
% writetable(ACI_r, [ num2str(year) num2str(month) num2str(day) ...
%     'ACI_r_statistics_bin' num2str(round(bin_size))]);
% ACI_N = struct2table(stats.ACI_N);
% writetable(ACI_N, [ num2str(year) num2str(month) num2str(day) ...
%     'ACI_N_statistics_bin' num2str(round(bin_size))]);
% ACI_Z = struct2table(stats.ACI_Z);
% writetable(ACI_Z, [ num2str(year) num2str(month) num2str(day) ...
%     'ACI_Z_statistics_bin' num2str(round(bin_size))]);

%% Plot time series of measurements
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) ' data time series'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'units','centimeters','Position',[2 50 8.3 15]);
% subplot(3,1,1)
%     p11 = pcolor(time,height,cloudnet.(date).Z');
%     set(p11, 'EdgeColor', 'none');
%     set(gca,'Ydir','normal');
%     ylabel('Height [km]')
%     set(gca, 'FontSize',8)
%     c=colorbar('north');
%     c.Label.String = 'Z [dBZ]';
%     c.Label.FontSize = 8;
% % xlabel('Time [UTC]')
%     t = title('Radar Reflectivity Factor', 'FontSize',10,'FontWeight','normal');
%     set(t, 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
%     h1 = get(t, 'position');
%     set(t, 'position', [0 h1(2) h1(3)])
%     grid on
%     xlim([0 24])
%     ylim([0 3000])
%     ax = gca;
%     axpos = ax.Position;
%     cpos = c.Position;
%     cpos(1) = axpos(1);
%     cpos(2) = cpos(2) + 0.045;
%     cpos(3) = axpos(3);
%     cpos(4) = 0.5*cpos(4);
%     c.Position = cpos;
%     ax.Position = axpos;
%     clear ax cpos c t
% 
%     subplot(3,1,2)
%     p12 = pcolor(time,height,cloudnet.(date).beta');
%     set(p12, 'EdgeColor', 'none');
%     set(gca,'Ydir','normal');
%     % set(gca,'ylim',[0 max_height])
%     % xlabel('Time [UTC]')
%     ylabel('Height [km]')
%     % set(gca,'YTick',500:1000:max_height)
%     % set(gca, 'YTickLabel', [1 1.5 2.5]); 
%     set(gca, 'FontSize',8)
%     t = title('Attenuated Backscatter Coefficient', 'FontSize',10,'FontWeight','normal');
%     set(t, 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
%     h1 = get(t, 'position');
%     set(t, 'position', [0 h1(2) h1(3)])
%     hold on
% %     a = I_max;
% %     a(isnan(a)) = 1;
% %     l3 = plot(time,height(a),'m','LineWidth',1.5);
% %     b = aerosol;
% %     b(isnan(b)) = 1;
% %     l4 = plot(time,height(b),'y','LineWidth',1.5);
% %     c = cloud;
% %     c(isnan(c)) = 1;
% %     l5 = plot(time,height(c),'c','LineWidth',1.5);
%     caxis([0 5*1.e-5])
%     c=colorbar('north');
%     c.Label.String = '\beta [m^{-1} sr^{-1}]';
%     set(c,'YTick',1*1.e-5:1*1.e-5:4*1.e-5)
%     set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%2.e')) )
%     c.Label.FontSize = 8;
%     grid on
%     xlim([0 24])
%     ylim([0 3000])
%     ax = gca;
%     axpos = ax.Position;
%     cpos = c.Position;
%     cpos(1) = axpos(1);
%     cpos(2) = cpos(2) + 0.045;
%     cpos(3) = axpos(3);
%     cpos(4) = 0.5*cpos(4);
%     c.Position = cpos;
%     ax.Position = axpos;
%     clear ax cpos c t
% 
%     subplot(3,1,3)
%     plot(time,cloudnet.(date).lwp,'.-')
%     hold on
% %     plot([time(1) time(end)],[min_lwp min_lwp],'r-')
% %     hold on
% %     plot([time(1) time(end)],[0 0],'k-')
%     axis tight
% %     hold on
% %     plot([time(1) time(end)],[max_lwp max_lwp],'r-')
%     ylabel('LWP [g m^{-2}]')
%     xlabel('Time [UTC]')
%     set(gca, 'FontSize',8)
%     set(gca, 'Color', 'none')
%     t = title('Liquid Water Path', 'FontSize',10,'FontWeight','normal');
%     set(t, 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
%     h1 = get(t, 'position');
%     set(t, 'position', [0 h1(2) h1(3)])
%     grid on
%      xlim([0 24])
%     clear t
% % 
% % fig_name = ([num2str(year) num2str(month) num2str(day) '_data_timeseries']);
% % export_fig(sprintf(fig_name), '-eps', '-transparent')
% % savefig(fig_name)
% % clear fig_name ax cpos 
% % 
%% Plot time series of retrievals
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) ' retrieval time series'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'units','centimeters','Position',[2 50 8.3 10]);
%     subplot(2,1,1)
%     p11 = pcolor(time,height,retrieval_knist.(date).re_hm');
%     set(p11, 'EdgeColor', 'none');
%       xlim([0 24])
%     ylim([0 2000])
%     % set(gca,'Ydir','normal');
%     % set(gca,'ylim',[0 max_height])
%     % set(gca,'YTick',500:1000:max_height)
% %     set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
%     % set(gca,'XTick',time_start+2:2:time_end-2)
%     ylabel('Height [km]')
%     set(gca, 'FontSize',8)
%     c=colorbar('north');
%     set(c,'YTick',4:2:12)
%     c.Label.String = 'r_e [\mum]';
%     c.Label.FontSize = 8;
%     ax = gca;
%     axpos = ax.Position;
%     cpos = c.Position;
%     cpos(1) = axpos(1);
%     cpos(2) = cpos(2) + 0.06;
%     cpos(3) = axpos(3);
%     cpos(4) = 0.5*cpos(4);
%     c.Position = cpos;
%     ax.Position = axpos;
%     clear ax cpos c 
%     t = title('Cloud Droplets Effective Radius','FontSize',10,'FontWeight','normal');
%     set(t, 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
%     h1 = get(t, 'position');
%     set(t, 'position', [0 h1(2) h1(3)])
%     grid on
%     clear t
% 
%     subplot(2,1,2)
%     plot(time,retrieval_knist.(date).N_hm,'*','markers',1.5)
%     ylabel('N_d [cm^{-3}] * 10^3')
%     xlabel('Time [UTC]')
%      xlim([0 24])
% %     ylim([0 2000])
%     % set(gca,'xlim',[time_start time_end])
%     % set(gca,'YTick',500:1000:2500)
% %     set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
%     % set(gca,'XTick',time_start+2:2:time_end-2)
%     set(gca, 'FontSize',8)
%     t = title('Cloud Droplets Number Concentration','FontSize',10,'FontWeight','normal');
%     set(t, 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
%     h1 = get(t, 'position');
%     set(t, 'position', [0 h1(2) h1(3)])
%     grid on
%     clear t
% 
% fig_name = ([num2str(year) num2str(month) num2str(day) '_wcp_timeseries']);
% export_fig(sprintf(fig_name), '-eps', '-transparent')
% % savefig(fig_name)
% clear fig_name

% clear data_* n_data
% [R,P] =  corrcoef(log(data_lidar_level),log(data_reff_level))