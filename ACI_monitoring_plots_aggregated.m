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
is_data_agg = find(nansum(agg_retrival.re_hm,2) ~= 0); % find only timesteps with retrieved data
% define level of aerosols to compare
for it = 1:n_time
aerosols(it) = I_max_agg(it) - 9; 
cloud(it) = I_max_agg(it) +2 ;
end

% define level of cloud to compare

aero_level = aerosols(is_data_agg) ; % Attenuated Backscatter
aero_level = aero_level(~isnan(aero_level));
cloud_level = cloud(is_data_agg) ;  % Radar Reflectivity Factor
cloud_level = cloud_level(~isnan(cloud_level));
aero_level = aero_level(~isnan(aero_level));
aero_level = aero_level(~isnan(cloud_level));
cloud_level = cloud_level(~isnan(aero_level));
is_data_agg = is_data_agg(~isnan(aero_level));
is_data_agg = is_data_agg(~isnan(cloud_level));

data_lidar = (int_beta_atten(is_data_agg,:));
data_radar = (int_Z(is_data_agg,:));
data_doppler = agg_cloudnet.v(is_data_agg,:);
data_reff = (agg_retrival.re_hm(is_data_agg,:));
data_ext_cloud = agg_retrival.ext_hm(is_data_agg,:);
data_lwc_cloud = agg_retrival.lwc_hm(is_data_agg,:);
data_Nd = agg_retrival.N_hm(is_data_agg);
data_tau_cloud = agg_retrival.tau_hm(is_data_agg);

% Find meteo variables in the analysed period
data_pressure = agg_cloudnet.pressure(is_data_agg,:);
data_temperature = agg_cloudnet.temperature(is_data_agg,:);
data_humidity = agg_cloudnet.specific_humidity(is_data_agg,:);


for i = 1:length(is_data_agg)
data_lidar_level_agg(i) = data_lidar(i,aero_level(i));
data_radar_level_agg(i) = data_radar(i,cloud_level(i));
data_doppler_level_agg(i) = data_doppler(i,cloud_level(i));
data_reff_level_agg(i) = data_reff(i,cloud_level(i));
data_ext_level_agg(i) = data_ext_cloud(i,cloud_level(i));
data_lwc_level_agg(i) = data_lwc_cloud(i,cloud_level(i));

data_pressure_interval(i) = data_pressure(i,aero_level(i));
data_temperature_interval(i) = data_temperature(i,aero_level(i));
data_humidity_interval(i) = data_humidity(i,aero_level(i)); 
end

data_lidar_level_agg = data_lidar_level_agg(~isnan(data_lidar_level_agg));
data_reff_level_agg = data_reff_level_agg(~isnan(data_reff_level_agg));
data_lidar_level_agg = data_lidar_level_agg(~isnan(data_reff_level_agg));
data_reff_level_agg = data_reff_level_agg(~isnan(data_lidar_level_agg));
data_pressure_interval = data_pressure_interval(~isnan(data_reff_level_agg));
data_pressure_interval = data_pressure_interval(~isnan(data_lidar_level_agg));
data_temperature_interval = data_temperature_interval(~isnan(data_reff_level_agg));
data_temperature_interval = data_temperature_interval(~isnan(data_lidar_level_agg));
data_humidity_interval = data_humidity_interval(~isnan(data_reff_level_agg));
data_humidity_interval = data_humidity_interval(~isnan(data_lidar_level_agg));

%% Define different LWP intervals
lwp_cloud_agg = agg_cloudnet.lwp(is_data_agg);
lwp_cloud_agg = lwp_cloud_agg(~isnan(data_lidar_level_agg));
lwp_cloud_agg = lwp_cloud_agg(~isnan(data_reff_level_agg));

if min(lwp_cloud_agg) < 30
    min_lwp = 30;
else min_lwp = min(lwp_cloud_agg) ;
end
if max(lwp_cloud_agg) > 90
max_lwp = 90;
else max_lwp = max(lwp_cloud_agg);
end
id_cloud_agg = find(lwp_cloud_agg > min_lwp & lwp_cloud_agg < max_lwp );

data_pressure_interval = data_pressure_interval(id_cloud_agg);
data_temperature_interval = data_temperature_interval(id_cloud_agg);
data_humidity_interval = data_humidity_interval(id_cloud_agg);
%% Define subplot net
sx = 3; % number of rows
sy = 2; %number of columns
n_bins_agg = sx*sy; % must be divided by 2
bin_size =  round((max_lwp - min_lwp)/n_bins_agg);
c_min = min_lwp;
c_max = max_lwp;

for i = 1:n_bins_agg
  
 id_lwp_agg = find(lwp_cloud_agg >= min_lwp + bin_size*(i-1) & ...
         lwp_cloud_agg <= min_lwp + bin_size*i);
 lwp_id_agg{i} = id_lwp_agg;
 % Divide data into bins based on the LWP
    for j = 1:length(lwp_id_agg{i})
      corr_lidar_agg(i,j) = data_lidar_level_agg(lwp_id_agg{i}(j));
      corr_radar_agg(i,j) = data_radar_level_agg(lwp_id_agg{i}(j));
      corr_reff_agg(i,j) = data_reff_level_agg(lwp_id_agg{i}(j));
    end
     
 % Calculate correlation coefficient for each scatter plot
  data_lidar_corr{i} = (corr_lidar_agg(i,:));
  data_radar_corr{i} = (corr_radar_agg(i,:));
  data_reff_corr{i} = (corr_reff_agg(i,:));
  
  data_reff_corr{i} = data_reff_corr{i}(~isnan(data_reff_corr{i}));
  data_lidar_corr{i} = data_lidar_corr{i}(~isnan(data_reff_corr{i}));
  data_lidar_corr{i} = data_lidar_corr{i}(~isnan(data_lidar_corr{i}));
  data_reff_corr{i} = data_reff_corr{i}(~isnan(data_lidar_corr{i}));

  data_Nd_corr{i} = (data_Nd(lwp_id_agg{i})) ;
  data_Nd_corr{i} = data_Nd_corr{i}(~isnan(data_reff_corr{i})); 
  data_Nd_corr{i} = data_Nd_corr{i}(~isnan(data_lidar_corr{i}));
  
  data_lwp_corr{i} = (lwp_cloud_agg(lwp_id_agg{i})) ;
  data_lwp_corr{i} = data_lwp_corr{i}(~isnan(data_reff_corr{i}));
  data_lwp_corr{i} = data_lwp_corr{i}(~isnan(data_lidar_corr{i}));
  
  data_time{i} = (time(lwp_id_agg{i}));
  data_time{i} = data_time{i}(~isnan(data_reff_corr{i}));
  data_time{i} = data_time{i}(~isnan(data_lidar_corr{i}));
  % Calculate correlations

  [correlations_reff_beta_agg{i}, pvalue_reff_beta_agg{i}] = ...
      corrcoef(log(data_lidar_corr{i}),...
     log(data_reff_corr{i}));
  [poly_reff_aci_agg{i}, stuct_aci_beta_agg{i}] = ...
      polyfit(log(data_lidar_corr{i}), log(data_reff_corr{i}),1);
  
  agg_ACI_corr.ACI_reff{i} = -1*poly_reff_aci_agg{i}(1);
clear id_lwp corr_*
end


% fig_name = ([num2str(year) num2str(month) num2str(day) '_wcp_timeseries']);
% export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
% clear fig_name

%% Plot scatter plot of ATB and Reff
TitleFigure=['Cabauw ACCEPT aggregated correlations beta vs Reff'];

figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units','centimeters','Position',[10 30 15 15]);
    for i = 1:n_bins_agg
     cx(i) = subplot(sx,sy,i) ; 
     scat_fig_R(i) = scatter(data_lidar_corr{i},(data_reff_corr{i}), s, ...
         data_lwp_corr{i} , 'fill');
     hold on
%     [corr_lidar_reff, gof_reff] = fit((data_lidar_corr{i}'),(data_reff_corr{i}'),'poly1');
     mdl_lidar_reff{i} = fitlm(log(data_lidar_corr{i}'),log(data_reff_corr{i}'));
%      plot(corr_lidar_reff)
     legend('off')
     set(gca, 'YScale', 'log')
     set(gca, 'XScale', 'log')
     set(gca, 'FontSize',8)
     caxis([c_min c_max])
     ylim([1 100])
     a=get(gca,'XLim');
     x=max(a)-(max(a)-1.25*min(a));
     b=get(gca,'YLim');
     y=max(b)-(max(b)-min(b))/1.5;
     text(x,y,['\it r = ',num2str(correlations_reff_beta_agg{i}(2),'%.3f'),...
      '\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.4f'),...
      '\newline','\itACI_{r} = ',num2str(agg_ACI_corr.ACI_reff{i},'%.3f'),...
      '\newline','\itp-value = ',num2str(mdl_lidar_reff{i}.Coefficients.pValue(2),'%.4f')], ...
      'FontSize',10);
     agg_ACI_corr.reff_r{i} = correlations_reff_beta_agg{i}(2);
     agg_ACI_corr.reff_r2{i} = mdl_lidar_reff{i}.Rsquared.Ordinary;
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
    
for ii=1:n_bins_agg
    fprintf('%8.0f < LWP < %8.0f & %8.2f & %8.2f & %8.2f & %8.0f \\\\ \n', ...
        round(min_lwp + (bin_size*(ii-1))), ...
        round(min_lwp + (bin_size*ii)), ...
        (poly_reff_aci_agg{ii}(1)*-1), ...
        correlations_reff_beta_agg{ii}(2), ...    
        mdl_lidar_reff{ii}.Rsquared.Ordinary, ...
        numel(data_lidar_corr{ii}))
    n_data(ii) = numel(data_lidar_corr{ii});
end 
data_points = sum(n_data)
   
% % 
% fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_r_' ...
%      num2str(round(bin_size))]);
% export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
% clear TitleFigure xlim ylim x y a b fig_name fig h2

%% Plot the scatter plot of ATB and N_d - ACI_N is not dependent on LWP, 
%  just one plot
  
[correlations_Nd_beta, pvalue_Nd_beta]= ...
      corrcoef(log(data_lidar_level_agg(id_cloud_agg)),log(data_Nd(id_cloud_agg)));
[poly_Nd_beta, stuct_Nd_beta, mu_Nd_beta] = ...
      polyfit(data_lidar_level_agg(id_cloud_agg), data_Nd(id_cloud_agg),1);
[poly_Nd_aci, stuct_Nd_aci] = ...
      polyfit(log(data_lidar_level_agg(id_cloud_agg)), log(data_Nd(id_cloud_agg)),1);

TitleFigure=['Cabauw ACCEPT aggregated correlations beta vs N_d'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 8.3 8.3]);
    scatter(data_lidar_level_agg(id_cloud_agg), data_Nd(id_cloud_agg), s, lwp_cloud_agg(id_cloud_agg) , 'fill')
    hold on
    mdl_lidar_Nd = fitlm(log(data_lidar_level_agg(id_cloud_agg)),log(data_Nd(id_cloud_agg)));
    legend('off')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
    set(gca,'FontSize',8)
%     xlim([0.1*1.e-3 3*1.e-3])
%     ylim([100 2500])
    a=get(gca,'XLim');
    x=max(a)-(max(a)-100* min(a));
    b=get(gca,'YLim');
    y=max(b)-(max(b)-min(b))/1.0005;
    text(x,y,[ '\it r = ',num2str(correlations_Nd_beta(2),'%.2f'), ...
        '\newline','\itr^2 = ',num2str(mdl_lidar_Nd.Rsquared.Ordinary,'%.4f'), ...
        '\newline','\itACI_{N} = ',num2str(poly_Nd_aci(1),'%.2f')],'FontSize',10);
    caxis([c_min c_max])
    title([num2str(round(min_lwp)) ' < LWP < ', num2str(round(max_lwp))],'FontSize',10, ...
         'FontWeight','normal')
    xlabel('\beta [sr^{-1}')
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
    agg_ACI_corr.ACI_N = poly_Nd_aci(1);
    agg_ACI_corr.N_r = correlations_Nd_beta(2);
    agg_ACI_corr.N_r2 = mdl_lidar_Nd.Rsquared.Ordinary;
    
    
%% Histograms of data considered
str = {['# days=' num2str(numel(days))],['# profiles=' num2str(data_points)]};
str_p_r = {['Effective Radius [\mum]']};
%histogram, mean and standard deviation of effective radius    
  figure
  hist_ok(data_reff_level_agg,50)
  title('Effective Radius')
  ylabel('counts')
  xlabel(str_p_r)
  grid on
  annotation('textbox',[0.6,0.45,0.1,0.1],...
                      'String',str,...
                      'Margin',4,...                
                      'edgecolor','w',...
                      'backgroundcolor','w',...
                      'LineWidth',1,...
                      'FitBoxToText','on')    
%histogram, mean and standard deviation of the attenuated backscatter    
str_p_atb = {['Integrated Attenuated Backscatter [sr^{-1}]']};
  figure
  hist_ok(data_lidar_level_agg,50)
  title('Integrated Attenuated Backscatter')
  ylabel('counts')
  xlabel(str_p_atb)
  grid on
  annotation('textbox',[0.6,0.45,0.1,0.1],...
                      'String',str,...
                      'Margin',4,...                
                      'edgecolor','w',...
                      'backgroundcolor','w',...
                      'LineWidth',1,...
                      'FitBoxToText','on')        
                  
%histogram, mean and standard deviation of the CDNC   
str_p_N = {['Droplet concentration [#/cm^3]']};
  figure
  hist_ok(data_Nd(id_cloud_agg),50)
  title('Cloud Droplet Number Concentration')
  ylabel('counts')
  xlabel(str_p_N)
  grid on
  annotation('textbox',[0.6,0.45,0.1,0.1],...
                      'String',str,...
                      'Margin',4,...                
                      'edgecolor','w',...
                      'backgroundcolor','w',...
                      'LineWidth',1,...
                      'FitBoxToText','on')     
%histogram, mean and standard deviation of LWP   
str_p_N = {['LWP [g/cm^2]']};
TitleFigure=('Liquid Water Path');
  figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 15 15]);
  hist_ok(lwp_cloud_agg(id_cloud_agg),50)
  title('Liquid Water Path')
  ylabel('counts')
  xlabel(str_p_N)
  grid on
  annotation('textbox',[0.6,0.45,0.1,0.1],...
                      'String',str,...
                      'Margin',4,...                
                      'edgecolor','w',...
                      'backgroundcolor','w',...
                      'LineWidth',1,...
                      'FitBoxToText','on') 
                  
%% Histogram, meteo data
TitleFigure=('Meteorological data - Pressure');
  figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 15 15]);
% subplot(1,3,1)
hist_ok(data_pressure_interval./100,10)
  title('Pressure')
  ylabel('counts')
  xlabel('hPa')
  grid on

TitleFigure=('Meteorological data - Temperature');
  figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 15 15]);  
% subplot(1,3,2)
hist_ok(data_temperature_interval,10)
  title('Temperature')
  ylabel('counts')
  xlabel('K')
  grid on

TitleFigure=('Meteorological data - specific humidity');
  figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 15 15]);  
% subplot(1,3,3)
hist_ok(data_humidity_interval,10)
  title('Specific Humidity')
  ylabel('counts')
  xlabel(' ')
  grid on

                  
                  
                  
                  
                  
                  
%% Scatter ACI against correlations anr r^2
ACI = cell2mat(agg_ACI_corr.ACI_reff);
reff_r = cell2mat(agg_ACI_corr.reff_r);
reff_r2 = cell2mat(agg_ACI_corr.reff_r2);
lwp_s = min_lwp:10:max_lwp;



TitleFigure=('ACI vs. Correlation Coefficient');
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 15 15]);
scatter(ACI, reff_r, 150, lwp_s(1:end-1) , 'fill')
caxis([min(lwp_s) max(lwp_s)])
title('ACI vs. Correlation Coefficient','FontSize',10, ...
         'FontWeight','normal')
xlabel('ACI')
ylabel('r ')
h2=colorbar;
xlabel(h2,'LWP [g m^{-2}]')
    set(h2,'Location','southoutside');
    ax = gca;
    axpos = ax.Position;
    cpos = h2.Position;
    cpos(4) = 0.5*cpos(4);
    h2.Position = cpos;
    ax.Position = axpos;
    
    
TitleFigure=('ACI vs. Coefficient of Determination');
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 15 15]);
scatter(ACI, reff_r2, 150, lwp_s(1:end-1) , 'fill')
caxis([min(lwp_s) max(lwp_s)])
title('ACI vs. Coefficient of Determination','FontSize',10, ...
         'FontWeight','normal')
xlabel('ACI')
ylabel('r^2 ')
h2=colorbar;
xlabel(h2,'LWP [g m^{-2}]')
    set(h2,'Location','southoutside');
    ax = gca;
    axpos = ax.Position;
    cpos = h2.Position;
    cpos(4) = 0.5*cpos(4);
    h2.Position = cpos;
    ax.Position = axpos;

                  
% % 
% fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_N']);
% export_fig(sprintf(fig_name), '-eps', '-transparent')
% savefig(fig_name)
% %% Divide by specigic humidity
% 
% %% Define subplot net
% sx2 = 1; % number of rows
% sy2 = 2; %number of columns
% bins_sh = sx*sy; % must be divided by 2
% sh_cloud = agg_cloudnet.specific_humidity(is_data_agg);
% max_sh = max(sh_cloud);
% min_sh = min(sh_cloud);
% bin_sh =  round((max_sh- min_sh)/bins_sh);
% c_min = min_lwp;
% c_max = max_lwp;
% 
% for i = 1:n_bins_agg
%   
%  id_sh = find(sh_cloud >= min_sh + bins_sh*(i-1) & ...
%          sh_cloud <= min_sh + bins_sh*i);
%  sh_id{i} = id_sh;
%  % Divide data into bins based on the LWP
%     for j = 1:length(sh_id{i})
%       corr_lidar_sh(i,j) = data_lidar_level_agg(sh_id{i}(j));
%       corr_radar_sh(i,j) = data_radar_level_agg(sh_id{i}(j));
%       corr_reff_sh(i,j) = data_reff_level_agg(sh_id{i}(j));
%     end
%      
%  % Calculate correlation coefficient for each scatter plot
%   data_lidar_corr_sh{i} = (corr_lidar_sh(i,:));
%   data_radar_corr_sh{i} = (corr_radar_sh(i,:));
%   data_reff_corr_sh{i} = (corr_reff_sh(i,:));
%   
%   
% %   data_reff_corr{i} = data_reff_corr{i}(~isnan(data_reff_corr{i}));
% %   data_lidar_corr{i} = data_lidar_corr{i}(~isnan(data_reff_corr{i}));
% %   data_radar_corr{i} = data_radar_corr{i}(~isnan(data_reff_corr{i}));
%   data_Nd_corr_sh{i} = (data_Nd(sh_id{i})) ;
% %   data_Nd_corr{i} = data_Nd_corr{i}(~isnan(data_reff_corr{i})); 
%   data_lwp_corr_sh{i} = (lwp_cloud_agg(sh_id{i})) ;
% %   data_lwp_corr{i} = data_lwp_corr{i}(~isnan(data_reff_corr{i}));
%   
%   data_time_sh{i} = (time(sh_id{i}));
% %   data_time{i} = data_time{i}(~isnan(data_reff_corr{i}));
%   % Calculate correlations
% 
%   [correlations_reff_beta_sh{i}, pvalue_reff_beta_sh{i}] = ...
%       corrcoef(log(data_lidar_corr_sh{i}),...
%       log(data_reff_corr_sh{i}));
%   [poly_reff_beta_sh{i}, stuct_reff_beta_sh{i}, mu_reff_beta_sh{i}] = ...
%       polyfit(data_lidar_corr_sh{i}, data_reff_corr_sh{i},1);
%   [poly_reff_aci_sh{i}, stuct_aci_beta_sh{i}] = ...
%       polyfit(log(data_lidar_corr_sh{i}), log(data_reff_corr_sh{i}),1);
%   
%   ACI_corr_sh.ACI_reff{i} = -1*poly_reff_aci_sh{i}(1);
% clear id_lwp corr_*
% end



% %% Plot scatter plot of ATB and Reff based on sh
% TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
%     num2str(year) '-' 'correlations beta vs Reff'];
% 
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[10 30 20 15]);
%     for i = 1:n_bins
%      cx(i) = subplot(sx,sy,i) ; 
%      scat_fig_R(i) = scatter(data_lidar_corr{i},(data_reff_corr{i}), s, ...
%          data_lwp_corr{i} , 'fill');
%      hold on
% %      [corr_lidar_reff, gof_reff] = fit((data_lidar_corr{i}'),(data_reff_corr{i}'),'poly1');
%      mdl_lidar_reff{i} = fitlm(log(data_lidar_corr{i}'),log(data_reff_corr{i}'));
% %      plot(corr_lidar_reff)
%      legend('off')
%      set(gca, 'YScale', 'log')
%      set(gca, 'XScale', 'log')
% %      set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
%      set(gca, 'FontSize',8)
%      caxis([c_min c_max])
%      xlim([0.1*1.e-3 5*1.e-3])
%      ylim([1 100])
%      a=get(gca,'XLim');
%      x=max(a)-(max(a)-1.5*min(a));
%      b=get(gca,'YLim');
%      y=max(b)-(max(b)-min(b))/1.5;
%      text(x,y,['\it r = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
%       '\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
%       '\newline','\itACI_{r} = ',num2str((poly_reff_aci{i}(1))*-1,'%.2f')], ...
%       'FontSize',10);
%      ACI_corr.(date).reff_r{i} = correlations_reff_beta{i}(2);
%      ACI_corr.(date).reff_r2{i} = mdl_lidar_reff{i}.Rsquared.Ordinary;
% %
%      title([num2str(round(min_lwp + (bin_size*(i-1)))) ' < LWP < ',  ... 
%        num2str(round(min_lwp + (bin_size*i)))],'FontSize',10,'FontWeight','normal')
%     end
%     H=labelEdgeSubPlots('\beta [sr^{-1} * 10^{-3}]','r_e [\mum]');
%     h2=colorbar;
%     set(h2, 'Position', [.925 .11 .0181 .81])
%     h2.Label.String = 'LWP \newline [g m^{-2}]';
%     h2.Label.FontSize = 8;
%     h2.Label.Rotation = 0; 
%     h2.Label.Position = [2 15 0];

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




