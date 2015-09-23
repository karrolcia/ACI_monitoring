%% Calculations with the adiabatic fraction - run after running ACI_scheme
% adiabatic LWC/LWP where the "categorization" data has diagnosed that 
% liquid water is present (only single layer)
[lwc_ad,lwp_ad]=adiabatic_cloud_properties_cloudnet(temperature,pressure,model_height,height,cb_layer,ct_layer);
% adiabatic fraction term
fr=lwp./lwp_ad';

fr_cloud = fr(is_data);

bin_size_fr = ((max(fr_cloud) - min(fr_cloud))/n_bins);

for i = 1:n_bins
     id_fr = find(fr_cloud >= 0.1 +bin_size_fr*i & ...
         fr_cloud <= 0.1 + bin_size_fr*(i+1)); 
     
      fr_id{i} = id_fr;
      
     for j = 1:length(fr_id{i})
         corr_lidar_fr(i,j) = (data_lidar(id_fr(j),aero_level(id_fr(j))));
         corr_radar_fr(i,j) = (data_radar(id_fr(j),cloud_level(id_fr(j))));
         corr_reff_fr(i,j) = (data_reff(id_fr(j),cloud_level(id_fr(j))));
     end
    
  data_lidar_corr_fr{i} = corr_lidar_fr(i,:);
  data_radar_corr_fr{i} = corr_radar_fr(i,:);
  data_reff_corr_fr{i} =  corr_reff_fr(i,:);
  data_fr_corr{i} = fr(id_fr) ;
  
  poly_Z_beta_fr{i}= polyfit(data_lidar_corr_fr{i}, data_radar_corr_fr{i},1);
  poly_reff_beta_fr{i}= polyfit(data_lidar_corr_fr{i}, data_reff_corr_fr{i},1);
  
  corr_coef_Z_beta_fr = corrcoef(data_lidar_corr_fr{i},data_radar_corr_fr{i});
  correlations_Z_beta_fr{i} = corr_coef_Z_beta_fr;
  
  clear id_fr corr_*
end

TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'Adiabatic fraction correlations beta vs Z'];
figure('name', TitleFigure, 'NumberTitle','off');
for i = 1:n_bins
 % Plot the scatter plot
 cx(i) = subplot(sx,sy,i) ; 
 scatter(data_lidar_corr_fr{i}, data_radar_corr_fr{i}, s,  data_fr_corr{i} , 'fill')
 hold on
 corr_lidar_reff_fr = fit(data_lidar_corr_fr{i}', data_radar_corr_fr{i}','poly1');
 plot(corr_lidar_reff_fr)
 legend('off')
 hold off
 str5 = [' \rho = ' num2str(correlations_Z_beta_fr{i}(2),'% 1.2f')];
 text(0.00055, 1.e-20, str5);
 str6 = [' slope = ' num2str(poly_Z_beta_fr{i}(1),'% 4.3e')];
 text(0.00055, 2.e-20, str6);
 hold on
 lsline
 set(gca, 'YScale', 'log')
 set(gca, 'XScale', 'log')
 caxis([min(fr_cloud) max(fr_cloud)])
 xlim([0.5*1.e-3 3*1.e-3])
 ylim([1*1.e-20 1*1.e-18])
%  xlim([-9 -5])
%  ylim([-46 -43])

 title(['Adiabatic fr ' num2str(((bin_size_fr*(i-1)))) '-',  ... 
   num2str((bin_size_fr*i))  ' g/m^2'],'FontSize',12,'FontWeight','bold') 
 xlabel('Attenuated Backscatter Coefficient [1/m/sr]','FontSize',9)
 ylabel('Reflectivity ','FontSize',9)
end
h2=colorbar('FontSize',12);
set(h2, 'Position', [.93 .11 .0181 .8150])
xlabel(h2,'g/m^2')
clear TitleFigure