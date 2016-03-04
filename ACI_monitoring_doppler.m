%% Plot scatter plot of ATB and Reff divided by Doppler velocity
% Define the max and min value of Doppler velocity

dop_max = 2;
dop_min = 0;
sx_dop = 2;
sy_dop = 2;
n_bins_dop = sx_dop*sy_dop;
dop_bin_size =  ((dop_max - dop_min)/n_bins_dop);

for ii = 1:n_bins_dop
  
 id_v = find(data_doppler_level >= dop_min + dop_bin_size*(ii-1)  & ...
         data_doppler_level <= dop_min + dop_bin_size*ii );
 dop_id{ii} = id_v;
 % Divide data into bins based on the LWP
    for j = 1:length(dop_id{ii})
         corr_lidar_dop(ii,j) = (data_lidar(id_v(j),aero_level(id_v(j))));
         corr_radar_dop(ii,j) = (data_radar(id_v(j),cloud_level(id_v(j))));
         corr_reff_dop(ii,j) = (data_reff(id_v(j),cloud_level(id_v(j))));
    end
     
 % Calculate correlation coefficient for each scatter plot
  data_lidar_corr_dop{ii} = double(corr_lidar_dop(ii,:));
  data_radar_corr_dop{ii} = double(corr_radar_dop(ii,:));
  data_reff_corr_dop{ii} = double(corr_reff_dop(ii,:));
  data_dop_corr{ii} = double(data_doppler_level(dop_id{ii})) ;

  data_time{ii} = double(time(id_v));
  
  % Calculate correlations
  corr_coef_Z_beta = corrcoef(data_lidar_corr_dop{ii},data_radar_corr_dop{ii});
  correlations_Z_beta{ii} = corr_coef_Z_beta;
  corr_coef_reff_beta = corrcoef((data_lidar_corr_dop{ii}),(data_reff_corr_dop{ii}));
  correlations_reff_beta{ii} = corr_coef_reff_beta;
  
  % Calculate slope and R^2
  [poly_Z_beta{ii}, stuct_Z_beta{ii}, mu_Z_beta{ii}] = ...
      polyfit(data_lidar_corr_dop{ii}, data_radar_corr_dop{ii},1);
  [poly_reff_beta{ii}, stuct_reff_beta{ii}, mu_reff_beta{ii}] = ...
      polyfit(data_lidar_corr_dop{ii}, data_reff_corr_dop{ii},1);
  [poly_reff_aci{ii}, stuct_aci_beta{ii}] = ...
      polyfit(log(data_lidar_corr_dop{ii}), log(data_reff_corr_dop{ii}),1);
   
%   yresid_Z_beta {i} =  data_radar_corr_dop{i} - yfit_Z_beta{i};
%   yresid_reff_beta {i} =  data_reff_corr_dop{i} - yfit_reff_beta{i};
%   SSresid_Z_beta{i} = sum(yresid_Z_beta{i}.^2);
%   SSresid_reff_beta{i} = sum(yresid_reff_beta{i}.^2);
%   SStotal_Z_beta{i} = (length(data_radar_corr_dop{i})-1) * var(data_radar_corr_dop{i});
%   SStotal_reff_beta{i} = (length(data_reff_corr_dop{i})-1) * var(data_reff_corr_dop{i});
%   rsq_Z_beta{i} = 1 - SSresid_Z_beta{i}./SStotal_Z_beta{i};
%   rsq_reff_beta{i} = 1 - SSresid_reff_beta{i}./SStotal_reff_beta{i};

clear corr_* dop_id
end



% Plot scatter plots
TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'correlations beta vs Reff divided by Doppler velocity'];
figure('name', TitleFigure, 'NumberTitle','off');
for ii = 1:n_bins_dop
 cx_dop(ii) = subplot(sx_dop,sy_dop,ii) ; 
 scatter((data_lidar_corr_dop{ii}),(data_reff_corr_dop{ii}), s,  data_dop_corr{ii} , 'fill')
 hold on
 [corr_lidar_reff_dop, gof_reff] = fit((data_lidar_corr_dop{ii}'),(data_reff_corr_dop{ii}'),'poly1');
 plot(corr_lidar_reff_dop)
 legend('off')
 lsline
 set(gca, 'YScale', 'log')
 set(gca, 'XScale', 'log')
 caxis([dop_min dop_max])
 xlim([0.5*1.e-3 3*1.e-3])
 ylim([1 100])
 a=get(gca,'XLim');
 x=max(a)-(max(a)-min(a))/2;
 b=get(gca,'YLim');
 y=max(b)-(max(b)-min(b))/2;
 text(x,y,['\rho = ',num2str(correlations_reff_beta{ii}(2),'%.2f'),...
    '\newline','\itm = ',num2str(poly_reff_beta{ii}(1),'%.2e') ...
    '\newline','\itACI_{r} = ',num2str(abs(poly_reff_aci{ii}(1)),'%.2f')]);
 title(['Doppler velocity ' num2str((dop_min + (dop_bin_size*(ii-1)))) ' to ',  ... 
   num2str((dop_min + (dop_bin_size*ii)))  ' m^{-1}'], 'FontSize',10,'FontWeight','bold') 
 xlabel('ATB [1/m/sr * m]')
 ylabel('R_e_f_f [\mum]')
 set(gca, 'FontSize',8)
end
h2=colorbar('FontSize',8);
set(h2, 'Position', [.93 .11 .0181 .8150])
xlabel(h2,'g/m^2')
clear TitleFigure xlim ylim x y a b
