%% Plot scatter plot of ATB and Reff
TitleFigure=[num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' 'correlations beta vs Reff'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','normalized','outerposition',[0 0 1 1]);
for i = 1:n_bins
 cx(i) = subplot(sx,sy,i) ; 
 scat_fig_R (i) = scatter((data_lidar_corr{i}),(data_reff_corr{i}), s,  data_lwp_corr{i} , 'fill');
 hold on
 [corr_lidar_reff, gof_reff] = fit((data_lidar_corr{i}'),(data_reff_corr{i}'),'poly1');
 mdl_lidar_reff{i} = fitlm(log(data_lidar_corr{i}'),log(data_reff_corr{i}'));
 plot(corr_lidar_reff)
 legend('off')
 set(gca, 'YScale', 'log')
 set(gca, 'XScale', 'log')
 caxis([c_min c_max])
 xlim([0.5*1.e-3 3*1.e-3])
 ylim([1 100])
 set(gca, 'XTickLabel', {'', '', ''})
 set(gca, 'YTickLabel', {'', '', ''})
 a=get(gca,'XLim');
 x=max(a)-(max(a)-min(a))/1.2;
 b=get(gca,'YLim');
 y=max(b)-(max(b)-min(b))/4;
 text(x,y,[...%'\itr = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
    ...%'\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
    ...%'\newline','\itm = ',num2str(poly_reff_beta{i}(1),'%.2f') ...
    '\newline','\itACI_{r} = ',num2str((poly_reff_aci{i}(1))*-1,'%.2f')],'FontSize',20);
 title([num2str(round(min_lwp + (bin_size*(i-1)))) ' < LWP < ',  ... 
   num2str(round(min_lwp + (bin_size*i)))],'FontSize',20,'FontWeight','bold')
%  xlabel('\beta [sr^{-1}]')
%  ylabel('r_e [\mum]')
 set(gca, 'FontSize',16)

end
H=labelEdgeSubPlots('\beta [sr^{-1}]','r_e [\mum]');
h2=colorbar('FontSize',16);
set(h2, 'Position', [.93 .11 .0181 .8150])
xlabel(h2,'LWP [g m^{-2}]')

c=get(gcf,'children');
c=c(strmatch('axes',get(c,'type')));
pos1=ones(length(c),4);
for ii=1:length(c)
    pos1(i,:)=get(c(i),'position');
    if pos1(ii,1)==min(pos1(:,1));
       set(c(ii),'XTickLabel', {'1', '2', '3'})
    end
    if pos1(ii,2)==min(pos1(:,2));
       set(c(ii),'YTickLabel', {'1','10', '100'}) 
    end
end

% %Add x labels
% h.xlabels=[];
% for ii=1:length(c)
%     if pos(ii,2)==min(pos(:,2));
%         h.xlabels=[h.xlabels,get(c(ii),'xlabel')];
%         set(h.xlabels(end),'string',xl)
%     end
% end
% 
% 
% %Add y labels
% h.ylabels=[];
% for ii=1:length(c)
%     if pos(ii,1)==min(pos(:,1));
%         h.ylabels=[h.ylabels,get(c(ii),'ylabel')];
%         set(h.ylabels(end),'string',yl)
%     end
% end


% fig_name = ([num2str(year) num2str(month) num2str(day) '_corr_beta_r_' ...
%     num2str(round(bin_size))]);
% print(fig_name, '-depsc','-r0')
% savefig(fig_name)
clear TitleFigure xlim ylim x y a b fig_name fig