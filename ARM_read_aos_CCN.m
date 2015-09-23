% This script is used to read data from the Aerosol Observing System.
% 
% The data is from the ARM Mobile Facility at the Graciosa Islnd, Azores,
% Portugal. 
%
%
% Karolina Sarna
% k.sarna@tdelft.nl
%
% January 2014



[FileName4,PathName]=uigetfile('*.cdf', 'Select Aerosol Observing System Data File', ...
    '/media/data/data/ARM_Azores/ARM_AOS_chosen_days/');
%    '/media/sarna/data/ARM_data/Azores_Graciosa/AOS_M1/chosen_days/');
if isequal(FileName4,0)
   disp('User selected Cancel')
   return
end

data_file4=[PathName FileName4];
% date_string_aos = FileName4 ;
date = FileName4(17:24) ;


ncid4 = netcdf.open(data_file4, 'NOWRITE');
%ncdisp(data_file)

%% Reading required variables from data file 

%% Basic data
% base_time - Base time in Epoch - seconds since 1970-1-1 0:00:00 0:00
% time_offset - Time offset from base_time - seconds since 

varid_base_time = netcdf.inqVarID(ncid4,'base_time');
 base_time_aos = netcdf.getVar(ncid4,varid_base_time);

varid_time_offset = netcdf.inqVarID(ncid4,'time_offset');
time_offset_aos = netcdf.getVar(ncid4,varid_time_offset);

%time_aos = datestr((time_offset_aos + double(base_time_aos)) / 86400 ...
%    + datenum(1970,1,1));

time_aos = time_offset_aos /60 /60 ; %change from seconds to hours

%% Aerosol properties
% Condensation nuclei concentration number (CNCN)
% Cloud droplet concentration number from summation of 
% condensation nuclei counter size bins

varid_N_CPC_1 = netcdf.inqVarID(ncid4,'N_CPC_1');
CNCN = netcdf.getVar(ncid4,varid_N_CPC_1);

CNCN (CNCN == -9999) = NaN;

varid_N_CCN_1 = netcdf.inqVarID(ncid4,'N_CCN_1');
CDCN = netcdf.getVar(ncid4,varid_N_CCN_1);

CDCN (CDCN == -9999) = NaN;


%% Plot Cloud Condensation Nuclei and Cloud Droplet Concentration

% plot(time_aos, CNCN, 'g', time_aos, CDCN, 'm')
% title([date, ' CNCN & CDCN'])
% set(gca, 'Ydir', 'normal')
% %axis([t_start,t_end,alt_min,alt_max])
% % set(gca,'XTick',t_start:t_end/12:t_end)
% %xlim([t_start t_end])
% %set(gca,'YTick',alt_min:alt_max/4:alt_max)
% legend('CNCN','CDCN')
% xlabel('Time [UTC]')
% ylabel('')