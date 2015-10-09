
% ACI Monitoring Scheme Data Processing program
%   
% This program filters data in accordance with the data selesction criteria
% of the ACI Monitoring Scheme. The criteria are as follows:
%
% 1. ACI monitoring scheme applies only to liquid water clouds on top of 
%    the boundary layer in well-mixed conditions. This is doneto ensures 
%    that the cloud is not decoupled from the boundary layer and the 
%    aerosol background below the cloud (Feingold et al., 2006).
% 2. Any form of precipitation, including drizzle, needs to be eliminated 
%    from the profile as it can obscure the formative stage of a cloud 
%    (Feingold et al., 2003). We use the Cloudnet categorisation data for
%    the classification of the observed targets. 
% 3. The assumption of ACI is that the variation in the aerosol concentration
%    affects the variation in the cloud properties. Thus, both aerosol 
%    and cloud parameters need to vary to observe ACI. 
% 4. The ACI monitoring scheme relies on measurements from three separate 
%    instruments. Only profiles where all three instruments provide good 
%    quality data are analysed.
% 5. We apply a constraint on LWP. In most cases we divide the data into 
%    bins of LWP of 10 gm-2.  LWP % should be above 30 gm-2 and below 
%    150 gm-2. Values below 30 gm-2 are disregarded because of the 
%    uncertainty of LWP calculated from MWR, which is around 15 gm-2 
%    (Turner et al., 2007). The values above 150 gm-2 170 are excluded to 
%    avoid precipitating clouds.
%
% After applying filtering data is further prepared for the analysis. This
% preparation includes converting Radar Reflectivity Factor to a linear
% values, integrating the Attenuated Backscatter Coefficient and Linera
% Radar Reflectivity Factor. If the depolarisation data is available, the
% Attenuated Backscatter Coefficient is divided between the parallel and
% perpendicular channel to detect the cloud base.
%
% Additional part of this program is the ACI_scheme_plots.m script which
% creates time series and scatter plots and performs the statistical
% analysis of the data.
% 
% version 1.2
% *AUTHOR*
% Karolina Sarna
% *DATE*
% 2015-07-16
% 
% k.sarna@tudelft.nl
%
 
clear all;
%% For the purpose of testing the file choosing is set to a one specific file - Cabauw 20141002

% function ACI_scheme(time_start,time_end,max_haight)
%% Set the start time and end time (default 0 to 24)
% if nargin<3,
    time_start = 0;
    time_end = 24;
    max_height = 3000;
% end;  


%% Step 1. 
% Read in the Cloudnet data file that you wish to process
% cloudnet categorization file

[FileName,PathName]=uigetfile('*.nc', 'Select Cloudnet categorization Data File', ...
    '/media/data/data/Cloudnet_Azores');
if isequal(FileName,0)
   disp('User selected Cancel')
   return
end

File=[PathName FileName];
% File='/media/data/data/ACCEPT/categorize/20141002_cesar_categorize.nc';
ncid = netcdf.open(File,'NC_NOWRITE');

% Attributes of the file
day = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'day');
month = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'month');
year = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'year');
location = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'location');

% % Array dimensions - time and height
% [dimname, n_time] = netcdf.inqDim(ncid,0);
% [dimname1, n_height] = netcdf.inqDim(ncid,1);
% clear dimname dimname1

id_time=netcdf.inqVarID(ncid,'time');
id_height=netcdf.inqVarID(ncid,'height');
id_model_height=netcdf.inqVarID(ncid,'model_height');
id_lwp=netcdf.inqVarID(ncid,'lwp');
id_lwp_error=netcdf.inqVarID(ncid,'lwp_error');
id_Z=netcdf.inqVarID(ncid,'Z');
id_Z_error=netcdf.inqVarID(ncid,'Z_error');
id_Z_bias=netcdf.inqVarID(ncid,'Z_bias');
id_doppler_v = netcdf.inqVarID(ncid,'v');
id_lidar_wavelength = netcdf.inqVarID(ncid,'lidar_wavelength');
id_ldr = netcdf.inqVarID(ncid,'ldr');
id_beta = netcdf.inqVarID(ncid,'beta');
id_beta_bias = netcdf.inqVarID(ncid,'beta_bias');
id_beta_error = netcdf.inqVarID(ncid,'beta_error');
id_temp=netcdf.inqVarID(ncid,'temperature');
id_p=netcdf.inqVarID(ncid,'pressure');
id_specific_humidity = netcdf.inqVarID(ncid,'specific_humidity');
id_specific_uwind = netcdf.inqVarID(ncid,'uwind');
id_specific_vwind = netcdf.inqVarID(ncid,'vwind');
id_rainrate = netcdf.inqVarID(ncid,'rainrate');
id_cat=netcdf.inqVarID(ncid,'category_bits');

% Variables needed for the data processing. Selection of the requested time
% period

time=netcdf.getVar(ncid,id_time);
time = time((time>=time_start)&(time<=time_end));

height=netcdf.getVar(ncid,id_height);
height = height((height <= max_height));
dh=height(2)-height(1);
model_height=netcdf.getVar(ncid,id_model_height);
model_height=model_height((model_height <= max_height));
lwp=netcdf.getVar(ncid,id_lwp);
lwp(lwp==-999)=NaN;
lwp = lwp((time>=time_start)&(time<=time_end));
Z=netcdf.getVar(ncid,id_Z)';
Z(Z==-999)=NaN;
Z(Z<=-60)=NaN;
Z = Z((time>=time_start)&(time<=time_end),(height <= max_height));
Z_error=netcdf.getVar(ncid,id_Z_error)';
Z_error(Z_error==-999)=NaN;
Z_error = Z_error((time>=time_start)&(time<=time_end),(height <= max_height));
Z_bias=netcdf.getVar(ncid,id_Z_bias)';
doppler_v = netcdf.getVar(ncid,id_doppler_v)';
doppler_v(doppler_v==-999)=NaN;
doppler_v = doppler_v((time>=time_start)&(time<=time_end),(height <= max_height));
lidar_wavelength = netcdf.getVar(ncid,id_lidar_wavelength);
lin_depol_ratio = netcdf.getVar(ncid,id_ldr)';
lin_depol_ratio(lin_depol_ratio==-999)=NaN;
lin_depol_ratio = lin_depol_ratio((time>=time_start)&(time<=time_end),(height <= max_height));
beta_atten = netcdf.getVar(ncid,id_beta)';
beta_atten(beta_atten==-999)=NaN;
beta_atten(beta_atten<=0)=NaN;
beta_atten(beta_atten<=5.e-7)=NaN;
beta_atten = beta_atten((time>=time_start)&(time<=time_end),(height <= max_height));
beta_atten_bias = netcdf.getVar(ncid,id_beta_bias);
beta_atten_error = netcdf.getVar(ncid,id_beta_error)';
temperature=netcdf.getVar(ncid,id_temp)';
temperature = temperature((time>=time_start)&(time<=time_end),(model_height <= max_height));
pressure=netcdf.getVar(ncid,id_p)';
pressure = pressure((time>=time_start)&(time<=time_end),(model_height <= max_height));
specific_humidity = netcdf.getVar(ncid,id_specific_humidity)';
uwind = netcdf.getVar(ncid,id_specific_uwind)';
uwind = uwind((time>=time_start)&(time<=time_end),(model_height <= max_height));
vwind = netcdf.getVar(ncid,id_specific_vwind)';
vwind = vwind((time>=time_start)&(time<=time_end),(model_height <= max_height));
rainrate = netcdf.getVar(ncid,id_rainrate);
rainrate = rainrate((time>=time_start)&(time<=time_end));
category_bits=netcdf.getVar(ncid,id_cat)';
category_bits = category_bits((time>=time_start)&(time<=time_end),(height <= max_height));
clear id_*

%% Step 2.
% Read data from the microphysical properties retrieval - C.Knist
%CKnist_micro_retrieval_data
[FileName2,PathName2]=uigetfile('*.nc', 'Select Water Cloud Properties Data File', ...
    '/media/data/data/Cloudnet_Azores/retrievals/');
if isequal(FileName2,0)
   disp('User selected Cancel')
   return
end
data_file2=[PathName2 FileName2];
% data_file2='/media/data/data/ACCEPT/retrievals/wcp_20141002_cesar_categorize.nc';
ncid2 = netcdf.open(data_file2, 'NOWRITE');

% Read in the variables
varid_time_cn = netcdf.inqVarID(ncid2,'time_cn');
varid_height_cn = netcdf.inqVarID(ncid2,'height_cn');
varid_lwp_cn = netcdf.inqVarID(ncid2, 'lwp_cn');
varid_lwp_ad = netcdf.inqVarID(ncid2, 'lwp_ad');
varid_cb_radar = netcdf.inqVarID(ncid2, 'cb_radar');
varid_ct_radar = netcdf.inqVarID(ncid2, 'ct_radar');

time_cn = netcdf.getVar(ncid2,varid_time_cn);
time_cn = time_cn((time_cn>=time_start)&(time_cn<=time_end));
height_cn = netcdf.getVar(ncid2,varid_height_cn);
height_cn = height_cn((height <= max_height));
lwp_cn = netcdf.getVar(ncid2, varid_lwp_cn);
lwp_cn = lwp_cn((time_cn>=time_start)&(time_cn<=time_end));
lwp_ad = netcdf.getVar(ncid2, varid_lwp_ad);
lwp_ad = lwp_ad((time_cn>=time_start)&(time_cn<=time_end));
cb_radar = netcdf.getVar(ncid2, varid_cb_radar);
cb_radar = cb_radar((time_cn>=time_start)&(time_cn<=time_end));
ct_radar = netcdf.getVar(ncid2, varid_ct_radar);
ct_radar = ct_radar((time_cn>=time_start)&(time_cn<=time_end));

% Cloud boundaries from the liquid water detected by the Cloudnet
% categorisation
varid_cb_layer = netcdf.inqVarID(ncid2, 'cb_layer');
varid_ct_layer = netcdf.inqVarID(ncid2, 'ct_layer');
varid_Z_cn = netcdf.inqVarID(ncid2, 'Z_cn');
varid_lwc_ad = netcdf.inqVarID(ncid2, 'lwc_ad');

cb_layer = netcdf.getVar(ncid2, varid_cb_layer);
cb_layer = cb_layer((time_cn>=time_start)&(time_cn<=time_end));
ct_layer = netcdf.getVar(ncid2, varid_ct_layer);
ct_layer = ct_layer((time_cn>=time_start)&(time_cn<=time_end));
Z_cn = netcdf.getVar(ncid2, varid_Z_cn)';
Z_cn = Z_cn((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
lwc_ad = netcdf.getVar(ncid2,varid_lwc_ad);
lwc_ad = lwc_ad((time_cn>=time_start)&(time_cn<=time_end));

% Results from HM - radar reflectivity-homogenous mixing retrieval (hm)
varid_lwc_hm = netcdf.inqVarID(ncid2, 'lwc_hm');
varid_re_hm = netcdf.inqVarID(ncid2, 're_hm');
varid_N_hm = netcdf.inqVarID(ncid2, 'N_hm');
varid_ext_hm = netcdf.inqVarID(ncid2, 'ext_hm');
varid_tau_hm = netcdf.inqVarID(ncid2, 'tau_hm');

lwc_hm = netcdf.getVar(ncid2,varid_lwc_hm)';
lwc_hm = lwc_hm((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
re_hm = netcdf.getVar(ncid2, varid_re_hm)';
re_hm = re_hm((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
N_hm = netcdf.getVar(ncid2, varid_N_hm);
N_hm = N_hm((time_cn>=time_start)&(time_cn<=time_end));
ext_hm = netcdf.getVar(ncid2, varid_ext_hm)';
ext_hm = ext_hm((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
tau_hm = netcdf.getVar(ncid2, varid_tau_hm);
tau_hm = tau_hm((time_cn>=time_start)&(time_cn<=time_end));

% Results from Mace retrieval
varid_lwc_mace = netcdf.inqVarID(ncid2, 'lwc_mace');
varid_re_mace = netcdf.inqVarID(ncid2, 're_mace');
varid_N_mace = netcdf.inqVarID(ncid2, 'N_mace');
varid_ext_mace = netcdf.inqVarID(ncid2, 'ext_mace');
varid_tau_mace = netcdf.inqVarID(ncid2, 'tau_mace');

lwc_mace = netcdf.getVar(ncid2,varid_lwc_mace)';
lwc_mace = lwc_mace((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
re_mace = netcdf.getVar(ncid2, varid_re_mace)';
re_mace = re_mace((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
N_mace = netcdf.getVar(ncid2, varid_N_mace);
N_mace = N_mace((time_cn>=time_start)&(time_cn<=time_end));
ext_mace = netcdf.getVar(ncid2, varid_ext_mace)';
ext_mace = ext_mace((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
tau_mace = netcdf.getVar(ncid2, varid_tau_mace);
tau_mace = tau_mace((time_cn>=time_start)&(time_cn<=time_end));

% Results from Scaled adiabatic stratified retrieval (sas)
varid_lwc_sas = netcdf.inqVarID(ncid2, 'lwc_sas');
varid_re_sas = netcdf.inqVarID(ncid2, 're_sas');
varid_N_sas = netcdf.inqVarID(ncid2, 'N_sas');
varid_ext_sas = netcdf.inqVarID(ncid2, 'ext_sas');
varid_tau_sas = netcdf.inqVarID(ncid2, 'tau_sas');

lwc_sas = netcdf.getVar(ncid2,varid_lwc_sas)';
lwc_sas = lwc_sas((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
re_sas = netcdf.getVar(ncid2, varid_re_sas)';
re_sas = re_sas((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
N_sas = netcdf.getVar(ncid2, varid_N_sas);
N_sas = N_sas((time_cn>=time_start)&(time_cn<=time_end));
ext_sas = netcdf.getVar(ncid2, varid_ext_sas)';
ext_sas = ext_sas((time_cn>=time_start)&(time_cn<=time_end),(height <= max_height));
tau_sas = netcdf.getVar(ncid2, varid_tau_sas);
tau_sas = tau_sas((time_cn>=time_start)&(time_cn<=time_end));

% Results from Verticaly Uniform retrieval (sas)
varid_lwc_vu = netcdf.inqVarID(ncid2, 'lwc_vu');
varid_re_vu = netcdf.inqVarID(ncid2, 're_vu');
varid_N_vu = netcdf.inqVarID(ncid2, 'N_vu');
varid_tau_vu = netcdf.inqVarID(ncid2, 'tau_vu');

lwc_vu = netcdf.getVar(ncid2,varid_lwc_vu);
lwc_vu = lwc_vu((time_cn>=time_start)&(time_cn<=time_end));
re_vu = netcdf.getVar(ncid2, varid_re_vu);
re_vu = re_vu((time_cn>=time_start)&(time_cn<=time_end));
N_vu = netcdf.getVar(ncid2, varid_N_vu);

N_vu = N_vu((time_cn>=time_start)&(time_cn<=time_end));
tau_vu = netcdf.getVar(ncid2, varid_tau_vu);
tau_vu = tau_vu((time_cn>=time_start)&(time_cn<=time_end));

clear varid_*
netcdf.close(ncid2)
netcdf.close(ncid)


%% Array dimensons
n_time = length(time);
n_height = length(height);

%% Step 3.
% Filtering the data based on the Cloudnet catgory bits
bits=double(category_bits);
% Bit 0: Small liquid droplets are present.
bit0 = bitget(bits,1);
% Bit 1: Falling hydrometeors are present - drizzle and/or rain; 
bit1 = bitget(bits,2);
% Bit 3: Melting ice particles are present.
bit3 = bitget(bits,4);
% Bit 4: Aerosol particles are present and visible to the lidar.
% bit4 = bitget(bits,5);
%Insects are present and visible to the radar.
%bit5=bitget(bits,6);

% % Removes the bit classification above 3500m
% bit1(:,109:end)=0;
% bit4(:,109:end)=0;

% Processing restriction in regard of LWP and Z values 
%(same as for the microretrievla - for consistency of the compared results)
z_max=-15;                %dBZ threshold for drizzle (see e.g. Frisch et al. 2000, Matrosov, 2004)
z_min=-60;
lwp_min=15;               %g/m^2 threshold for MWR-LWP retrieval
lwp_max=700;

% Detection of the data that can be used for processing: single layer,
% liquid clouds.
water_data=bit0;
% only use cases when LWP is larger than 15g/m^2
water_data(lwp<lwp_min,:)=0;
water_data(lwp>lwp_max,:)=0;
% selection criteria: no drizzle and single layer

% Declare parameters used for filtering
Z_water=NaN(size(Z(:,:)));    %single layer water cloud Z-profile for retrieval input
% beta_atten_sel = NaN(size(beta_atten));
n_liq=NaN(n_time,1);      %number of detected liquid water cloud bins per profile
n_z=NaN(n_time,1);        %number of detecetd Z-bins per profile 

% This is based on the selesction from the microphysical properties
% retrieval
for i=1:n_time;
    
    %Find only use single layer clouds
     is_single=find(water_data(i,:)==1);
     if ~isempty(is_single);
       %no double layer clouds
         double_layer=find(diff(is_single)>1, 1);
       if ~isempty(double_layer)
            water_data(i,:)=0;
       end
     end
         
     il=find(bit0(i,:)==1, 1);
      if ~isempty(il)
%        cb_drizzle(i)=height(il(1));
        dl=find(bit1(i,il(1))==1, 1); 
      if length(dl)<2 
        bit1(i,:)=0;
      end
      end
     
    % Find cases without drizzle/ice
    %(even when drizzle/ice is detected below cloud base, the cloud layer Z-profile is not considered)
    is_drizzle=find(bit1(i,:)==1, 1);
      if ~isempty(is_drizzle)
     water_data(i,:)=0;
      end
    is_ice=find(bit3(i,:)==1, 1);
      if ~isempty(is_ice)
     water_data(i,:)=0;
      end
    % Length of the liquid layer detected
    n_liq(i)=length(find(water_data(i,:)==1));
end
clear id ice is d io

% single layer liquid water cloud Z-profile
% further selection criteria based on observed Z profile
for i=1:n_time;
     is_water=find(water_data(i,:)==1);
     Z_water(i,is_water)=Z(i,is_water);
     
     % Apply the threshold for the Z max/min
     is_maxmin=find(Z_water(i,:)>z_max, 1);
     if ~isempty(is_maxmin)
     Z_water(i,:)=NaN;
     end
     clear is_maxmin
     is_maxmin=find(Z_water(i,:)<z_min, 1);
     Z_water(i,is_maxmin)=NaN;
     
     %only use cases where Z-profile within cloud layer is continous
     is_cont=find(~isnan(Z_water(i,:)));
     if ~isempty(is_cont);
         double_layer=find(diff(is_cont)>1, 1);
       if ~isempty(double_layer)
           Z_water(i,:)=NaN;
       end
     end
    n_z(i)=ct_radar(i)-cb_radar(i);
end    
clear ic d im

is_refl = find(nansum(Z_water,2) ~= 0);
beta_atten_sel = NaN(size(beta_atten(:,:)));
beta_atten_sel(is_refl,:) = beta_atten(is_refl,:);

%%  Convert dbZ to Z
Z_linear = 10.^((Z_water-180)./10);
iZ_w=(nansum(Z_linear)')*(height(2)-height(1));
Z_linear_temp = Z_linear ;
Z_linear_temp(isnan(Z_linear_temp))=0;
temp = cumsum(Z_linear_temp,2);
int_Z=temp.*(height(2)-height(1));
clear Z_linear_temp 
% int_Z = cumsum(temp,2);
% int_Z_cloud = nansum(temp,2);

%% Integrate attenuated_backscatter
temp_beta = beta_atten_sel;
temp_beta(isnan(temp_beta)) = 0;
int_beta_atten = cumsum(temp_beta,2).*(height(2)-height(1));
int_beta_atten(int_beta_atten<=0) = NaN;
clear temp_beta

%% Divide the signal into the paralel and perpendicular 
%- only is the depolarisation ratio is available
[Peak_P, I_max] = max(beta_atten, [], 2) ;

id_cb_lidar = NaN(n_time,1);
id_cb_radar = NaN(n_time,1);
id_ct_radar = NaN(n_time,1);

 for i = 1:n_time
   if ~isnan(lin_depol_ratio(i,:)) ~=0
    para = beta_atten ./ (lin_depol_ratio + 1) ;
    perp = beta_atten - para;
    [Peak_perp, perp_max] = max(perp, [], 2) ; 
    if perp_max(i) > 1
        id_cb_lidar(i) = find(perp(i,1:perp_max(i)) > (Peak_perp(i)/10), 1,'first');
    else
        id_cb_lidar(i) = NaN;  
    end
    end  
    id_cloud=find(~isnan(Z_linear(i,:)));
     if ~isempty(id_cloud);
          id_cb_radar(i)=(id_cloud(1));
          id_ct_radar(i)=(id_cloud(end));
     end
 end

%% Step 4. 
% Plots
ACI_scheme_plots
 
%% Calculate correlation in data divided by adiabatic fractio (optional)
% adiabatic_fraction

%% Save all variables to .mat file
filename = [num2str(location) '-' num2str(day) '-' num2str(month) '-' ...
    num2str(year) '-' num2str(time_start) 'to' num2str(time_end)];
save(filename)