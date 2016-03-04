%% Processing data daily data
close all; clear all;
ACI_monitoring_netcdf_read
% agg_cloudnet - aggregated data from Cloudnet for the chosen days
% agg_retrieval - aggregated retrieved data (based on Cloudnet) for the
% chosen days
% cloudent - daily data from Cloudnet for the chosen days
% retrieval_knist - daily retrieved data for the chosen days

%% Select data
%% Make a loop over the fields in cloudnet.(date)
days = (fieldnames(cloudnet));

for id_f = 1: numel(days)
date =  char(days(id_f));  

day = cloudnet.(date).day;
month = cloudnet.(date).month ;
year = cloudnet.(date).year ;
location = cloudnet.(date).location ;

time = cloudnet.(date).time;
height = cloudnet.(date).height;


%% Array dimensons
n_time = size(cloudnet.(date).time,2);
n_height = size(cloudnet.(date).height,2);
height_res = (cloudnet.(date).height(2)-cloudnet.(date).height(1));

%% Step 3.
% Find peak of the Attenuated backscatter - if peak is a below 6 range exclude
% the profile - this distance is necessary to be able to comapre aerosols
% below the cloud
[Peak_P, I_max] = max(cloudnet.(date).beta, [], 2) ;
I_max(isnan(I_max))=1;
low_cloud = (find(height > 2000, 1, 'first'));
filter_beta = find(I_max(I_max ==1));
var_to_nan = {'beta', 'Z', 'v', 'lwp'};
var_to_nan2 ={'re_hm', 'N_hm'};

for ih = 1:n_time
    if I_max(ih) < 16
        for jj=1:numel(var_to_nan)
            var = var_to_nan{jj};
            cloudnet.(date).(var)(ih)=NaN;
        end
        for ji = 1:numel(var_to_nan2)
            var2 = var_to_nan2{ji};
            retrieval_knist.(date).(var2)(ih)=NaN;
        end
    end
end

I_max(I_max < 16)=NaN;
I_max(I_max > low_cloud)=NaN;
% Filtering the data based on the Cloudnet catgory bits
bits=cloudnet.(date).category_bits;
% Bit 0: Small liquid droplets are present.
bit0 = bitget(bits,1);
% Bit 1: Falling hydrometeors are present - drizzle and/or rain; 
bit1 = bitget(bits,2);
% Bit 3: Melting ice particles are present.
bit3 = bitget(bits,4);
% Bit 4: Aerosol particles are present and visible to the lidar.
bit4 = bitget(bits,5);
%Insects are present and visible to the radar.
bit5=bitget(bits,6);
 
%% Detection of the data that can be used for processing: single layer,
% liquid clouds.
% Processing restriction in regard of LWP and Z values 
% (same as for the microretrievla - for consistency of the compared results)
z_max=-15; %dBZ threshold for drizzle (see e.g. Frisch et al. 2000, Matrosov, 2004)
z_min=-60;
lwp_min=20;%g/m^2 threshold for MWR-LWP retrieval
lwp_max=150;
water_data=bit0;
% only use cases when LWP is larger than 15g/m^2
water_data(cloudnet.(date).lwp<lwp_min,:)=0;
water_data(cloudnet.(date).lwp>lwp_max,:)=0;
% selection criteria: no drizzle and single layer

%% Declare parameters used for filtering 
Z_water=NaN(size(cloudnet.(date).Z(:,:)));    %single layer water cloud Z-profile for retrieval input

n_liq=NaN(n_time,1);      %number of detected liquid water cloud bins per profile
n_z=NaN(n_time,1);        %number of detecetd Z-bins per profile 
% 
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
%    if ~isempty(il)
% %   cb_drizzle(i)=height(il(1));
%      dl=find(bit1(i,il(1))==1, 1); 
%     if length(dl)<2 
%         bit1(i,:)=0;
%     end
%     end
     
 % Find cases without drizzle/ice
 % (even when drizzle/ice is detected below cloud base, 
 % the cloud layer Z-profile is not considered)
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
     Z_water(i,is_water)=cloudnet.(date).Z(i,is_water);
     
     % Apply the threshold for the Z max/min
     is_maxmin=find(Z_water(i,:)>z_max, 1);
     if ~isempty(is_maxmin)
     Z_water(i,:)=NaN;
     end
     clear is_maxmin
     is_maxmin=find(Z_water(i,:)<z_min, 1);
     Z_water(i,is_maxmin)=NaN;
     
     % only use cases where Z-profile within cloud layer is continous
     is_cont=find(~isnan(Z_water(i,:)));
     if ~isempty(is_cont);
         double_layer=find(diff(is_cont)>1, 1);
       if ~isempty(double_layer)
           Z_water(i,:)=NaN;
       end
     end
 end    
clear ic d im

 
%%  Convert dbZ to Z
Z_linear = 10.^((Z_water-180)./10);

Z_linear_temp = Z_linear ;
Z_linear_temp(isnan(Z_linear_temp))=0;
temp = Z_linear_temp .* height_res;
int_Z = (cumsum(temp,2));
int_Z(isnan(Z_linear))=NaN;
is_refl = find(nansum(int_Z,2) ~= 0);

%% Integrate attenuated_backscatter
beta_sel = NaN(size(cloudnet.(date).beta(:,:)));
beta_sel(is_refl,:) = cloudnet.(date).beta(is_refl,:);
beta_sel(:,1:4)=NaN; % First four - overlap
temp_beta = beta_sel;
temp_beta(isnan(temp_beta)) = 0;
int_beta_atten = cumsum(temp_beta,2).*height_res;
int_beta_atten(int_beta_atten<=0) = NaN;
% int_beta_atten(int_beta_atten< 1.0e-5) = NaN;
% int_beta_atten(int_beta_atten>0.0025) = NaN;

clear temp_beta
Azores_plots_daily
ACI_reff =  ACI_corr.(date).ACI_reff;
reff_corr = ACI_corr.(date).reff_r;
reff_cr2 = ACI_corr.(date).reff_r2;
end