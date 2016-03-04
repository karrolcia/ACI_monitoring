%% Read in multiple netcdf files
% input = directory with selected Cloudnet files
% Create data structure from multiple netcdf files. Saves all days
% separetely in 'data' and a combined data serires of all days of interest
% in 'ACI' (both for measurements and retrieval data (for retrieval the
% names are 'data_ret' and 'ACI_ret' respectively. 


%% Step 1. 
% Read in the Cloudnet data file that you wish to process
% cloudnet categorization file

% reading cloudnet categorization file from user selected directory
folder_name = ('');
ncfiles = dir([folder_name '/*categorize.nc']);
% 
% if isempty(ncfiles) 
% disp('No necessary Cloudnet categorization data file exists')
% return
% end

for file = ncfiles'

filename=[folder_name '/' file.name];

ncid = netcdf.open(filename,'NC_NOWRITE');
day = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'day');
month = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'month');
year = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'year');
location = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'location');

date = (strcat(num2str(year),'-',num2str(month),'-',num2str(day)));
date = datestr(date,'mmmdd_yyyy');
% create a cell array of variables to load
variables_instruments = {'time', 'height', 'lwp', 'lwp_error', 'Z', ...
    'Z_error','Z_bias','v','lidar_wavelength', 'ldr', 'beta', ...
    'beta_bias', 'beta_error', 'rainrate','category_bits', 'model_height'};
variables_model = {'temperature','pressure', ...
    'specific_humidity', 'uwind','vwind'};
% loop over the variables
 for j=1:numel(variables_instruments)
    % extract the jth variable (type = string)
    var = variables_instruments{j};
    % use dynamic field name to add this to the structure
    cloudnet.(date).(var) = ncread(filename,var)';
    % convert from single to double, if that matters to you (it does to me)
    if isa(cloudnet.(date).(var),'single')
        cloudnet.(date).(var) = double(cloudnet.(date).(var));
    end
    cloudnet.(date).(var)(cloudnet.(date).(var)==-999)=NaN;
 end
  for j=1:numel(variables_model)
    % extract the jth variable (type = string)
    var_m = variables_model{j};
    % use dynamic field name to add this to the structure
    cloudnet.(date).(var_m) = ncread(filename,var_m)';
    % convert from single to double, if that matters to you (it does to me)
    if isa(cloudnet.(date).(var_m),'single')
        cloudnet.(date).(var_m) = double(cloudnet.(date).(var_m));
    end
    % Change model height to height of the other instruments
    cloudnet.(date).(var_m) = interp1q(cloudnet.(date).model_height', ...
       cloudnet.(date).(var_m)',cloudnet.(date).height')';
    cloudnet.(date).(var_m)(cloudnet.(date).(var_m)==-999)=NaN;
  end
% Create a time variable with the datestamp
date2 = datestr(date,'yyyy-mm-dd');
serdate = cloudnet.(date).time./24;
cloudnet.(date).time2 = strcat(date2 ,'_', datestr(serdate, 'HH:MM:SS'))';
cloudnet.(date).day = day;
cloudnet.(date).month = month;
cloudnet.(date).year = year;
cloudnet.(date).location = location;

end

days = (fieldnames(cloudnet));
variables_hor = {'lwp','time2','day','month','year','location'}; % one dimensional variables
variables_vert = {'Z','v','ldr', 'beta','temperature','pressure', ...
    'specific_humidity', 'uwind','vwind','category_bits' };
other = {'height', 'model_height'}; %two dimensional variables (and time2)

for jj=1:numel(variables_vert)
    var = variables_vert{jj};
    agg_cloudnet.(var) = [];
  for i = 1:length(days)
    name = char(days(i));
    agg_cloudnet.(var) = vertcat(agg_cloudnet.(var), cloudnet.(name).(var));
  end
end
for jj=1:numel(variables_hor)
    var = variables_hor{jj};
    agg_cloudnet.(var) = [];
  for i = 1:length(days)
    name = char(days(i));
    agg_cloudnet.(var) = horzcat(agg_cloudnet.(var), cloudnet.(name).(var));
  end
end
agg_cloudnet.(other{1}) = cloudnet.(name).(other{1});
agg_cloudnet.(other{2}) = cloudnet.(name).(other{2});

%% Step 2. 
% Read in the water cloud properties retrieval data
% cloudnet categorization file

% reading cloudnet categorization file from user selected directory
folder_name2 = ('/media/data/data/Cloudnet_Azores/ACI@Azores_ret/');
ncfiles2 = dir([folder_name2 '/wcp_*']);

% if isempty(ncfiles2) 
% disp('No necessary retrievalsn data file exists')
% return
% end

for i = 1:size(ncfiles2,1)

file = ncfiles2(i).name;
filename=[folder_name2 '/' file];
date2 = char(days(i));
variables_retrieval = {'lwc_hm','re_hm','N_hm','ext_hm','tau_hm','cb_layer',...
    'ct_layer','Z_cn','lwc_ad'};
% loop over the variables
 for j=1:numel(variables_retrieval)
    % extract the jth variable (type = string)
    var_r = variables_retrieval{j};
    % use dynamic field name to add this to the structure
    retrieval_knist.(date2).(var_r) = ncread(filename,var_r)';
    % convert from single to double, if that matters to you (it does to me)
    if isa(retrieval_knist.(date2).(var_r),'single')
        retrieval_knist.(date2).(var_r) = double(retrieval_knist.(date2).(var_r));
    end 
 end
end
variables_hor = {'N_hm', 'tau_hm','cb_layer','ct_layer'}; % one dimensional variables
variables_ver = {'lwc_hm','re_hm','ext_hm','Z_cn','lwc_ad'}; %two dimensional variables

for jj=1:numel(variables_ver)
    var = variables_ver{jj};
    agg_retrival.(var) = [];
  for i = 1:length(days)
    name = char(days(i));
    agg_retrival.(var) = vertcat(agg_retrival.(var), retrieval_knist.(name).(var));
  end
end
for jj=1:numel(variables_hor)
    var = variables_hor{jj};
    agg_retrival.(var) = [];
  for i = 1:length(days)
    name = char(days(i));
    agg_retrival.(var) = horzcat(agg_retrival.(var), retrieval_knist.(name).(var));
  end
end


clearvars('-except', 'agg_cloudnet', 'agg_retrival', 'retrieval_knist', 'cloudnet') 
