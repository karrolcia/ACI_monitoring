%% Read in multiple netcdf files


%% Step 1. 
% Read in the Cloudnet data file that you wish to process
% cloudnet categorization file

[FileName,PathName]=uigetfile('*.nc', 'Select Cloudnet categorization Data File', ...
    '/media/data/data/ACCEPT/categorize');
if isequal(FileName,0)
   disp('User selected Cancel')
   return
end

filename=[PathName FileName];
% File='/media/data/data/ACCEPT/categorize/20141002_cesar_categorize.nc';
ncid = netcdf.open(filename,'NC_NOWRITE');
day = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'day');
month = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'month');
year = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'year');
location = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'location');

date = FileName(1:8);
% create a cell array of variables to load
variables_to_load_Cloudnet = {'time', 'height', 'model_height', 'lwp', 'lwp_error',...
    'Z','Z_error','Z_bias','v','lidar_wavelength', 'ldr', 'beta', ...
    'beta_bias', 'beta_error', 'temperature','pressure','specific_humidity', ...
    'uwind','vwind','rainrate','category_bits'};


% loop over the variables
for j=1:numel(variables_to_load_Cloudnet)
    % extract the jth variable (type = string)
    var = variables_to_load_Cloudnet{j};

    % use dynamic field name to add this to the structure
    data_struct.(var) = ncread(filename,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(data_struct.(var),'single')
        data_struct.(var) = double(data_struct.(var));
    end
end

for k = 1:numel(data_struct.time)
(genvarname([date datestr])) = data_struct.time(k);
pause(60)
end
