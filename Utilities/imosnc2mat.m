function s = imosnc2mat(ncfile,down_up)
% function s = imosnc2mat(ncfile)
% ncfile = file name
% down_up = orientation (1 = up, 0 = down)
% Takes IMOS netcdf format data and converts to s structure:
% s = 
% 
%             serial: '14168'
%               name: 'LR75            '
%              depth: [11458x1 double]
%              orien: [11458x1 double]
%               time: [11458x1 double]
%             bdepth: [11458x32 double]
%           pressure: [11458x1 double]
%          bpressure: [11458x32 double]
%           latitude: -8.8550
%          longitude: 127.1918
%              pitch: [11458x1 double]
%               roll: [11458x1 double]
%               head: [11458x1 double]
%        temperature: [11458x1 double]
%               vbat: [11458x1 double]
%                  u: [11458x32 double]
%                  w: [11458x32 double]
%                erv: [11458x32 double]
%                 qc: [1x4 struct]
%              units: {20x2 cell}
%               type: 'RDI ADCP '
%            mooring: 'Timor North'
%      planned_depth: 111
%           time_int: 3600
%          tbath_cal: 0
%             magdec: 2.3651
%            time_qc: [11458x1 double]
%     temperature_qc: [11458x1 double]
%        pressure_qc: [11458x1 double]
%           depth_qc: [11458x1 double]
%        latitude_qc: 1
%       longitude_qc: 1
%               u_qc: [11458x32 double]
%               w_qc: [11458x32 double]
%           pitch_qc: [11458x1 double]
%            roll_qc: [11458x1 double]
%            head_qc: [11458x1 double]
%          bdepth_qc: [11458x32 double]
%       bpressure_qc: [11458x32 double]
%           time_raw: [11458x1 double]
%        time_offset: 0
%           qcthresh: [1x1 struct]
%             brange: [1x32 double]
%         beam_angle: 20
%Bec Cowley, June, 2013
%updated March, 2016.

%open the file
nc =netcdf.open(ncfile,'NC_NOWRITE');

%get the size information
[ndims,nvars,ngatts,unlimdimid]= netcdf.inq(nc);

%get dimension information
dimname = {};dimlen = NaN*ones(ndims,1);
for a = 1:ndims
    [dimname{a}, dimlen(a)] = netcdf.inqDim(nc,a-1); %zero-based indexing
end

%get the variable information:
dimids = NaN*ones(nvars,ndims);
varname = {};
for a = 1:nvars
    [varname{a},xtype,dims,natts] = netcdf.inqVar(nc,a-1); %zero-based indexing
    dimids(a,1:length(dims)) = dims;
end

%get the global attributes information:
attname = {};
for a = 1:ngatts
    attname{a} = netcdf.inqAttName(nc,netcdf.getConstant('NC_GLOBAL'),a-1);
end

%is depth measured or inferred?
varid = netcdf.inqVarID(nc,'DEPTH');
try
    comm = netcdf.getAtt(nc,varid,'comment');
catch
    comm = [];
end
if ~isempty(findstr('Depth inferred',comm)) %1= measured, 0=inferred
    s.meas_inf = 0;
else
    s.meas_inf = 1;
end

%now get the data that we want and map to structure:
[varname,names] = get_s_names(varname,s);
[attname,gnames] = get_s_names(attname,s);

for a = 1:length(names)
    varid = netcdf.inqVarID(nc,varname{a});
    dat = netcdf.getVar(nc,varid);
    dat = squeeze(dat);
    %clean out the fill values
    try
        fillval = netcdf.getAtt(nc,varid,'_FillValue');
        ibad = dat == fillval;
        dat(ibad) = NaN;
    catch
    end
    %orient it the right way
    [m,n]=size(dat);
    if n > m
        dat = dat';
    end
    eval(['s.' names{a} '= dat;'])
    
end

for a = 1:length(gnames)    
    % Get value of global attribute.
    gattval = netcdf.getAtt(nc,netcdf.getConstant('NC_GLOBAL'),attname{a});
    eval(['s.' gnames{a} '= gattval;'])
end

    %put the instrument information in from the file global atts
    s.type = [];
 %trim the instrument name to get type:
ins = read_ins_info('instrument_info.csv');

% get the instrument types:
ii = strmatch(s.serial,ins.serial);
s.type = ins.type(ii,:);

%fix the time to matlab format:
s.time = s.time + datenum('1950-01-01 00:00:00');
%and get the time to middle of bin if it's there and offset the time:
try
    toff = ncreadatt(ncfile,'TIME','seconds_to_middle_of_measurement');
    s.time = s.time + toff/60/60/24;
    disp('Applied middle of bin time offset')
catch ME
end

if isfield(s,'u')
    %put u and v into complex format:
    u = complex(s.u,s.v);
    s.u = u;
    s = rmfield(s,'v');
    try
        s = rmfield(s,'v_qc');
    catch
    end
    %make up the depth matrix, remove bad data first!
    if isfield(s,'depth')
        ii = s.depth_qc > 2;
        dd = s.depth;
    elseif isfield(s,'depth_inferred')
        ii = s.depth_inferred_qc > 2;
        dd = s.depth_inferred;
    end
    dd(ii) = NaN;
    %get orientation
    if ~isempty(strmatch('HEIGHT_ABOVE_SENSOR',varname))
        if isfield(s,'brange')
            dat = s.brange;
        else
            varid = netcdf.inqVarID(nc,'HEIGHT_ABOVE_SENSOR');
            dat = netcdf.getVar(nc,varid);
        end
        s.bdepth = repmat(dd,1,length(dat')) - repmat(dat',length(dd),1);
        if all(dat < 0) %downward looking, add the bin depths to transducer depth
            s.orien = zeros(size(s.time));
        else %upward looking
            s.orien = ones(size(s.time));
        end
    end
end

%close the netcdf file
netcdf.close(nc);

end

function [varname,names] = get_s_names(varname,s)
%get the structure names that match what is in the netcdf file:
fid = fopen('imosParametersMap.txt');
txt = textscan(fid,'%s%s%*s%*s%*s%*s%*s%*s%*s%*s\n','delimiter',',');
fclose(fid);

%match the varnames and get equivalents:
names = cell(size(varname));
for a = 1:length(varname)
    ii = strmatch(varname{a},txt{1},'exact');
    if ~isempty(ii)
        names{a} = char(txt{2}(ii));
    end
    if length(ii) > 1
        %DEPTH and DEPTH_INFERRED MATCH
%         if s.meas_inf == 1
            names{a} = 'depth';
%         else
%             names{a} = 'depth_inferred';
%         end 
    end
    %now add the qc fields:
    try
        ii = strfind(varname{a},'quality_control');
        if ~isempty(ii)
            vn = varname{a}(1:ii-2);
            ii = strmatch(vn,txt{1},'exact');
            if length(ii) > 1
                %DEPTH and DEPTH_INFERRED MATCH
%                 if s.meas_inf == 1
                    names{a} = 'depth_qc';
%                 else
%                     names{a} = 'depth_inferred_qc';
%                 end
            elseif ~isempty(ii)
                names{a} = [char(txt{2}(ii)) '_qc'];
            end
        end
    catch
    end
end
%clean it up:
b=1;idel=[];
for a = 1:length(varname)
    if isempty(names{a})
        idel(b) = a;
        b=b+1;
    end
end
varname(idel) = [];
names(idel) = [];
end
