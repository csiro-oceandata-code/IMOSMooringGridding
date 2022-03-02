% import the IMOS netcdf files and take a look:
clear

moorn = 'EAC3200';
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/data_processing/';
inputdir = ['/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/data_processing/IMOSnetcdf/' moorn '/'];
outputdir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/data_processing/matdata_qcd_toolbox/';

% homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/othermooring/NSI/data_processing/';
% m = dir(homedir);
% moorn = {};
% for a = 3:length(m)
%     if m(a).isdir
%         moorn = [moorn,m(a).name];
%     end
% end

% for fold = 3%length(moorn)
% inputdir = [homedir  moorn{fold} '/'];
% outputdir = [homedir moorn{fold} '/matdata_qcd/'];
if exist(outputdir,'dir') ~= 7
    system(['mkdir ' outputdir])
end
fn = dir([inputdir '*FV01*.nc']);
% fn = dir([inputdir '*.nc']);
% ext = {'1105','1201','1207','1212','1306','1401'};

for a = 1:length(fn)
    disp(fn(a).name)
    s = imosnc2mat([inputdir fn(a).name]);
    s.fromToolbox = 1;
%     s.bot_depth = 63;%for NSI mooring
%     if ~isempty(findstr('RDI',s.name))
%         s.type = 'RDI';
%     elseif ~isempty(findstr('WQM',s.name))
%         s.type = 'WQM';
%     else
%         s.type = 'SBE37';
%     end
    
    %BIN AVERAGE HERE for single-ping ADCP data
    if contains('RDI',s.name)
        if s.time(5) - s.time(4) < 0.01
            %ensure that the bdepth flags inherit the u_qc
            s.bdepth_qc = s.u_qc;
            s = imosRDIbinAveragingCSIRO(s);
            %remove bdepths without data in u
            inan = isnan(s.u);
            s.bdepth_qc(inan) = 4;
        end
    end
    
%         s.serial = [s.serial '_' ext{b}];
    %edit time_in and time_out fields:
    ii = findstr('T',s.time_in);
    ij = findstr('Z',s.time_in);
    if ~isempty(ii)
        s.time_in(ii) = ' ';
        s.time_in(ij) = ' ';
        ii = findstr('T',s.time_out);
        ij = findstr('Z',s.time_out);
        s.time_out(ii) = ' ';
        s.time_out(ij) = ' ';
        s.starttime = datenum(s.time_in,'yyyy-mm-dd HH:MM:SS');
        s.endtime = datenum(s.time_out,'yyyy-mm-dd HH:MM:SS');
        s.time_in = datestr(datenum(s.time_in,'yyyy-mm-dd HH:MM:SS'),'dd/mm/yyyy HH:MM:SS');
        s.time_out = datestr(datenum(s.time_out,'yyyy-mm-dd HH:MM:SS'),'dd/mm/yyyy HH:MM:SS');
    else
        s.starttime = s.time_in;
        s.endtime = s.time_out;
        s.time_in = datestr(s.time_in,'dd/mm/yyyy HH:MM:SS');
        s.time_out = datestr(s.time_out,'dd/mm/yyyy HH:MM:SS');        
    end
    %save the file:
    save([outputdir s.serial '.mat'],'s')
    
%     %extract the mooring metadata to create the instrument_info.csv file
%     fid = fopen([inputdir 'instrument_info.csv'],'a');
%     try
%     fprintf(fid,'\n%s',[s.name ',' s.name ',' s.mooring ',' s.serial ',' num2str(s.time_int) ',' ...
%         num2str(s.planned_depth) ',,' num2str(s.bot_depth) ',0,' num2str(s.latitude) ',' num2str(s.longitude) ',0,0,' ...
%         s.time_in ',' s.time_out ',01/01/2001 00:00,01/01/2001 00:00'])
%     catch
%     fprintf(fid,'\n%s',[s.name ',' s.name ',' s.mooring ',' s.serial ',' num2str(s.time_int) ',' ...
%         num2str(s.planned_depth) ',,,0,' num2str(s.latitude) ',' num2str(s.longitude) ',0,0,' ...
%         s.time_in ',' s.time_out ',01/01/2001 00:00,01/01/2001 00:00'])
%     end
%     fclose(fid)
    
end
% end
return
% %% make up a list of the mooring names:
% inp = '/home/cowley/work/moorings/anmn_timorS/matdata/';
% fn = dir([inp '*.mat']);
% clear ins
% for a = 1:length(fn)
%     load([inp fn(a).name])
%     ins.name{a} = s.name;
%     if isfield(s,'type')
%         ins.type{a} = s.type;
%     else
%         ins.type{a} = [];
%     end
%     ins.serial{a} = s.serial;
%     ins.lat(a) = s.latitude;
%     ins.lon(a) = s.longitude;
%     ins.mooring{a} = s.mooring;
%     ins.time_int(a) = s.time_int;
%     ins.planned_depth(a) = s.planned_depth;
%     ins.meas_inf(a) = s.meas_inf;
%     ins.filen{a} = fn(a).name;
% end
% % save instrument_info_TSth.mat ins
%% do some plots
fn = dir([outputdir '*.mat']);
    figure(4);clf;hold on
    coast

for a = 1:length(fn)
    load([outputdir fn(a).name])
    try
        figure(1);clf;hold on
        plot(s.u,s.bdepth,'.')
        axis ij; grid
    catch
    end
    
    figure(2);clf;hold on
    plot(s.time,s.temperature)
    grid
    
    try
        figure(3);clf;hold on
        plot(s.time,s.pitch)
        plot(s.time,s.roll,'r')
    catch
    end
    
    figure(5);clf;hold on
    try
    plot(s.time,s.depth)
    catch
        plot(s.time,s.depth_inferred,'ko')
    end
    grid
    
    s = clean_data(s);
    
    try
        figure(1)
        plot(s.u,s.bdepth,'ko')
        
        
        figure(3)
        plot(s.time,s.pitch,'ko')
        plot(s.time,s.roll,'go')
    catch
    end
    
    figure(2)
    plot(s.time,s.temperature,'ko')
    datetick('x','myyyy','keeplimits')
    title(s.serial)

    figure(4)
    plot(s.longitude,s.latitude,'x')
    title(s.serial)
    
    
    figure(5);
    try
        plot(s.time,s.depth,'ko')
    catch
        plot(s.time,s.depth_inferred,'ko')
    end
    axis ij;grid
    datetick('x','myyyy','keeplimits')
    title(s.serial)
    
    pause

end
