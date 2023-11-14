function ins = read_ins_info(fname,moorn)
% function ins = read_ins_info(fname)
% read the instrument information file from ss2012_v05
% optional argument: moorn = mooring name to just return info from one
% mooring
% file is in csv format.
% ins = 
% 
%                 type: [84x9 char]
%                 name: [84x16 char]
%                  lat: [84x1 double]
%                  lon: [84x1 double]
%              mooring: [84x11 char]
%               serial: [84x7 char]
%             time_int: [84x1 double]
%        planned_depth: [84x1 double]
%     max_depth_rating: [84x1 double]
%            tbath_cal: [84x1 logical]
%                units: {4x2 cell}
% Bec Cowley, 29th Sept, 2012


%open the file:
[fid,msg] = fopen(fname);
if fid < 0
    disp(msg)
    disp(['fid = ' num2str(fid)])
    ins = [];
    return
end
%read it
c = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',','headerlines',1);
fclose(fid);
if nargin == 2
    ii = strmatch(moorn,char(c{3}),'exact');
else
    ii = 1:length(char(c{3}));
end
if isempty(ii)
    disp(['Mooring ' moorn ' not found. instrument information not read.'])
    ins = [];
    return
end
ins.type = char(c{1}(ii));
ins.name = char(c{2}(ii));
ins.mooring = char(c{3}(ii));
ins.lat = str2num(char(c{10}(ii)));
ins.lon = str2num(char(c{11}(ii)));
ins.serial = c{4}(ii);
ins.time_int = str2num(char(c{5}(ii)));
ins.time_offset = str2num(char(c{12}(ii)));
ins.time_offset_end = str2num(char(c{13}(ii)));
ins.planned_depth = str2num(char(c{6}(ii)));
ins.max_depth_rating = str2num(char(c{7}(ii)));
ins.bot_depth = str2num(char(c{8}(ii)));
ins.tbath_cal = logical(str2num(char(c{9}(ii))));
ins.tbath_in = datenum(char(c{16}(ii)),'dd/mm/yy HH:MM');
ins.tbath_out = datenum(char(c{17}(ii)),'dd/mm/yy HH:MM');
ins.start = datenum(char(c{14}(ii)),'dd/mm/yy HH:MM');
ins.end = datenum(char(c{15}(ii)),'dd/mm/yy HH:MM');
ins.units = {'time_int','seconds';
    'planned_depth','m';
    'max_depth_rating','m';
    'tbath_cal','1=yes, 0=no';
    'time_offset','seconds'};

%make up a array of measured/inferred depths information
% %get the number of ins on each mooring by type
% load ins_counts 
% 
% for a = 1:length(ins.lat)
%     ij = strmatch(ins.type(a,:),ins_counts.type);
%     
%     ins.meas_inf(a) = ins_counts.ins_measures(ij,3);
% end


end