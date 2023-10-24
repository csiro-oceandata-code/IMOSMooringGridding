% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
function mooring_EAC4200
dist =  600;
moorn = 'EAC4200';
dirn = moorn;
homedir='/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/';
inputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/data_processing/matdata_qcd_toolbox/'];
inputdir2=['/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/data_processing/matdata_qcd_toolbox/'];
doutputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/stacked/'];
poutputdir = '/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/data_processing/plots/';
%% Go through each instrument, starting with the ADCP to get the time base
cd([homedir 'data_processing'])
%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

% Use RDI at 120m 
load([inputdir '16374'])
s = clean_data(s);
% Create hourly time base with the correct start and end times
%         tstep = 1.0/24.0;
%         tbase = s.time(1):tstep:s.time(end);
%         temp = interp1(s.time,s.temperature,tbase);
%         dep = interp1(s.time,s.depth,tbase);
%         tbase = double(tbase(:)); dep = double(dep(:));temp = double(temp(:));
tbase = s.time;
temp = s.temperature;
dep = s.depth;
%trim time ends:
inwater = tbase>s.starttime & tbase < s.endtime;
tbase = tbase(inwater);
temp = temp(inwater);
dep = dep(inwater);

%% now use the common stacking code:
stack_mooring

