% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
function mooring_EAC4200
dist =  600;
moorn = 'EAC4200';
dirn = moorn;
homedir = '/work/observations/oceanobs_data/EACdata/mooring/EAC1204_1308/';
inputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1204_1308/data_processing/matdata_qcd_toolbox/'];
inputdir2=['/work/observations/oceanobs_data/EACdata/mooring/EAC1204_1308/data_processing/matdata_qcd_toolbox/'];
doutputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1204_1308/stacked/'];
poutputdir = '/work/observations/oceanobs_data/EACdata/mooring/EAC1204_1308/data_processing/plots/';
%% Go through each instrument, starting with the ADCP to get the time base
cd([homedir 'data_processing'])
%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

% Use 16432
load([inputdir '16432'])
s = clean_data(s);
tbase = s.time;
temp = s.temperature;
dep = s.depth;

%% now use the common stacking code:
stack_mooring

