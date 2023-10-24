% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
function mooring_EAC4200
dist =  600;
moorn = 'EAC4200';
dirn = moorn;
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1505_1611/';
inputdir=['/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1505_1611/data_processing/matdata_qcd_toolbox/'];
inputdir2=['/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1505_1611/data_processing/matdata_qcd_toolbox/'];
doutputdir=['/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1505_1611/stacked/'];
poutputdir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1505_1611/data_processing/plots/';
%% Go through each instrument, starting with the ADCP to get the time base
cd([homedir 'data_processing'])
%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

% Use 16429
load([inputdir '17054'])
%keep the temperature data for this exercise:
s = clean_start_end(s);
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
