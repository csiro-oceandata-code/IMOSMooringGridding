% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
function mooring_EAC4700
dist =  600;
moorn = 'EAC4700';
dirn = moorn;
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1611_1805/';
inputdir=[homedir 'data_processing/matdata_qcd_toolbox/'];
inputdir2=[homedir 'data_processing/matdata_qcd_toolbox/'];
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'data_processing/plots/'];
%% Go through each instrument, starting with the ADCP to get the time base
cd([homedir 'data_processing'])
%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

%load([inputdir '9175'])
load([inputdir '16374'])
s = clean_data(s);
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