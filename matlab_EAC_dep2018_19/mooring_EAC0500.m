% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
function mooring_EAC0500
dist =  600;
moorn = 'EAC0500';
dirn = moorn;
homedir='/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/';
inputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/data_processing/matdata_qcd_toolbox/'];
inputdir2=['/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/data_processing/matdata_qcd_toolbox/'];
doutputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/stacked/'];
poutputdir = '/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/data_processing/plots/';
%% Go through each instrument, starting with the ADCP to get the time base
cd([homedir 'data_processing'])
%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

% % Use 300khz RDI. This data has been bin averaged and has time stamp in
% middle of bin.
load([inputdir '14254'])
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
