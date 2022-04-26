% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
clear
dist =  600;
moorn = 'EAC0500';
dirn = moorn;
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1611_1805/';
inputdir=[homedir 'data_processing/matdata_qcd_toolbox/'];
inputdir2=[homedir 'data_processing/matdata_qcd_toolbox/'];
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'data_processing/plots/'];
%% Go through each instrument, starting with the ADCP to get the time base

%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

%set up some variables:
t_lag = [];t_cor = [];t=[];dept = [];namet={};named={};
names = {}; nameu = {};sal=[];deps=[];u=[];depu=[];
pdept=[];pdepu=[];pdepd=[];pdeps=[];

% % Use 9167 to get start and end times. one ADCP failed, other had short
% record
load([inputdir '9167'])
s = clean_data(s);
% Create hourly time base with the correct start and end times
tstep = 1.0/24.0;
tbase = s.time(1):tstep:s.time(end);
temp = interp1(s.time,s.temperature,tbase);
dep = interp1(s.time,s.depth,tbase);
tbase = double(tbase(:)); dep = double(dep(:));temp = double(temp(:));
%trim time ends:
inwater = tbase>s.starttime & tbase < s.endtime;
tbase = tbase(inwater);
temp = temp(inwater);
dep = dep(inwater);

%% now use the common stacking code:
stack_mooring
