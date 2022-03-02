% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
clear
dist =  600; 
moorn = 'EAC4800';
dirn = moorn;
inputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/data_processing/matdata_qcd/'];
inputdir2=['/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/data_processing/matdata_qcd/'];
doutputdir=['/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/stacked/'];
poutputdir = '/work/observations/oceanobs_data/EACdata/mooring/EAC1805_1909/data_processing/plots/';
%% Go through each instrument, starting with the ADCP to get the time base

%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

%set up some variables:
t_lag = [];t_cor = [];t=[];dept = [];namet={};named={};
names = {}; nameu = {};sal=[];deps=[];u=[];depu=[];
pdept=[];pdepu=[];pdepd=[];pdeps=[];
load([inputdir '14890']) %
s = clean_data(s);
tbase = s.time;
temp = s.temperature;
dep = s.depth;
%trim time ends:
inwater = tbase>s.time_in & tbase < s.time_out;
tbase = tbase(inwater);
temp = temp(inwater);
dep = dep(inwater);

%% now use the common stacking code:
stack_mooring