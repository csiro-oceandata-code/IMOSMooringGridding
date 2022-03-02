% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
clear
% moorn = 'SEQ400';
moorn = 'EAC0500';
dist =  60;
dirn = moorn;
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/othermooring/SEQ/';
inputdir=[homedir 'matdata_qcd/'];
inputdir2=[homedir 'matdata_qcd/'];
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'plots/'];
%% Go through each instrument, starting with the ADCP to get the time base

%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

%set up some variables:
t_lag = [];t_cor = [];t=[];dept = [];namet={};named={};
names = {}; nameu = {};sal=[];deps=[];u=[];depu=[];
pdept=[];pdepu=[];pdeps=[];
% Use 17055
load([inputdir '17055'])
s = clean_data(s);
tbase = s.time;
temp = s.temperature;
dep = s.depth;

%nominate instruments for inferring pressure:
%No partial pressures required:
parprs = {'','',''};

%No AQDs with pressure offsets
offprs = {''};
    
%format is instrument with no press, pressure above, pressure below
% or instrument with no press, pressure closest, pressure second closest
prs = {'4100','9336','9337' %now starmon minis.
    '4101','9337','17055'
    '4102','9339','9340'
    };
%% now use the common stacking code:
stack_mooring
