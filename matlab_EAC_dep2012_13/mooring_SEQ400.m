% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
function mooring_SEQ400
% moorn = 'SEQ400';
moorn = 'EAC0500';
dist =  60;
dirn = moorn;
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/SEQ/';
inputdir=[homedir 'matdata_qcd/'];
inputdir2=[homedir 'matdata_qcd/'];
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'plots/'];
%% Go through each instrument, starting with the ADCP to get the time base
cd(homedir)
%get some information about the mooring:
ins = read_ins_info('instrument_info.csv',moorn);

% Use 17055
load([inputdir '17055'])
s = clean_data(s);
tbase = s.time;
temp = s.temperature;
dep = s.depth;

%% now use the common stacking code:
stack_mooring


%% lets fix these old mat files that don't have the same format as recent
% constructions from IMOS files. I QC'd these after finding issues in the
% toolbox versions, so will stick with them, but adjust.
% one off. Bec Cowley, 20 October, 2023
% clear
% moorn = 'EAC0500';
% homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/SEQ/';
% inputdir=[homedir 'matdata_qcd/'];
% %% Go through each instrument, starting with the ADCP to get the time base
% cd(homedir)
% %get some information about the mooring:
% ins = read_ins_info('instrument_info.csv',moorn);
% for a = 1:size(ins.serial,1)
%     load([inputdir ins.serial{a}])
%     if isfield(s,'depth_inferred')
%         s.meas_inf = 0;
%         s.depth = s.depth_inferred;
%     else
%         s.meas_inf = 1;
%     end
%     save([inputdir ins.serial{a}], 's');
% end
