%% run all the stacking code for each deployment in the EAC 
% Run NSI independently as it needs coaching along
clear
depdir = {'EAC1204_1308','EAC1505_1611','EAC1611_1805',...
    'EAC1805_1909','EAC1909_2105','EAC2105_2207'};% 'othermooring/NSI',
cdir = {'EAC_dep2012_13','EAC_dep2015_16','EAC_dep2016_18',...
    'EAC_dep2018_19','EAC_dep2019_21','EAC_dep2021_22'};%'NSI',
indir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/';
codedir = '/tube1/cow074/Documents/moorings/IMOSMooringGridding/matlab_';


mdir = {'SEQ400','EAC1520','EAC0500','EAC2000','EAC3200','EAC4200','EAC4700','EAC4800'};%'NSI',

for istack = 1:length(depdir)

    for im = 5:length(mdir)
        disp([cdir{istack} '/mooring_' mdir{im}])
        stackc=[codedir cdir{istack} '/mooring_' mdir{im} '.m'];
        if ~isfile(stackc)
            continue
        end
        stackc = [codedir cdir{istack} '/mooring_' mdir{im}];
        run(stackc)
    end
end