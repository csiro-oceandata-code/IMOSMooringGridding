% compiles all data from single mooring onto one time base
% put u onto depth , put t onto tdepth
clear
dist =  10;
homedir = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/othermooring/NSI/';
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'plots/'];
m = dir([homedir 'data_processing/']);
moor = {};
for ifn = 3:length(m)
    if m(ifn).isdir
        moor = [moor,m(ifn).name];
    end
end
%% need to smash together WQM data for the 201204, 201209, 201302 deployments
%201204: dep 1, 044 = 60m, dep2, 044_2 = 60m; dep 1, 173 = 20m, dep2, 174 = 20m
%201209: dep1, 173 = 60m, dep 2, 044 = 60m; dep1 = 176 = 20m, dep2 176_2 = 20m
%201302: dep1 175 = 20m, dep2 176 = 20m; dep1 174 = 60m, dep2, 044 = 60m
%do this once and save the files, then move the others out of the folder:
% cd /Volumes/DWmoorings1/othermoorings/NSI/data_processing/201204/matdata_qcd
% d20 = {'173','174'};
% d60 = {'044','044_2'};
% cd /Volumes/DWmoorings1/othermoorings/NSI/data_processing/201209/matdata_qcd
% d20 = {'176','176_2'};
% d60 = {'173','044'};
% cd /Volumes/DWmoorings1/othermoorings/NSI/data_processing/201302/matdata_qcd
% d20 = {'175','176'};
% d60 = {'174','044'};
% flds ={'time','pressure','temperature','depth','salinity'};
% deps = {'d20','d60'};
% for moor = 1%:3
%     for depn = 1:2
%         for b = 1:2
%             eval(['d = ' deps{depn} ';']);
%             load(d{b})
%             if b == 1
%                 snew = s;
%             end
%             if b == 2
%                 snew.endtime = s.endtime;
%                 snew.serial = [snew.serial '_' s.serial];
%                 for c = 1:length(flds)
%                     snew.(flds{c}) = [snew.(flds{c}); s.(flds{c})];
%                 end
%                 s = snew;
%                 save([s.serial '.mat'],'s');
%             end
%         end
%     end
% end
%%
for fold =20:23%1%:length(moor)
    clear ins
    dirn = moor{fold};
    moorn = moor{fold};
    inputdir=[homedir 'data_processing/' dirn '/matdata_qcd/'];
    inputdir2=inputdir;
    
    %% Go through each instrument, starting with the ADCP to get the time base
    
    %get the serial number informations for the stacking:
    fn = dir([inputdir '*.mat']);
    for ifn = 1:length(fn)
        load([inputdir fn(ifn).name])
        ins.serial(ifn,:) = '           ';
        ins.serial(ifn,1:length(fn(ifn).name(1:end-4))) = fn(ifn).name(1:end-4);%s.serial;
%         ins.serial(ifn,1:length(s.serial)) = s.serial;            
        if isfield(s,'u')
            s = clean_data(s);
            tbase = s.time;
            temp = s.temperature;
            dep = s.depth;
            %trim time ends:
            inwater = tbase>s.starttime & tbase < s.endtime;
            tbase = tbase(inwater);
            temp = temp(inwater);
            dep = dep(inwater);
        end
    end
        if contains(dirn,'201302')
            %need to jam two temperature records from the WQM at 60m
            %together
            clear temp tbase dep
            load([inputdir '174_044.mat'])
            tbase = s.time;
            temp = s.temperature;
            dep = s.depth;
            instart = s.starttime;
            outend = s.endtime;
            inwater = tbase>instart & tbase < outend;
            tbase = tbase(inwater);
            temp = temp(inwater);
            dep = dep(inwater);
       end
    %set up some variables:
    t_lag = [];t_cor = [];t=[];dept = [];namet={};named={};
    names = {}; nameu = {};sal=[];deps=[];u=[];depu=[];
    pdept=[];pdepu=[];pdepd=[];pdeps=[];
    
    %nominate instruments for inferring pressure:
    %No partial pressures required:
    parprs = {''};
    
    %No AQDs and others with pressure offsets
    offprs = {};
    
    %format is instrument with no press, pressure above, pressure below
    % or instrument with no press, pressure closest, pressure second closest
    prs = {''
        };
    %% now use the common stacking code:
    try
        stack_mooring
        if contains(dirn,'201204')
            %nan out gap
            dept(4862:8253,1) = NaN;
            t(4862:8253,1) = NaN;
            sal(4862:8253,2) = NaN;
            deps(4862:8253,2) = NaN;
            save([doutputdir dirn 'mooring_allbins.mat'],'tbase','u', 'name*','moorn',...
    'dep*','t','sal','tdepth_m','t_cor','t_lag','xmoor','ymoor','botdepth')
            save([doutputdir dirn 'mooring.mat'],'tbase','u', 'name*','moorn',...
    'dep*','t','sal','tdepth_m','t_cor','t_lag','xmoor','ymoor','botdepth')

        end
    catch Me
        disp(Me)
        continue
    end
    
end
return
%% fix the salinity in 201209 - WQM 173 bad salinity
clear
homedir = '/Volumes/DWmoorings1/othermoorings/NSI/';
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'plots/'];
load([doutputdir '201209mooring_allbins.mat'])
sal(1:4694,1)=NaN;
save([doutputdir '201209mooring_allbins.mat'])
save([doutputdir '201209mooring.mat'])
%and 201803
load([doutputdir '201803mooring_allbins.mat'])
sal(:,2)=NaN;
save([doutputdir '201803mooring_allbins.mat'])
save([doutputdir '201803mooring.mat'])

%% check all the data stacks:
clear
homedir = '/Volumes/DWmoorings1/othermoorings/NSI/';
doutputdir=[homedir 'stacked/'];
poutputdir = [homedir 'plots/'];
labl=[];
filn = dir([doutputdir '*allbins.mat'])
figure(1);clf;hold on
figure(2);clf;hold on
figure(3);clf;hold on
figure(4);clf;hold on
for a = 1:length(filn)
    load([doutputdir filn(a).name])
    try
        figure(1)
        for b = 1:length(names)
            ii = strmatch(names{b},namet);
            plot(sal(:,b),t(:,ii),'.')
        end
        figure(2)
        plot(sal,deps,'.')
        figure(3)
        plot(tbase,sal,'.-')
        figure(4)
        plot(tbase,dept,'.-')
        labl = [labl, moorn];
    catch
    end
end
% legend(labl)
% datetick('x','mm/yy')

%% and velocity check

figure(1);clf;hold on
figure(2);clf;hold on
figure(3);clf;hold on
figure(4);clf;hold on
for a = 1:length(filn)
    load([doutputdir filn(a).name])
    figure(1)
    pcolor(tbase,depu',real(u)');shading flat;axis ij;caxis([-1,1]*0.5)
    figure(2)
    pcolor(tbase,depu',imag(u)');shading flat;axis ij;caxis([-1,1]*0.5)
    figure(3)
    plot(real(u),depu,'.'),axis ij,grid,title('U')
    figure(4)
    plot(imag(u),depu,'.'),axis ij,grid,title('V')
end
