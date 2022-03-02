%% this script is the main stacking code that is common to all moorings and deployments.
% It is required that the setup script for the individual mooring is run
% first, then this script gets called. Eg: mooring_EAC4200.m
% Updated February, 2021 to remove redundant code after QC processes were
% transferred to the IMOS toolbox in 2020.
% Bec Cowley, CSIRO Oceans and Atmosphere

%% go through each instrument on the mooring:
%set up some variables
tdepth_m = [];bindex = [];isup=[];

% For each instrument on the mooring (as read in from the csv file of
% metadata)
for a = 1:size(ins.serial,1)
    %clear time lag variables in case re-running.
    tlag = NaN;maxcor = NaN; %time lag = tlag, maximum correlation = maxcor
    bdepth = []; % bin depth variable for ADCPs
    try
        %load the mat file from the conversion of the FV01 IMOS files
        load([inputdir strtrim(ins.serial(a,:))])
    catch
        disp(['No file for ' ins.serial(a,:)])
        continue
    end
    
    disp([s.serial ', ' s.name])
        
    %fix single if it's a problem:
    try
        s.temperature = double(s.temperature);
    catch
    end
            
    %make the name label
    nam = {[strtrim(s.name) ' ' s.serial]};
    
    % get the planned depth info for later sorting:
    pdep = s.planned_depth;
    
    %     clean up the data - qc has already been done. Apply the flags.
    s = clean_data(s);
    
    %skip this instrument if no valid data in any field
    if isfield(s,'temperature') %all instruments have temperature
        % check for bad depth, all have depth either calculated or inferred
        if isfield(s,'depth') 
            if sum(~isnan(s.depth)) == 0 & sum(~isnan(s.temperature)) == 0
                continue
            end
        end
    end
    
    % check for the presence of each field:
    %this also gets the time lag information
    if isfield(s,'temperature') %all instruments have temperature
        % look for time lag using pressure as this is the most accurate
        % way of getting a time lag
        if isfield(s,'pressure') 
            if ~isnan(nansum(s.pressure))
                newp =match_timebase(tbase,s.time,s.pressure);
                ins.meas_inf(a) = 1; %record that the instrument measured pressure
            else
                ins.meas_inf(a) = 0; %record that this is an inferred depth
                %the pressures have been flagged as bad, rescue them for
                %time lag calculation only!
                % won't work for every one as some are really crap
                sp = s;
                load([inputdir strtrim(ins.serial(a,:))])
                s = clean_start_end(s);
                newp =match_timebase(tbase,s.time,s.pressure);
                s=sp;
            end
            %might crash here
            try
                [tlag,maxcor]=find_timeoffset(tbase,dep,newp);
                clear newp
                if abs(tlag) > 0.2 %if there is still bad pressure data causing problems, don't use it
                    tlag = NaN;
                end
            catch
                disp(['bad pressures for ' s.serial ', ' s.name])
            end
        else
            ins.meas_inf(a) = 0; % no pressure, inferred
        end
        %now try using temperature to calculate tlag if we haven't done it
        %above
        %Also, interpolate the temperature onto the time base (tbase)
        if ~isnan(nansum(s.temperature))
            newt =match_timebase(tbase,s.time,s.temperature);
            if isnan(tlag) %if we haven't already calculated the tlag with pressure in previous loop
                ins.meas_inf(a) = 0;
                if abs(s.planned_depth - nanmedian(dep)) > dist
                    %try using temperature from the instrument above (already
                    %on tbase). gives a better match of temperatures.
                    try
                        ig = find(~isnan(nanmean(t)));
                        [tlag,maxcor]=find_timeoffset(tbase,t(:,ig(end)),newt);
                    catch
                        [tlag,maxcor] = deal(NaN);
                    end
                else
                    [tlag,maxcor]=find_timeoffset(tbase,temp,newt);
                end
            end
            t= [t,newt]; %concatenate the interpolated temperature to the temp matrix
            clear newt
        else %no valid temperature data, fill the space with NaNs
            %Might not get here if we skipped due to no valid data in
            %previous check
           t = [t,NaN*t(:,end)];
           tlag = NaN;maxcor=NaN;
        end
        
        if isnan(tlag)
            %didn't do any timebase matching
            %this happens sometimes and that's Ok, just means we don't have
            %the information about time lags. Data is still interpolated
            %and stacked.
            disp([strtrim(s.name) ', ' strtrim(s.serial) ': Did not match timebase with either pressure or temperature.'])
        end
        
        %concatenate the time lag, correlation, instrument name and planned
        %depth to previous instrument info
        t_lag = [t_lag,tlag];
        t_cor = [t_cor,maxcor];
        namet(end+1) = nam; %names of all the temperature instruments
        pdept = [pdept,pdep];
    else
        %oops, no temperature, shouldn't get here at all, but a backup
        disp('Error, no temperature field!')
        break
    end
    
    %Now the depth data to be interpolated
    if isfield(s,'depth')
        if ~isnan(nansum(s.depth))
            newp =match_timebase(tbase,s.time,s.depth);
            dept = [dept,newp(:)]; %this is the depth for every temperature record
            named(end+1) = nam; %names of all the instruments with depth
        end
    else
        %oops, no depth, shouldn't get here at all, but a backup
        disp('Error, no depth field!')
        break
    end
    
    if isfield(s,'salinity') %no need for an else in this one, not every instrument has salinity
        if ~isnan(nansum(s.salinity))
            news =match_timebase(tbase,s.time,s.salinity);
            sal= [sal,news];
            deps = [deps,newp(:)];
            names(end+1) = nam;
            pdeps = [pdeps,pdep];
        end
    end
    
    %now velocities. u is complex (so also contains v)
    if isfield(s,'u')
        %RDI's here as they are profiling and have multiple bins of data
        if ~isempty(strmatch('RDI',s.type)) 
            bsize = size(s.bdepth,2);
            %RDI, has binned data, have to loop through each bin:
            sn.bdepth = NaN*ones(length(tbase),bsize);
            sn.u = sn.bdepth*(1+NaN*i);sn.w = sn.bdepth;
            for xx = 1:bsize
                sn.bdepth(:,xx) =match_timebase(tbase,s.time,s.bdepth(:,xx));
                sn.u(:,xx) =match_timebase(tbase,s.time,s.u(:,xx));
                sn.w(:,xx) =match_timebase(tbase,s.time,s.w(:,xx));
            end
            
            %just record bins with some good data:
            idat = ~isnan(nansum(sn.u,1));
            u = [u,sn.u(:,idat)];
            depu = [depu,sn.bdepth(:,idat)];
            nameu(end+1:end+sum(idat)) = repmat(nam,sum(idat),1);
            %up or down looking?
            pdepu = [pdepu,pdep-s.brange(idat)'];
            if s.brange(2) < 0
                isup = [isup,repmat(0,1,length(s.brange(idat)))];
            else
                isup = [isup,repmat(1,1,length(s.brange(idat)))];
            end
            %record which bins are closest to the instrument
            bindex = [bindex,1:length(s.brange(idat))];
        else
            %Nortek instruments are not profiling for the EAC moorings
            newu =match_timebase(tbase,s.time,s.u);
            u = [u,newu];
            depu = [depu,newp(:)];
            nameu(end+1) = nam;
            pdepu = [pdepu,pdep];
            bindex = [bindex,1];
            isup = [isup,1];
        end
        
    end
    %keep track of measured or inferred depth information to help with
    %diagnosis if the depths look odd
    tdepth_m = [tdepth_m,ins.meas_inf(a)];
    
end

%let's exclude the aquadopp and RDI temperatures from the temp matrix
%these are low precision and can cause temperature inversions in the
%stacked data where seabird instruments are immediately adjacent on the
%mooring.
irem = find(cellfun(@isempty,strfind(upper(namet),'AQUA'))==0);
t(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'RDI'))==0);
t(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'LR'))==0);
t(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'WH'))==0);
t(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'BB'))==0);
t(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];


% reorder all the matrices, using planned depth
%temperature matrices
[val,is]=sort(pdept);

dept = dept(:,is);
dept = double(dept);
t = t(:,is);
namet = namet(is);
named = named(is);
t_cor = t_cor(is);
t_lag = t_lag(is);
tdepth_m = tdepth_m(is);

%salinity matrices
[val,is]=sort(pdeps);
deps = deps(:,is);
s = s(:,is);
names = names(is);

% and for u matrices
[val,is]=sort(pdepu);
depu = depu(:,is);
u = u(:,is);
nameu = nameu(is);
pdepu = pdepu(is);
isup = isup(is);
bindex = bindex(is);

xmoor = s.longitude;
ymoor = s.latitude;
botdepth = s.bot_depth;
%save the full stacked data, before bin removal
save([doutputdir dirn 'mooring_allbins.mat'],'tbase','u', 'name*','moorn',...
    'dep*','t','sal','tdepth_m','t_cor','t_lag','xmoor','ymoor','botdepth')

% outputdir = [poutputdir dirn '/'];
outputdir = [poutputdir '/'];

% find and remove ADCP bin that "see through" deeper instruments
% only if the data is RDI ADCP. Otherwise, remove the aquadopp if that is
% the overlap. Temporary fix till source of phase errors between
% instruments is found and hopefully fixed.

%first work out the difference between the planned depths
%then look for diffs that are less than the bin lenght of ADCP instruments.
%if the overlap is with AQD, remove the AQD record
%if the overlap is with another ADCP, keep bins that are closest to
%whichever ADCP that produced them

vel_name = unique(nameu,'stable');
irem = logical(zeros(1,size(u,2)));
if length(vel_name) > 1
    nbins = NaN*ones(1,length(vel_name));blen = nbins;
    %identify profiling instruments
    for ii = 1:length(vel_name)
        ij = strmatch(vel_name{ii},nameu);
        nbins(ii) = length(ij);
        blen(ii) = nanmean(diff(pdepu(ij))); %bin length
    end
    %find the ADCPs by number of bins
    iadcp = find(nbins > 1);
    for ii = 1:length(iadcp)
        ij = strmatch(vel_name{iadcp(ii)},nameu);
        %find where any other instrument is in between the bins:
        iother = ij(find(diff(ij)>1)+1); %other instrument index
        for b = 1:length(iother)
            if irem(iother(b)) %already flagged for removal
                continue
            end
            if ~isempty(findstr('AQUA',upper(nameu{iother(b)-1})))
                %if there is less than 75% coverage in the ADCP bin, keep
                %the aquadopp
                if sum(~isnan(real(u(:,iother(b)))))/size(u,1) > 0.75
                    %remove aquadopp, enough data in RDI
                    irem(iother(b)-1) = 1;
                else
                    if ~isup(iother(b)) %only do this for down looking
                        %remove the RDI data from this bin and others beyond it
                        irem2=strmatch(vel_name{iadcp(ii)},nameu(iother(b):end));
                        irem(irem2+iother(b)-1) = 1;
                    else %remove the aquadop that is among up-looking adcp
                        irem(iother(b)-1) = 1;
                    end
                end
            else
                %check for which RDI bins are closest to origin and keep
                if bindex(iother(b))<bindex(iother(b)-1) 
                    %remove the shallower one
                    irem(iother(b)-1) = 1;
                else
                    %remove the deeper one
                    irem(iother(b)) = 1;
                end
            end
        end
    end
    
    %perform removal of data:
    u(:,irem) = [];
    nameu(irem) = [];
    pdepu(irem) = [];
    depu(:,irem) = [];
    
end
%save version with bins removed
save([doutputdir dirn 'mooring.mat'],'tbase','u', 'name*','moorn',...
    'dep*','t','sal','tdepth_m','t_cor','t_lag','xmoor','ymoor','botdepth')
%% time correlation check - if any are low, stop and review why!!
icor = find(t_cor < .5);
if ~isempty(icor)
    disp('Time lag correlations less than 50% for instruments:')
    disp(namet(icor))
    disp('Stacking has stopped!! Possible problem with times of these instruments.')
    figure(2);clf;hold on
    nn = 1:length(t_cor);
    plot(nn,t_lag*24,'bd')
    plot(nn,t_cor,'kx-')
    try
    plot(nn(~logical(ins.meas_inf)),t_cor(~logical(ins.meas_inf)),'ro')
    catch
    end
    title('time correlations')
    legend('Time lag, hours','Correlation','Correlated with temperature')
    grid
    pause
end

%% some plots to check the data:
rawdata_plots
%pause
%% Check angles between all current meters on the moorings.
check_angles(tbase,u,depu,nameu,moorn,outputdir,dirn)

%% Check for temp drift and pressure drift between instruments
figure(2);clf
[tt,toff]= check_t_drift(tbase,t,namet,moorn);
print('-dpng',[outputdir dirn '_tdrift.png'])
figure(3);clf
[tt,poff]= check_p_drift(tbase,dept,namet,tdepth_m,moorn);
print('-dpng',[outputdir dirn '_pdrift.png'])
figure(4);clf
[tt,toff]= check_s_drift(tbase,sal,names,moorn);
return
%%
plot_u_comparison
print('-dpng',[outputdir 'ucomp.png'])
plot_v_comparison
print('-dpng',[outputdir 'vcomp.png'])

%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


