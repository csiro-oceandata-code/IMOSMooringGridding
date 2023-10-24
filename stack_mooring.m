%% this script is the main stacking code that is common to all moorings and deployments.
% It is required that the setup script for the individual mooring is run
% first, then this script gets called. Eg: mooring_EAC4200.m
% Updated February, 2021 to remove redundant code after QC processes were
% transferred to the IMOS toolbox in 2020.
% Bec Cowley, CSIRO Oceans and Atmosphere
%
% Oct 2023: Include uncertainty estimates:
%   SBE37 calibration certificiates, can source this information for each
%   instrument, but for now, use these values which are from a sub-set of
%   calibration certificates:
%       PSAL: 0.0018 (converted from conductivity = 0.002 at STP using 
%               gsw_SP_from_C, use cond=42.918 - 42.916 mS/cm )
%       TEMP: 0.0015 deg C
%       PRES: 0.5 dbar + 0.01 % of reading (translate to DEPTH)
%   Star Oddis:
%       TEMP: 0.01 deg C (can be better or worse than this from calibration certs)
%       PRES: inferred, inherit SBE PRES uncertainties + 0.5 dbar for
%       measurement uncertainty (translate to DEPTH)
%   SBE39 calibration certificiates, can source this information for each
%   instrument, but for now, use these values which are from a sub-set of
%   calibration certificates:
%       TEMP: 0.0015 deg C
%       PRES: 0.5 dbar + 0.01 % of reading (translate to DEPTH)
%   Aquadopps from setup log file, instrumental uncertainty:
%       TEMP: 0.1 deg C; NOT USED IN PRODUCTS
%       PRES: 0.5 % of reading; BUT MOST HAVE LARGE BIASES AND ONLY GOOD
%           DATA USED IN PRODUCTS, INFERRED DEPTHS WILL INHERIT SBE PRES
%           UNCERTAINTIES + 0.5 dbar for measurement uncertainty (translate to DEPTH)
%       Currents: 0.009 m/s
%   RDI:
%       Currents: from standard deviation of error velocity, calculated. Could also
%           add the  instrumental uncertainty, from this formula:
%           sd(cm/s) = 1.6 x 10^7/(F * I * B * sqrt(P))
%           F = ADCP acoustic frequency in Hz, get from 'instrument' in global atts
%               153600 for a 150kHz, 307200 for a 300kHz, 76800 for a 75kHz
%           I = Transmit ppulse length in meters, get from diff(bdepth)
%           B = Beam angle coefficient (1 for 30°; 0.684 for 20°)
%           P = number of pings per ensemble
%               for 75kHz with on-board averaging and single ping averaged to hourly:
%                   3:45s intervals (225s), hourly averaging = 16 pings per ensemble
%               for 300kHz & 150kHz with on-board averaging and single ping averaged to hourly:
%                   1:12s intervals (72s), hourly averaging = 50 pings per ensemble
%       TEMP: 0.4 deg C; NOT USED IN PRODUCTS
%       PRES: 0.1 % of reading (translate to DEPTH)
%   WQM:
%       Use same values for TEMP, PSAL, PRES as SBE37
%
% Note: PRES/DEPTH uncertainties have not been checked and do not carry through
% to the data products at this stage (October, 2023)


%% go through each instrument on the mooring:
%set up some variables
tdepth_m = [];bindex = [];isup=[];
t_lag = [];t_cor = [];t=[];dept = [];namet={};named={};
names = {}; nameu = {};sal=[];deps=[];u=[];depu=[];
tunc = []; sunc= []; uunc = []; dunc = [];
pdept=[];pdepu=[];pdepd=[];pdeps=[];

% For each instrument on the mooring (as read in from the csv file of
% metadata)
for a = 1:size(ins.serial,1)
    %clear time lag variables in case re-running.
    tlag = NaN;maxcor = NaN; %time lag = tlag, maximum correlation = maxcor
    bdepth = []; % bin depth variable for ADCPs
    try
        %load the mat file from the conversion of the FV01 IMOS files
        load([inputdir strtrim(ins.serial{a})])
    catch
        disp(['No file for ' ins.serial{a}])
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
    % first use the U/V flags for ERV so we only use good data for
    % uncertainty estimates:
    if isfield(s,'u')
        s.erv_qc = s.u_qc;
    end
    s = clean_data(s);
    
    %skip this instrument if no valid data in any field
    if isfield(s,'temperature') %all instruments have temperature
        % check for bad depth, all have depth either calculated or inferred
        if isfield(s,'depth') 
            if sum(~isnan(s.depth)) == 0 & sum(~isnan(s.temperature)) == 0
                continue
            end
        end
        % create uncertainties for S, T, P variables here
        if contains(s.name, 'Star')
            s.temperature_unc = repmat(0.01, length(s.time),1);
        elseif contains(s.name,'Aqua')
            s.temperature_unc = repmat(0.1, length(s.time),1);
            s.pressure_unc = 0.5/100*s.pressure + 0.5;
        elseif contains(s.name, 'RDI')
            s.temperature_unc = repmat(0.4, length(s.time),1);
            s.pressure_unc = 0.1/100*s.pressure;   
            if contains(s.name, '150')
                F = 153600;
                P = 50;
            elseif contains(s.name, '300')
                F = 307200;
                P = 50;
            elseif contains(s.name, '75')
                F = 76800;
                P = 16;
            else
                disp('Frequency can''t be determined')
                return
            end
            I = abs(diff(s.brange(1:2)));
            B = 0.684;
            sd = 1.6 * 10^7/(F * I * B * sqrt(P));
            % convert to m/s
            sd = sd/100;
        elseif contains(s.name, 'SBE') | contains(s.name, 'WQM')
            if isfield(s,'pressure')
                s.pressure_unc = 0.01/100*s.pressure + 0.5;
            end
            s.temperature_unc = repmat(0.0015, length(s.time),1);
            if isfield(s,'salinity')
                s.salinity_unc = repmat(0.0018, length(s.time),1);
                s.salinity_unc(isnan(s.salinity)) = NaN;
            end
        end
        % clean up uncertainty arrays with qc. Pressure is based on already
        % cleaned data, salinity done above
        s.temperature_unc(isnan(s.temperature)) = NaN;
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
                load([inputdir strtrim(ins.serial{a})])
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
            newunc = match_timebase(tbase,s.time,s.temperature_unc);
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
            tunc = [tunc, newunc];
            clear newt newunc
        else %no valid temperature data, fill the space with NaNs
            %Might not get here if we skipped due to no valid data in
            %previous check
           t = [t,NaN*t(:,end)];
           tunc = [tunc, NaN*t(:,end)];
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
            if ins.meas_inf(a) == 0 %inferred depths
               dunc = [dunc, -gsw_z_from_p(0.01/100*newp + 0.5 + 0.5, s.latitude)];
            else % measured from pressure
                newpunc = match_timebase(tbase,s.time,s.pressure_unc);
                dunc = [dunc, -gsw_z_from_p(newpunc, s.latitude)];
            end
        end
    else
        %oops, no depth, shouldn't get here at all, but a backup
        disp('Error, no depth field!')
        break
    end
    
    if isfield(s,'salinity') 
        if ~isnan(nansum(s.salinity))
            news =match_timebase(tbase,s.time,s.salinity);
            newsunc = match_timebase(tbase,s.time,s.salinity_unc);
            sal= [sal,news];
            sunc = [sunc, newsunc];
            deps = [deps,newp(:)];
            names(end+1) = nam;
           pdeps = [pdeps,pdep];
        else %keep a space if the instrument has no salinity to help with masking later in gridding step
            sal = [sal,NaN*sal(:,end)];
            sunc = [sunc,NaN*sal(:,end)];
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
            sn.u = sn.bdepth*(1+NaN*i);sn.w = sn.bdepth;sn.erv = sn.bdepth;
            for xx = 1:bsize
                sn.bdepth(:,xx) =match_timebase(tbase,s.time,s.bdepth(:,xx));
                sn.u(:,xx) =match_timebase(tbase,s.time,s.u(:,xx));
                sn.w(:,xx) =match_timebase(tbase,s.time,s.w(:,xx));
                sn.erv(:,xx) = match_timebase(tbase,s.time,s.erv(:,xx));
            end
            % where there is a u value, but NaN erv, replace with zero:
            inan = isnan(sn.erv) & ~isnan(sn.u);
            sn.erv(inan) = 0;
            %just record bins with some good data:
            idat = ~isnan(nansum(sn.u,1));
            u = [u,sn.u(:,idat)];

            depu = [depu,sn.bdepth(:,idat)];
            % calculate uncertainties on a moving window of +/- 3 hours
            ut = sn.erv(:,idat);
            mstd = movstd(ut,7,0,2,"omitnan");
            % make a matrix of just instrument uncertainty
            sdm = repmat(sd,size(sn.erv(:,idat)));
            uunc = [uunc, mstd+sdm]; % uncertainty estimate
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
            uunc = [uunc, repmat(0.009,length(tbase),1)];% hard code to instrument setup log output for uncertainty (m/s)
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
% tidy uncertainties where NaNs in the data
uunc(isnan(u)) = NaN;
tuunc(isnan(t)) = NaN;
sunc(isnan(sal)) = NaN;
dunc(isnan(dept)) = NaN;

%let's exclude the aquadopp and RDI temperatures from the temp matrix
%these are low precision and can cause temperature inversions in the
%stacked data where seabird instruments are immediately adjacent on the
%mooring.
irem = find(cellfun(@isempty,strfind(upper(namet),'AQUA'))==0);
t(:,irem) = [];
tunc(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'RDI'))==0);
t(:,irem) = [];
tunc(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'LR'))==0);
t(:,irem) = [];
tunc(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'WH'))==0);
t(:,irem) = [];
tunc(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];
irem = find(cellfun(@isempty,strfind(upper(namet),'BB'))==0);
t(:,irem) = [];
tunc(:,irem) = [];
namet(irem) = [];
dept(:,irem)=[];
pdept(irem) = [];


% reorder all the matrices, using planned depth
%temperature matrices
[val,is]=sort(pdept);

dept = dept(:,is);
dept = double(dept);
dunc = dunc(:,is);
t = t(:,is);
tunc = tunc(:,is);
namet = namet(is);
named = named(is);
t_cor = t_cor(is);
t_lag = t_lag(is);
tdepth_m = tdepth_m(is);

%salinity matrices
[val,is]=sort(pdeps);
deps = deps(:,is);
sal = sal(:,is);
sunc = sunc(:,is);
names = names(is);

% and for u matrices
[val,is]=sort(pdepu);
depu = depu(:,is);
u = u(:,is);
uunc = uunc(:,is);
nameu = nameu(is);
pdepu = pdepu(is);
isup = isup(is);
bindex = bindex(is);

xmoor = s.longitude;
ymoor = s.latitude;
botdepth = s.bot_depth;
%save the full stacked data, before bin removal
save([doutputdir dirn 'mooring_allbins.mat'],'tbase','u', '*unc','name*','moorn',...
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
    uunc(:,irem) = [];
    
end
%save version with bins removed
save([doutputdir dirn 'mooring.mat'],'tbase','u', 'name*','moorn', '*unc',...
    'dep*','t','sal','tdepth_m','t_cor','t_lag','xmoor','ymoor','botdepth')
return
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
    % pause
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
print('-dpng',[outputdir dirn '_sdrift.png'])
return
%%
plot_u_comparison
print('-dpng',[outputdir 'ucomp.png'])
plot_v_comparison
print('-dpng',[outputdir 'vcomp.png'])

%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


