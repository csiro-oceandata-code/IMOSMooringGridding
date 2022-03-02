%vert_extrap_mooring_EAC.m
% Script initially developed by Susan Wijffels and edited many times by Bec
% Cowley.
% This code will interpolate, filter and plot the raw mooring data for the
% EAC products.
%
% Comment in/out parts of the first cell as required to set up for
% different moorings
% Bec Cowley, February, 2022

clear all

  % ALL EAC deployments. SEQ400m mooring is now included in EAC0500 mooring
 isnsi = 0;isseq = 0;
depn = {'SEQ','EAC1204_1308','EAC1505_1611','EAC1611_1805','EAC1805_1909','EAC1909_2105'};
ndep = length(depn);

inputdir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/';
outputdir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/';
outputdirplots='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/plots/';

%Note that the SEQ400m mooring is now included in the EAC dataset and
%is labelled as EAC0500.
mor = {'EAC1520','EAC0500','EAC2000','EAC3200','EAC4200','EAC4700','EAC4800'};

%add an offset for the interpolation step
depoff = [20,20,20,20,20,20,20];

%FOR THE NSI mooring
% isnsi = 1;isseq = 0;
% depn = {'201204' '201209' '201302' '201311' '201405' '201410'	'201503' '201509'	'201602'	'201606' '201610' 	'201702'	'201706'...
%     '201710' '201803' '201806'	'201810'   '201902'	'201907', '201912','202010','202103','202106'};
% ndep = length(depn); 
% inputdir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/othermooring/NSI/';
% outputdir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/';
% outputdirplots='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/plots/';
% mor = {'NRSNSI'};
% 
% %add an offset for the interpolation step
% depoff = repmat(5,ndep);

% %FOR THE SEQ mooring, just the 200m one. 400m mooring is now part of the
% %EAC500 (see above).
% %deal with EAC deployments since May 2015 - Sept 2019 (3 deployments)
% isnsi = 0;isseq = 1;
% depn = {'SEQ'};
% ndep = length(depn); 
% inputdir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/';
% outputdir='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/';
% outputdirplots='/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/plots/';
% mor = {'SEQ200'};
% 
% %add an offset for the interpolation step
% depoff = [5];


%% Get the time stamps, an hourly one and a daily one
%builds matrices of information that are mooring (rows) * deployment
%(columns)
[xfinal,yfinal,botd,mintime,maxtime,maxd] = deal(NaN*ones(length(mor),max(ndep)));
for b = 1:length(mor)
    for a = 1:ndep
        if ~isnsi
            try
                load([inputdir '/' depn{a} '/stacked/' mor{b} 'mooring.mat' ]);
            catch
            disp(['No file: ' inputdir '/' depn{a} '/stacked/' mor{b} 'mooring.mat' ]);
            continue
            end
        else
            load([inputdir '/stacked/' depn{a} 'mooring.mat' ]);
%             disp(['No file: ' inputdir '/' depn{a} '/stacked/' mor{b} 'mooring.mat' ]);
%             continue
        end
            xfinal(b,a) = xmoor;
            yfinal(b,a) = ymoor;
            botd(b,a) = botdepth;
            mintime(b,a) = min(tbase); %first in
            maxtime(b,a) = max(tbase); %last out
            maxd(b,a) = max([max(dept),max(depu)]);
    end
end

dt = 1/24;
time_hrly = floor(min(min(mintime))):dt:max(max(maxtime));
time_hrly = time_hrly(:); %hourly grid
time_daily = floor(time_hrly(1))+0.5:1:ceil(time_hrly(end))+0.5;%daily grid
maxd = max(maxd,[],2,'omitnan');
%% Now start the interpolations
warning('on','all')
orig_state=warning; %turn off warnings

% cycle through each mooring
 for im=1:length(mor)
    moorn = mor{im};
    disp(moorn)

    % now put all deployments onto the common depth grid:
    if exist('maxd','var')
        if ~isnsi & ~isseq
            di = [0:10:400,420:20:maxd(im)];
        else
            di = [0:10:maxd(im)]; %for NSI and seq
        end
        ext = '_maxD';
    else
        if ~isnsi & ~isseq
            di = [0:10:400,420:20:max(botd(im,:))];
        else 
            di = [0:10:max(botd(im,:))];%FOR NSI
        end
        ext = '_bottom';
    end
   
    for idep =1:ndep
        disp(['Deployment ' depn{idep}])
        if ~isnsi & ~isseq
            if ~exist([inputdir '/' depn{idep} '/stacked/' moorn 'mooring.mat' ])
                continue
            else
                load([inputdir '/' depn{idep} '/stacked/' moorn 'mooring.mat' ]);
            end
        elseif isseq
            load([inputdir '/' depn{idep} '/stacked/' moorn 'mooring.mat' ]);
        else
            load([inputdir '/stacked/' depn{idep} 'mooring.mat' ]);
        end
       
        %set up empty matrices
        ui = NaN*ones(length(tbase),length(di))*[1+i];
        [ti,si ] = deal(NaN*ones(length(tbase),length(di)));
        
        
        [nt,nu]=size(u);
        [nt,ntp]=size(t);
                
        
        % set up for masking missing values. First, need to figure out the
        % depths that are missing, and re-fill them
        %find missing values
        filleddepu = fillmissing(depu,'linear',2,'endvalues','nearest');
        filleddept = fillmissing(dept,'linear',2,'endvalues','nearest');
        filleddeps = fillmissing(deps,'linear',2,'endvalues','nearest');
        
       for j=1:length(tbase)
            
            if any(filleddepu(j,:) > 0)
                jdep = filleddepu(j,:);
                ju = u(j,:);
                %let's straight away just use unique values as filleddeps might have duplication at the edges
                [jdep,iunique] = unique(jdep);
                ju = ju(iunique);
                
                %trim it to data that has depth > 0
                ig = jdep < 0;
                jdep(ig) = [];
                ju(ig) = [];
                %should be no nans in jdep (as it was filled), so just need
                %to identify data with values in u and > depth 0
                ig = find(~isnan(real(ju)) & jdep > 0 );
                
                if length(ig) > 2
                    %                     ui(j,:) = akima(depu(j,ig)',u(j,ig)',di');
                    % put flow = 0 at the bottom:
                    if botdepth < 600
                        idd = find(di <= botdepth);
                    else
                        idd = find(di < [botdepth - 60]);
                    end
                    % this is an alternative interpolation that was done in
                    % early iterations and was built by Susan Wijffels.
                    % ui(j,idd) = akima([jdep(ig)';jdep(ig(end))'+depoff(im);botd(im,idep)],[ju(ig),ju(ig(end)),[0+i*0]].',di(idd)');
                    dd = [jdep(ig)';jdep(ig(end))'+depoff(im)];
                    jju = [ju(ig),ju(ig(end))].';
                    % simple linear interpolation
                    ui(j,idd) = interp1(dd,jju,di(idd)');
                    
                    %mask out missing data chunks 
                    igz = ~isnan(ju);
                    %then we can do the mask using our filled depth values:
                    mask = interp1(jdep,double(igz),di);
                    ib = mask > 0.2;
                    ui(j,~ib) = (1+i)*NaN;
                                        
                    [msg,warnID]=lastwarn;
                    warnStruct = warning('off',warnID);
                end
            end
            
            if any(dept(j,:) > 5)
                jdep = filleddept(j,:);
                jt = t(j,:);
                %let's straight away just use unique values as filleddeps might have duplication at the edges
                [jdep,iunique] = unique(jdep);
                jt = jt(iunique);
                
                %trim it to data that has depth > 0
                ig = jdep < 0;
                jdep(ig) = [];
                jt(ig) = [];
                
                ig = find(~isnan(jt));
                if length(ig) >= 2
                    ti(j,:) = interp1q([jdep(ig)';jdep(ig(end))'+depoff(im)],[jt(ig)';jt(ig(end))],di')';
                    %then we can do the mask using our filled depth values:
                    igz = ~isnan(jt);
                    mask = interp1(jdep,double(igz),di);
                    ib = mask > 0.5;
                    ti(j,~ib) = NaN;
                end
            end
            
            %salinity
            if ~isempty(deps)
            if any(deps(j,:) > 5) 
                jdep = filleddeps(j,:);
                js = sal(j,:);
                %let's straight away just use unique values as filleddeps might have duplication at the edges
                [jdep,iunique] = unique(jdep);
                js = js(iunique);
                
                %trim it to data that has depth > 0
                ig = jdep < 0;
                jdep(ig) = [];
                js(ig) = [];
                
                ig = find(~isnan(js));
                
                if length(ig) > 2
                    si(j,:) = interp1q([jdep(ig)';jdep(ig(end))'+depoff(im)],[js(ig)';js(ig(end))],di')';
                    %then we can do the mask using our filled depth values:
                    igz = ~isnan(js);
                    mask = interp1(jdep,double(igz),di);
                    ib = mask > 0.99;
                    si(j,~ib) = NaN;
                end
            end
            end
        end  % j loop
  
        warning(orig_state);

        
        %
        
        eval(['t' num2str(idep) ' = ti;'])
        eval(['u' num2str(idep) ' = ui;'])
        eval(['s' num2str(idep) ' = si;'])
        eval(['ndate' num2str(idep) ' = tbase;'])
    end % end deployment loop
     
    % merge deployments on hourly time grid:
    
    ui = NaN*ones(length(time_hrly),length(di))*[1+i];
    [ti,si ] = deal(NaN*ones(length(time_hrly),length(di)));
    ttim = [];uu = []; tt = [];ss = [];utim = [];stim = [];
    for idep = 1:ndep
        if isnan(xfinal(im,idep))
            continue
        end
        eval(['ttim = [ttim; ndate' num2str(idep) '];'])
        eval(['utim = [utim; ndate' num2str(idep) '];'])
        eval(['stim = [stim; ndate' num2str(idep) '];'])
        eval(['uu = [uu; u' num2str(idep) '];'])
        eval(['tt = [tt; t' num2str(idep) '];'])
        eval(['ss = [ss; s' num2str(idep) '];'])
    end
    
    %interpolate across time, at each depth
    % get rid of no data times (mask)
    
    for id =1:length(di)
        
        z = [uu(:,id)];
        ig = ~isnan(z);
         if sum(ig) > 1*30*24 % need 1 months of data
            ui(:,id)=interp1(utim(ig),z(ig),time_hrly(:));
            mask = interp1(utim ,double(ig),time_hrly(:));
            ib = find(mask < 0.2);
            ui(ib,id) = (1+i)*NaN;
            
         end
        z = [tt(:,id)];
        ig = ~isnan(z);
         if sum(ig) > 1*30*24 % need 1 months of data
            ti(:,id)=interp1(ttim(ig),z(ig),time_hrly(:));
            mask = interp1(ttim,double(ig),time_hrly(:));
            
            ib = find(mask < 0.2);
            ti(ib,id) = NaN;
         end
        
        z = [ss(:,id)];
        ig = ~isnan(z);
         if sum(ig) > 1*30*24 % need 1 months of data
            si(:,id)=interp1(stim(ig),z(ig),time_hrly(:));
            mask = interp1(stim,double(ig),time_hrly(:));
            
            ib = find(mask < 0.2);
            si(ib,id) = NaN;
         end
    end
    
    % Set Nan for any time periods outside the original deployment times (between deployments)
    if ndep > 1;
        for indep=1:ndep-1
            outofrange = time_hrly > maxtime(im,indep) & time_hrly < mintime(im,indep+1);
            ti(outofrange,:) = NaN;
            si(outofrange,:) = NaN;
            ui(outofrange,:) = NaN+NaN*i;
        end
    end
   
    %record mooring lat/long information
    xmoor = xfinal(im,1:ndep);
    ymoor = yfinal(im,1:ndep);

    
    save([outputdir mor{im} '_vert_extrap' ext '.mat'],'di','time_hrly','ui','ti','si','moorn','xmoor','ymoor')
    
 
 

    % do a spectra as function of depth:
    %optional plotting, set 'dospec' to 1 if you want it to run
    dospec = 0;
    if dospec
        nsampfreq = 1./[diff(time_hrly(1:2))*24];
        % U,V and T
        
        time_hrly = time_hrly(:);
        for iv = 1:3
            switch iv
                case 1
                    zz = real(ui);
                case 2
                    zz = imag(ui);
                case 3
                    zz =ti;
            end
            
            % fill gaps and detrend
            for j=1:length(di);
                ig = (~isnan(zz(:,j))');
                if sum(ig) > length(time_hrly)*0.7
                    zz(ig,j) = detrend(zz(ig,j));
                    zz(~ig,j) = interp1q(time_hrly(ig), zz(ig,j),time_hrly(~ig));
                else
                    zz(:,j) = NaN*zz(:,j);
                end
            end
            
            
            
            idg = find(~all(isnan(zz)));
            ig = (all(~isnan(zz(:,idg))'));
            if sum(ig) > length(time_hrly)*0.7
                
                %[P,freq] = specvp_mt(detrend(zz(ig,idg(1))'),nsampfreq,ymoor);
                nsamp = floor(sum(ig)/2);
                [P,freq] = specvp(detrend(zz(ig,idg(1))'),nsampfreq,nsamp,floor(nsamp/2),ymoor);
                junk = NaN*ones(length(freq),length(di));
                
                figure(1)
                % clf
                for j=idg;
                    clf
                    % [P,freq] = specvp_mt(detrend(zz(ig,j)),nsampfreq,ymoor);
                    [P,freq] = specvp(detrend(zz(ig,j)),nsampfreq,nsamp,floor(nsamp/2) ,ymoor);
                    junk(:,j) = P(:,1);
                end
                
                
            else
                P = NaN*ones(length(time_hrly/2),length(di));
            end   % ig test
            
            switch  iv
                case 1;
                    PU = junk; freqU = freq;
                case 2
                    PV = junk;
                case 3
                    PT = junk; freqT= freq;
            end
            
            
        end % iv loop
        
%         eval(['save   ',moorn,'_interp_spectra_all2.mat mtime freqU PU PV PT moorn xmoor ymoor freqT'])
        
%%        % plot spectra
        
        idg = all( ~isnan(PU)  );
        ipf = 2:length(freqU);
        
        figure(2),clf
        clmap(7)
        subplot(311),
        pcolor(log10(freqU(ipf)),di(idg),PU(ipf,idg)')
        caxis([0.,0.003])
        xlab = {'120';'60';'30';'10';'5';'2.5';'1';'12hr';'6hr';'3hr'};
        xtck = log10(1./[[120,60,30,10,5,2.5,1]*24,12,6,3]);
        set(gca,'Xtick',xtck)
        set(gca,'XtickLabel',xlab)
        axis ij, shading flat
        title(['UP Spectra for U at ',moorn])
        
        subplot(312),
        pcolor(log10(freqU(ipf)),di(idg),PV(ipf,idg)')
        caxis([0.,0.003])
        xlab = {'120';'60';'30';'10';'5';'2.5';'1';'12hr';'6hr';'3hr'};
        xtck = log10(1./[[120,60,30,10,5,2.5,1]*24,12,6,3]);
        set(gca,'Xtick',xtck)
        set(gca,'XtickLabel',xlab)
        axis ij, shading flat
        title(['VP Spectra for V at ',moorn])
        
        idg = all( ~isnan(PT)  );
        ipf = 2:length(freqT);
        subplot(313),
        pcolor(log10(freqT(ipf)),di(idg),PT(ipf,idg)')
        caxis([0.,.1])
        xlab = {'120';'60';'30';'10';'5';'2.5';'1';'12hr';'6hr';'3hr'};
        xtck = log10(1./[[120,60,30,10,5,2.5,1]*24,12,6,3]);
        set(gca,'Xtick',xtck)
        set(gca,'XtickLabel',xlab)
        axis ij, shading flat
        title(['VP Spectra for T at ',moorn])
        
        orient tall
        
%         print('-dpng',[outputdir,moorn,'_interp_spectra.png'])
        %
    end % dospec
    
    % now low pass
    clear range
    
    % now filter the data at nf days: choose 5 days
    uiff = [NaN+ i*NaN]*ones(length(time_hrly),length(di));
    tiff = NaN*ones(length(time_hrly),length(di));
    siff = NaN*ones(length(time_hrly),length(di));
    nf = ceil(5/diff(time_hrly(1:2)));
    for j=1:length(di);
        ig = ~isnan(ui(:,j));
        if sum(ig) > 2*nf;
            uiff(ig,j) = filt_ends(hamming(nf),ui(ig,j),1);
        end
        ig = ~isnan(ti(:,j));
        if sum(ig) > 2*nf;
            tiff(ig,j) = filt_ends(hamming(nf),ti(ig,j),1);
        end
        ig = ~isnan(si(:,j));
        if sum(ig) > 2*nf;
            siff(ig,j) = filt_ends(hamming(nf),si(ig,j),1);
        end
    end
    
    
    % interpolate onto a daily timebase
    time_daily = time_daily(:);
        
    uif = interp1q(time_hrly,uiff,time_daily);
    
    tif = interp1q(time_hrly,tiff,time_daily);
    
    sif = interp1q(time_hrly,siff,time_daily);

    %fill the NaN+0i in uif with NaN+NaNi
    inan = isnan(uif);
    uif(inan) = NaN+NaN*i;
    
    %check here for temperature inversions and remove them
    dtif = diff(tif,1,2);
    iinvers = dtif>0;
    iinvers = [false(size(tif,1),1), iinvers];
    tif(iinvers) = NaN;

    save([outputdir,mor{im},'_daily_interp' ext '.mat'],'di','moorn','time_daily','uif','tif','sif','xmoor','ymoor')
    
    % make a daily average version:
    uif = NaN*ones(length(time_daily),length(di))*[1+i];tif = uif;sif = uif;
    for iday = 1:length(time_daily)
        %+/-12 hours each side of the timestamp
        iav = find(time_hrly > time_daily(iday) - 0.5 & time_hrly < time_daily(iday) + 0.5);
        uif(iday,:) = nanmean(uiff(iav,:));
        tif(iday,:) = nanmean(tiff(iav,:));
        sif(iday,:) = nanmean(siff(iav,:));
    end
    %fill the NaN+0i in uif with NaN+NaNi
    inan = isnan(uif);
    uif(inan) = NaN+NaN*i;
    save([outputdir,mor{im},'_daily_mean' ext '.mat'],'di','moorn','time_daily','uif','tif','sif','xmoor','ymoor')
    
    
    mtimef = time_daily;
    
    %idec = 1:length(mtimef);
    idec = ~all(isnan(real(uif')));
    % plot interpolated data
    
    figure(3)
%     clmap(23)
    clf
    subplot(311)
    contourf(mtimef,di,real(uif)'*100,[-100:10:100]),caxis([-50,50]),shading flat,
    
    hold on
%     contour(mtimef(idec),di,real(uif(idec,:))'*100,[-50:10:50],'w'),
    [c,h]=contour(mtimef,di,real(uif)'*100,[0:0.01:0.01],'w');
    axis ij,axis( [range(mtimef),0.,max(di)+20]),datetick('x','m','keeplimits')
    set(h,'linewidth',1,'linestyle','-'),
    title(['Low pass (5 days) U at ',mor{im}])
%     colorbar
    
    subplot(312)
    contourf(mtimef,di,imag(uif)'*100,[-100:10:100]),caxis([-50,50]),shading flat,
    hold on
%     contour(mtimef(idec),di,imag(uif(idec,:))'*100,[-50:10:50],'w'),
    [c,h]=contour(mtimef,di,imag(uif)'*100,[0:0.01:0.01],'w');
    axis ij,axis( [range(mtimef),0.,max(di)+20]),datetick('x','m','keeplimits')
    set(h,'linewidth',1,'linestyle','-'),
    title(['Low pass (5 days) V at ',mor{im}])
%     colorbar
    
    idec = ~all(isnan(tif'));
    subplot(313)
    contourf(mtimef,di,tif',[0:1:30]),caxis([3,30]),shading flat,
    hold on
    contour(mtimef,di,tif',[0:5:30],'w'),
    [c,h]=contour(mtimef,di,tif',[0:10:30],'w');
    axis ij,axis( [range(mtimef(idec)),0.,max(di)+20]),datetick('x','m','keeplimits')
    set(h,'linewidth',2),
    title(['Low pass (5 days) T at ',mor{im}])
%     colorbar
    
    orient landscape
    print('-djpeg',[outputdirplots '/' mor{im} ,'_interp_lowpass' ext '.jpg'])
    
    
    % .........................................................................
   ................
    
end  % mooring loop

