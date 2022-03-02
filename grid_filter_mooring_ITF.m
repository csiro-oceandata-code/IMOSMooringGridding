% extrapolate and plot the raw mooring data vertically:
clear
% Go through each instrument, starting with the ADCP to get the time base
% %EAC
% ndep = [1,1,1,1,1,1,1];
% inputdir=['/home/cowley/work/moorings/SS2013_V05/moorings/'];
% outputdir='/home/cowley/work/moorings/SS2013_V05/moorings/plots/';
% mn = {'EAC_M1'  'EAC_M2'  'EAC_M3' 'EAC_M4'  'EAC_M5'  'SEQ200'  'SEQ400'};
% mor = mn;
% omn = mn;
% depdir = {'/home/cowley/work/moorings/SS2013_V05/moorings/',...
%     '/home/cowley/work/moorings/SS2013_V05/moorings/'};

% ITF
ndep = [3,2,3,1,9,9,9,8];
depoff = [20,20,20,20,0,0,0,0];
%depth grid: set for ITF Ombai and Timor Sill to 1500m. Will fail when
%trying to run for shelf moorings. Need to code in. Remove this line if
%interpolating to the bottom.
% maxd = [1500,1200,1500,1500];
inputdir=['/Users/cow074/Documents/work_mac/moorings/ITF_stacked_files/'];
outputdir='/Users/cow074/Documents/work_mac/moorings/ITFanmn_data/plots/';
mor = {'Timor Sill','Timor North','Ombai','Timor North Slope','Timor South','JBG','MHB','FTB'};
mn = {'Timor Sill';'Timor North';'Ombai';'Timor North Slope'
    {'TimorSth_1105';'TimorSth_1201';'TimorSth_1207';'TimorSth_1212';
    'TimorSth_1306';'TimorSth_1401';'TimorSth_1409';'TimorSth_1502';'TimorSth_1508'}
    {'JBG_1105';'JBG_1201';'JBG_1207';'JBG_1212';'JBG_1306';'JBG_1401';'JBG_1409';'JBG_1502';'JBG_1508'}
    {'MHB_1105';'MHB_1201';'MHB_1207';'MHB_1212';'MHB_1306';'MHB_1401';'MHB_1409';'MHB_1502';'MHB_1508'}
    {'FTB_1105';'FTB_1207';'FTB_1212';'FTB_1306';'FTB_1401';'FTB_1409';'FTB_1502';'FTB_1508'}};
omn = {'TimorSill','TimorNorth','Ombai','TimorNorthSlope','TimorSth','JBG','MHB','FTB'};
depdir = {'/Users/cow074/Documents/work_mac/moorings/ITF_2012/'
    '/Users/cow074/Documents/work_mac/moorings/ITF_2014/'
    '/Volumes/DWmoorings1/Solander_2015/2014_recoveries/'
    '/Users/cow074/Documents/work_mac/moorings/ITFanmn_data/'};

% %EAC dep2
% ndep = [1,1,1,1,2,1,2,2,2];
% depoff = [20,20,20,20,20,20,20,20,20];
% inputdir=['U:\moorings\IN2016_V06\moorings\EAC_stacked_files\'];
% outputdir='U:\moorings\IN2016_V06\moorings\plots\';
% mor = {'SEQ200','SEQ400','EAC0500','EAC_M1','EAC2000','EAC3200','EAC4200','EAC4700','EAC4800'};
% mn = {'SEQ200';'SEQ400';'EAC0500';'EAC_M1';
%     {'EAC_M2';'EAC2000'}
%     'EAC3200';
%     {'EAC_M3';'EAC4200'}
%     {'EAC_M4';'EAC4700'}
%     {'EAC_M5';'EAC4800'}};
% omn = mor;
% depdir = {'U:\moorings\IN2016_V06\moorings\EAC_stacked_files\dep1_info\', ...
%     'U:\moorings\IN2016_V06\moorings\EAC_stacked_files\dep2_info\'};

%%
[xfinal,yfinal,botd,mintime,maxtime] = deal(NaN*ones(length(mor),max(ndep)));
for b = 5:length(mor)
    for a = 1:ndep(b)
        % get positions and set indices in position vectors
        if b == 4
            ins = read_ins_info([depdir{3} 'instrument_info.csv']);
        elseif b > 4
            ins = read_ins_info([depdir{4} 'instrument_info.csv']);
        else
            ins = read_ins_info([depdir{a} 'instrument_info.csv']);
        end
        %get some information about the mooring:
        [um,m,n] = unique(ins.mooring,'rows');
%         if b > 4 %FOR ITF ONLY AT THIS STAGE
%             ii = strmatch(mn{b}{a},um,'exact');
%         else
         if ndep(b) > 1
            ii = strmatch(mn{b}{a},um,'exact');
         else
            ii = strmatch(mn{b},um,'exact');
         end
        if ~isempty(ii)
            xfinal(b,a) = ins.lon(m(ii));
            yfinal(b,a) = ins.lat(m(ii));
            botd(b,a) = ins.bot_depth(m(ii));
            
            % preset the timebase for all moorings to be put onto:
            %EAC
%             mintime = min(ins.start(m(1:5))); %first in
%             maxtime = max(ins.end(m(1:5))); %last out
            %ITF
            mintime(b,a) = min(ins.start(m(ii))); %first in
            maxtime(b,a) = max(ins.end(m(ii))); %last out
            %     mintime(a) = min(ins.start(m(1:3))); %last in
            %     maxtime(a) = max(ins.end(m(1:3))); %first out
        end
    end
end

dt = 1/24;
mtime = min(min(mintime)):dt:max(max(maxtime));
mtime = mtime(:);
idate = min(mintime(:,1)):1:max(max(maxtime));%ITF
%     idate = min(mintime(:,1)):1:max(maxtime(:,2));%ITF, just looking at the dates our deepwater moorings cover
%     idate = mintime:1:maxtime;%EAC
%%
 for im=5:length(mor)
    moorn = mor{im};

    xmoor = xfinal(im,:);
    ymoor = yfinal(im,:);

    % now put all deployments onto the common depth grid:
    if exist('maxd','var')
        di = [5:10:405,425:20:maxd(im)];
        ext = '_maxD';
    else
        if max(botd(im,:)) < 400
            di = [5:10:max(botd(im,:))+5];
        else
            di = [5:10:405,425:20:max(botd(im,:))];
        end
        ext = '_bottom';
    end
    %%
    for idep =1:ndep(im)

%         if ndep(im) < 2;
            load([inputdir,omn{im} 'mooring_dep',int2str(idep)]);
%         else
%             load([inputdir,mn{im}{idep} 'mooring_dep',int2str(idep)]);
%         end
        
        ui = NaN*ones(length(tbase),length(di))*[1+i];
        [ti,si ] = deal(NaN*ones(length(tbase),length(di)));
        
        
        [nt,nu]=size(u);
        [nt,ntp]=size(t);
        
        % now use to grid velocity data
        %       di = [0:5:300, 325:25:max([depu(:);dept(:)])];
        
        
        % ti,ui are interpolated values, tw and uw are weights to track
        % 'holes'
        
        
        for j=1:length(tbase);
            if any(depu(j,:) > 0);
                ig = find(~isnan(real(u(j,:))) & ~isnan(depu(j,:)) & depu(j,:) > 0 );
                
                if length(ig) > 2
                    if im > 4 %anmn moorings
                        ui(j,:) = akima(depu(j,ig)',u(j,ig)',di');
                    else
                        % put flow = 0 at the bottom:
                        idd = find(di < [botd(im)-60]);
                        ui(j,idd) = akima([depu(j,ig)';depu(j,ig(end))'+depoff(im);botd(im,idep)],[u(j,ig),u(j,ig(end)),[0+i*0]].',di(idd)');
                    end
                    %                     if j > 14200
%                     clf,plot(real(ui(j,:)),di,real(u(j,ig)),depu(j,ig),'rx'),axis ij,grid
%                     pause
%                     end
                end
            end
            
            if any(dept(j,:) > 5);
                
                ig = find(~isnan(t(j,:)) & ~isnan(dept(j,:)) );
                if length(ig) > 1
                    ti(j,:) = interp1q([dept(j,ig)';dept(j,ig(end))'+depoff(im);botd(im,idep)],[t(j,ig)';t(j,ig(end));t(j,ig(end))],di')';
%                     clf,plot(di,ti(j,:),'r-o');hold on;plot(dept(j,:),t(j,:),'x')
%                     pause
                end
            end
            
            %salinity
            if ~isempty(deps)
            if any(deps(j,:) > 5) ;
                
                ig = find(~isnan(sal(j,:)) & ~isnan(deps(j,:)) );
                if length(ig) > 0
                    si(j,:) = interp1q([deps(j,ig)';deps(j,ig(end))'+depoff(im);botd(im,idep)],[sal(j,ig)';sal(j,ig(end));sal(j,ig(end))],di')';
%                     clf,plot(di,ti(j,:),'r-o');hold on;plot(dept(j,:),t(j,:),'x')
%                     pause
                end
            end
            end
        end  % j loop
        
        %
        
        
        eval(['t' num2str(idep) ' = ti;'])
        eval(['u' num2str(idep) ' = ui;'])
        eval(['s' num2str(idep) ' = si;'])
        eval(['ndate' num2str(idep) ' = tbase;'])
    end % end deployment loop
    
    %% merge deployments and put on common   time grid:
    
    ui = NaN*ones(length(mtime),length(di))*[1+i];
    [ti,si ] = deal(NaN*ones(length(mtime),length(di)));
    ttim = [];uu = []; tt = [];ss = [];utim = [];stim = [];
    for idep = 1:ndep(im)
        eval(['ttim = [ttim; ndate' num2str(idep) '];'])
        eval(['utim = [utim; ndate' num2str(idep) '];'])
        eval(['stim = [stim; ndate' num2str(idep) '];'])
        eval(['uu = [uu; u' num2str(idep) '];'])
        eval(['tt = [tt; t' num2str(idep) '];'])
        eval(['ss = [ss; s' num2str(idep) '];'])
    end
    % get rid of no data times
    
    ibv = find(all(isnan(real(uu)')));
    utim = ttim;utim(ibv) = [];
    uu(ibv,:) = [];
    
    ibt = find(all(isnan(tt')));
    ttim(ibt) = [];
    tt(ibt,:) = [];
    
    ibt = find(all(isnan(ss')));
    stim(ibt) = [];
    ss(ibt,:) = [];

    for id =1:length(di)
        
        
        z = [uu(:,id)];
        ig = ~isnan(z);
        if sum(ig) > 3*30*24 % need 3 months of data
            ui(:,id)=interp1(utim(ig),z(ig),mtime(:));
            mask = 0.*ones(size(z));
            mask = interp1(utim ,double(ig) ,mtime(:));
            ib = find(mask < 0.5);
            ui(ib,id) = (1+i)*NaN;
            
        end
        z = [tt(:,id)];
        ig = ~isnan(z);
        if sum(ig) > 3*30*24 % need 3 months of data
            ti(:,id)=interp1(ttim(ig),z(ig),mtime(:));
            mask = interp1(ttim,double(ig),mtime(:));
            
            ib = find(mask < 0.5);
            ti(ib,id) = NaN;
        end
        
        z = [ss(:,id)];
        ig = ~isnan(z);
        if sum(ig) > 3*30*24 % need 3 months of data
            si(:,id)=interp1(stim(ig),z(ig),mtime(:));
            mask = interp1(stim,double(ig),mtime(:));
            
            ib = find(mask < 0.5);
            si(ib,id) = NaN;
        end
    end
    
    % Set Nan for any time periods outside the original deployment times (between deployments)
    if ndep(im) > 1;
        for indep=1:ndep(im)-1
            outofrange = mtime > maxtime(im,indep) & mtime < mintime(im,indep+1);
            ti(outofrange,:) = NaN;
            si(outofrange,:) = NaN;
            ui(outofrange,:) = NaN + NaN*i;
        end
    end
   
    
    save([inputdir,omn{im},'_vert_extrap' ext '.mat'],'di','mtime','ui','ti','si','moorn','xmoor','ymoor')
    
    
    %% do a spectra as function of depth:
    
    dospec = 0
    if dospec
        nsampfreq = 1./[diff(mtime(1:2))*24];
        % U,V and T
        
        for iv = 1:3
            switch iv
                case 1
                    zz = real(ui);
                case 2
                    zz = imag(ui);
                case 3
                    zz =ti;
            end
            
            mtime = mtime(:);
            % fill gaps and detrend
            for j=1:length(di);
                ig = (~isnan(zz(:,j))');
                if sum(ig) > length(mtime)*0.7
                    zz(ig,j) = detrend(zz(ig,j));
                    zz(~ig,j) = interp1q(mtime(ig), zz(ig,j),mtime(~ig));
                else
                    zz(:,j) = NaN*zz(:,j);
                end
            end
            
            
            
            idg = find(~all(isnan(zz)));
            ig = (all(~isnan(zz(:,idg))'));
            if sum(ig) > length(mtime)*0.7
                
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
                P = NaN*ones(length(mtime/2),length(di));
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
        
%         print('-dpng',[outputdir,omn{im},'_interp_spectra.png'])
        %
    end % dospec
    
    %% now low pass
    clear range
    
    % now filter the data at nf days: choose 5 days
    uif = [NaN+ i*NaN]*ones(length(mtime),length(di));
    tif = NaN*ones(length(mtime),length(di));
    sif = NaN*ones(length(mtime),length(di));
    nf = ceil(5/diff(mtime(1:2)));
    for j=1:length(di);
        ig = ~isnan(ui(:,j));
        if sum(ig) > 2*nf;
            uif(ig,j) = filt_ends(hamming(nf),ui(ig,j),1);
        end
        ig = ~isnan(ti(:,j));
        if sum(ig) > 2*nf;
            tif(ig,j) = filt_ends(hamming(nf),ti(ig,j),1);
        end
        ig = ~isnan(si(:,j));
        if sum(ig) > 2*nf;
            sif(ig,j) = filt_ends(hamming(nf),si(ig,j),1);
        end
    end
    
    
    % decimate: this is where we need to interpolate onto a common timebase
    idate = idate(:);
    
    nf = ceil(5/diff(mtime(1:2)));
    idec = 1:ceil(nf/2):length(mtime);
    
    uif = interp1q(mtime,uif,idate);
    
    tif = interp1q(mtime,tif,idate);
    
    sif = interp1q(mtime,sif,idate);


    save([inputdir,omn{im},'_daily_interp' ext '.mat'],'di','moorn','idate','uif','tif','sif','xmoor','ymoor')
    
    
    mtimef = idate;
    
    %idec = 1:length(mtimef);
    idec = ~all(isnan(real(uif')));
    %% plot interpolated data
    
    figure(4)
    clmap(23)
    clf
    subplot(311)
    contourf(mtimef(idec),di,real(uif(idec,:))'*100,[-100:10:100]),caxis([-50,50]),shading flat,
    
    hold on
%     contour(mtimef(idec),di,real(uif(idec,:))'*100,[-50:10:50],'w'),
    [c,h]=contour(mtimef(idec),di,real(uif(idec,:))'*100,[0:0.01:0.01],'w');
    axis ij,axis( [range(mtimef(idec)),0.,max(di)+20]),datetick('x','m','keeplimits')
    set(h,'linewidth',1,'linestyle','-'),
    title(['Low pass (5 days) U at ',moorn])
%     colorbar
    
    subplot(312)
    contourf(mtimef(idec),di,imag(uif(idec,:))'*100,[-100:10:100]),caxis([-40,40]),shading flat,
    hold on
%     contour(mtimef(idec),di,imag(uif(idec,:))'*100,[-50:10:50],'w'),
    [c,h]=contour(mtimef(idec),di,imag(uif(idec,:))'*100,[0:0.01:0.01],'w');
    axis ij,axis( [range(mtimef(idec)),0.,max(di)+20]),datetick('x','m','keeplimits')
    set(h,'linewidth',1,'linestyle','-'),
    title(['Low pass (5 days) V at ',moorn])
%     colorbar
    
    subplot(313)
    contourf(mtimef(idec),di,tif(idec,:)',[0:1:30]),caxis([3,30]),shading flat,
    hold on
    contour(mtimef(idec),di,tif(idec,:)',[0:5:30],'w'),
    [c,h]=contour(mtimef(idec),di,tif(idec,:)',[0:10:30],'w');
    axis ij,axis( [range(mtimef(idec)),0.,max(di)+20]),datetick('x','m','keeplimits')
    set(h,'linewidth',2),
    title(['Low pass (5 days) T at ',moorn])
%     colorbar
    
    orient landscape
    print('-djpeg',[outputdir omn{im} '/' omn{im} ,'_interp_lowpass' ext '.jpg'])
    
    
    % .........................................................................
    %% ................
    
end  % mooring loop

