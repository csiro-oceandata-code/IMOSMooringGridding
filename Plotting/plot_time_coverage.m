% Load up all the stacked files and highlight the time gaps
clear all

%list the deployment names first:
depn = {'EAC1505_1611','EAC1611_1805','EAC1805_1909'}; %'EAC1204_1308',

%list the sites:
siten = {'EAC0500','EAC2000','EAC3200','EAC4200','EAC4700','EAC4800'};
leg = {};
%color by deployment:
cc = {'b','k','r','g','m','c'};

% now get the plots ready. Need one for each site, each type = 6*3 = 18
% figures
for b = 5:9
    figure(b);clf;hold on
end
fign = 5;

for a = 6%1:length(siten)
    leg1 = cell(1,length(depn));leg2 = leg1;leg3 = leg1;
    for idep = 1:length(depn)
        % Start in the  stacked directory
        d = ['/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/' depn{idep} '/stacked/'];
        fn = [d siten{a} 'mooring_allbins.mat'];
        
        if ~exist(fn,'file')
            fn = [d siten{a} 'mooring.mat'];
            if ~exist(fn,'file')
                continue
            end
        end
        
        %load the file
        load(fn)
        
        %will interpolate to daily
        tday = tbase(1):tbase(end);
        
        %temperature
        d = dept;
        d(isnan(t)) = NaN;
        tempd = interp1(tbase,d,tday);
        figure(fign)
        h1 = plot(tday,tempd,'-','color',cc{idep});
        leg1{idep} = h1(1);
        
        %salinity
        d = deps;
        d(isnan(sal)) = NaN;
        sald = interp1(tbase,d,tday);
        figure(fign+1)
        h2 = plot(tday,sald,'-','color',cc{idep});
        leg2{idep} = h2(1);
        
        %currents
        d = depu;
        d(isnan(u)) = NaN;
        ud = interp1(tbase,d,tday);
        figure(fign+2)
        h3 = plot(tday,ud,'-','color',cc{idep});
        leg3{idep} = h3(1);
    end
    
    %clean up the legend
    hh1 = {};leglabel1=[];
    for x = 1:length(leg1)
        if ~isempty(leg1{x})
            hh1 = [hh1,leg1{x}(1)];
            leglabel1 = [leglabel1,depn(x)];
        end
    end
    %clean up the legend
    hh2 = {};leglabel2=[];
    for x = 1:length(leg2)
        if ~isempty(leg2{x})
            hh2 = [hh2,leg2{x}(1)];
            leglabel2 = [leglabel2,depn(x)];
        end
    end
    %clean up the legend
    hh3 = {};leglabel3=[];
    for x = 1:length(leg3)
        if ~isempty(leg3{x})
            hh3 = [hh3,leg3{x}(1)];
            leglabel3 = [leglabel3,depn(x)];
        end
    end
    
    %make the plots pretty
    figure(fign);axis ij;grid;datetick('x','mm/yy')
    title([siten{a} ' Temperature data coverage'])
    xlabel('Time')
    ylabel('Depth')
    l=legend(hh1,leglabel1,'location','southeast');
    set(l,'Interpreter','none')
    
    figure(fign+1);axis ij;grid;datetick('x','mm/yy')
    title([siten{a} ' Salinity data coverage'])
    xlabel('Time')
    ylabel('Depth')
    legend(siten)
    l=legend(hh2,leglabel2,'location','se');
    set(l,'Interpreter','none')
    
    figure(fign+2);axis ij;grid;datetick('x','mm/yy')
    title([siten{a} ' Currents data coverage'])
    xlabel('Time')
    ylabel('Depth')
    l=legend(hh3,leglabel3,'location','se');
    set(l,'Interpreter','none')
    
    fign = fign + 3;

end
