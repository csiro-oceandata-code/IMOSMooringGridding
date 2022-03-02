%% plot vert_extrap_mooring.m output for checking
clear
%ITF
ndep = [2,2,2,6,6,6,5];
inputdir=['/home/cowley/work/moorings/ITF_stacked_files/'];
outputdir='/home/cowley/work/moorings/ITF_stacked_files/';
mor = {'Timor Sill','Timor North','Ombai','Timor South','JBG','MHB','FTB'};
mn = {'Timor Sill';'Timor North';'Ombai',
    {'TimorSth_1105';'TimorSth_1201';'TimorSth_1207';'TimorSth_1212';
    'TimorSth_1306';'TimorSth_1401'},
    {'JBG_1105';'JBG_1201';'JBG_1207';'JBG_1212';'JBG_1306';'JBG_1401'},
    {'MHB_1105';'MHB_1201';'MHB_1207';'MHB_1212';'MHB_1306';'MHB_1401'},
    {'FTB_1105';'FTB_1207';'FTB_1212';'FTB_1306';'FTB_1401'}};
omn = {'TimorSill','TimorNorth','Ombai','TimorSth','JBG','MHB','FTB'};
depdir = {'/home/cowley/work/moorings/s1205moor/',
    '/home/cowley/work/moorings/ITF_dep2/2012_recoveries/',
    '/home/cowley/work/moorings/ITFanmn_data/'};

%%
%load up each deployment stack and plot
for a = 1:length(mor) %for each mooring
    %prepare the figures
    figure(1);clf;hold on
    title(['Temperature ' mor{a}])
    xlabel('Time')
    ylabel('Temperature')
    figure(2);clf;hold on
    title(['Velocity ' mor{a}])
    xlabel('Time')
    ylabel('Velocity - U')

    for b = 1:ndep(a)
        load([inputdir,omn{a} 'mooring_dep',int2str(b)]);
        
        %plot the stacked data
        figure(1)
        plot(tbase,t)
        
        figure(2)
        plot(tbase,real(u(:,1:5))')
                
    end
    %overlay the vert_extrap output
        figure(1)
    load([inputdir,omn{a} '_daily_interp.mat'])
    plot(idate,tif,'k-')
    datetick('x','mmmyy')
    grid
    
    figure(2)
    plot(idate,real(uif(:,1:20))','k-')
    datetick('x','mmmyy')
    grid
    
    figure(1)
    load([inputdir,omn{a} '_vert_extrap.mat'])
    plot(mtime,ti,'g-')
        
    pause
end

%%
cc = jet(length(mor));
figure(3);clf;hold on
figure(4);clf;hold on
figure(5);clf;hold on
for a = 1:length(mor) %for each mooring
    load([inputdir,omn{a} '_daily_interp.mat'])
    %mean velocity profile
    mu = nanmean(real(uif));
    mv = nanmean(imag(uif));
    figure(3)
    plot(mu,di,'k','linewidth',2,'color',cc(a,:))
    figure(4)
    plot(xmoor,ymoor,'.','color',cc(a,:),'markersize',20)
    figure(5)
    plot(mv,di,'k','linewidth',2,'color',cc(a,:))    
end
figure(3)
xlabel('Mean U velocity (m/s)')
ylabel('Depth')
legend(mor,'location','sw')
axis ij
grid
print('-dpng','/home/cowley/work/moorings/ITF_stacked_files/plots/ITFmeanUvel.png')
ylim([0 500])
print('-dpng','/home/cowley/work/moorings/ITF_stacked_files/plots/ITFmeanUvel_zoom.png')

figure(4)
xlabel('Longitude')
ylabel('Latitude')
grid
coast('k')
axis equal
xlim([120 135])
ylim([-16 -4])
legend(mor,'location','ne')
print('-dpng','/home/cowley/work/moorings/ITF_stacked_files/plots/ITFlocs.png')

figure(5)
xlabel('Mean V velocity (m/s)')
ylabel('Depth')
legend(mor,'location','sw')
axis ij
grid
print('-dpng','/home/cowley/work/moorings/ITF_stacked_files/plots/ITFmeanVvel.png')
ylim([0 500])
print('-dpng','/home/cowley/work/moorings/ITF_stacked_files/plots/ITFmeanVvel_zoom.png')
