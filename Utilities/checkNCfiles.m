%% Check exported netcdf files make some sort of sense!

inpd = '/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC_joined/';
cd(inpd)
fils = dir([inpd 'EAC*daily-depth*-2022*.nc']);
var = {'TEMP','PSAL','UCUR','VCUR'};
for b = 3%:length(var)
    for a = 1:length(fils)
        dat = nc2struct(fils(a).name);
        figure(1);clf
        subplot(2,1,1)
        pcolor(dat.TIME,dat.DEPTH,dat.(var{b}));shading flat; axis ij;colorbar
        title(var{b})
        subplot(2,1,2)
        pcolor(dat.TIME,dat.DEPTH,dat.([var{b} '_uncertainty']));shading flat; axis ij;colorbar
        title(fils(a).name)
        figure(2);clf
        subplot(2,1,1)
        pcolor(dat.TIME,dat.DEPTH,dat.([var{b} '_FILLED']));shading flat; axis ij;colorbar
        title([var{b} ' FILLED'])
        subplot(2,1,2)
        pcolor(dat.TIME,dat.DEPTH,dat.([var{b} '_FILLED_uncertainty']));shading flat; axis ij;colorbar
        title(fils(a).name)
        figure(3);clf
        pcolor(dat.TIME,dat.DEPTH,dat.(var{b}) ./ dat.([var{b} '_uncertainty']));shading flat; axis ij;colorbar
        caxis([-0.01 0.01])
        pause
    end
end
    
