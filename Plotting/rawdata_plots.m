[nt,nd]=size(u);

  %%% cumsum plot
 jp = [1:nd];

 clear range
 clf
subplot(2,1,1)
pcolor(tbase,depu',real(u)');caxis([-1 1]);shading flat;colorbar
hold on
% plot(tbase,depu,'k')
datetick('keeplimits');axis ij
ylabel('Bin')
xlabel('Time')
title(['Zonal velocity after stacking - ' moorn])
caxis([-0.5 0.5])
subplot(2,1,2)
pcolor(tbase,depu',imag(u)');caxis([-.5,.5]);shading flat;colorbar
hold on
% plot(tbase,depu,'k')
caxis([-1 1])
datetick('keeplimits');axis ij
title(['Meridional velocity after stacking - ' moorn])
ylabel('Depth')
xlabel('Time')
eval(['print -dpng ',[outputdir dirn],'_veloc.png'])
% pause

 clf
pcolor(tbase,dept',t');caxis([10 25]);shading flat;colorbar
datetick('keeplimits')
axis ij
title(['Temperature after stacking - ' moorn])
ylabel('Depth')
xlabel('Time')
eval(['print -dpng ',[outputdir dirn],'_Temp.png'])


 clf
%  subplot(121)
%  h=plot(cumsum(repnan(u(:,jp),0.)) - diff(range(tbase))*i*ones(size(u(:,jp)))*diag(nanmean(depu(:,jp)))/100);
% %  legend(h,nameu(jp),-1)
%  title(['Prog vector plot - ' moorn]);% axis equal

%  subplot(122)
 plot(t,dept,'o')
 axis ij
 ylabel('Depth (m)')
 xlabel('temperature')
 orient portrait
 grid
 % print -dpng timor_sslope_raw.png
eval(['print -dpng ',[outputdir dirn],'_Temperature.png'])

figure(1);clf
 plot(sal,deps,'o')
 axis ij
 ylabel('Depth (m)')
 xlabel('salinity')
 orient portrait
 grid
 % print -dpng timor_sslope_raw.png
eval(['print -dpng ',[outputdir dirn],'_sal_dep.png'])
% pause

clf
h=plot(tbase,real(u) + ones(length(tbase),1)*[0:(nd-1)]/3);
%  legend(nameu,'location','best')
 title([moorn ' - U waterfall plot'])
 xlabel('Time')
 ylabel('No value - U')
datetick('x','m','keeplimits')
orient tall
eval(['print -dpng ',[outputdir dirn],'_U_time_raw.png'])

figure(1)
clf
hold on ,
h=plot(tbase,imag(u) + ones(length(tbase),1)*[0:(nd-1)]/3);
%  legend(nameu,'location','best')
 ylabel('No value - V')
 xlabel('Time')
 title([moorn ' - V waterfall plot'])
datetick('x','m','keeplimits')
orient tall
eval(['print -dpng ',[outputdir dirn],'_V_time_raw.png'])
% pause
figure(1)
clf
h=plot(tbase,dept );
axis ij
 legend(h,named,'location','eo')
 title([moorn '- Depth with time'])
 xlabel('Time')
 ylabel('Depth')
datetick('x','m','keeplimits')
orient landscape
eval(['print -dpng ',[outputdir dirn],'_p_time_raw.png'])

figure(1)
clf
h=plot(tbase,sal );
axis ij
 legend(h,names,'location','eo')
 title([moorn '- Salinity'])
 xlabel('Time')
 ylabel('Salinity')
datetick('x','m','keeplimits')
orient landscape
eval(['print -dpng ',[outputdir dirn],'_s_time_raw.png'])
% pause

clf
h=plot(tbase,t);
 legend(h,namet,'location','eo')
title([moorn '- T'])
datetick('x','m','keeplimits')
ylabel('temperature')
 xlabel('Time')
orient landscape
eval(['print -dpng ',[outputdir dirn],'_t_time_raw.png'])


clf
plot(real(u),depu,'.','markersize',2),axis ij
text(0.3*ones(1,length(nameu)),nanmean(double(depu)),nameu)
title([moorn ' - U'])
xlabel('U')
 ylabel('Mean Depth')
grid
orient tall
eval(['print -dpng ',[outputdir dirn],'_U_p_raw.png'])



clf
plot(imag(u),depu,'.','markersize',2),axis ij
text(0.3*ones(1,length(nameu)),nanmean(double(depu)),nameu)
title([moorn ' - V'])
xlabel('V')
 ylabel('Depth')
orient tall
grid
eval(['print -dpng ',[outputdir dirn],'_V_p_raw.png'])
% pause
%t-s plots
clf;hold on
for a = 1:length(names)
    ii = strmatch(names{a},namet);
    if ~isempty(ii)
        plot(sal(:,a),t(:,ii),'.')
    end
end
grid
xlabel('Salinity')
ylabel('Temperature')
legend(names,'location','best')
% pause
eval(['print -dpng ',[outputdir dirn],'t-s.png'])

% %check some of the temperatures for t-offset correction:
% for a = 1:1000:length(tbase)
% figure(4);clf
%     tof = tmp_off ~= 0;
%     plot(t(a,:),dept(a,:),'x-')
%     hold on
%     plot(t(a,(tof)),dept(a,(tof)),'ro')
%     axis ij
%     grid
%     pause
% end
% 
%% comparitive plots for more than one instrument on a mooring.
% check_angles

% [tt,poff]= check_p_drift(tbase,tdepth,snt,tdepth_m,moorn);
% orient landscape
% eval(['print -dpng ',[outputdir moorn],'_P_drift.png'])
% 
% [tt,toff]= check_t_drift(tbase,t,snt,moorn);
% orient landscape
% eval(['print -dpng ',[outputdir moorn],'_T_drift.png'])
% 
% plot_u_comparison
% orient landscape
% eval(['print -dpng ',[outputdir moorn],'_Ucomp.png'])
% 
% plot_v_comparison
% orient landscape
% eval(['print -dpng ',[outputdir moorn],'_Vcomp.png'])

