function check_angles(tbase,u,depu,nameu,moorn,outputdir,dirn)
% check velocity offsets

%  input : u, nameu - velocity data in complex vector form on same time base
% output: plots and estimates of amplitude and angle offsets
% assumes velocities have been depu sorted.

%first, let's find adjacent bins WITH DATA from different instruments
snan = sum(~isnan(u))/length(tbase);
%remove bins with <25% data
irem = snan < 0.025;
u(:,irem) = [];
nameu(irem) = [];
depu(:,irem)=[];
%find adjacent instruments
[vel_name,ia,ib] = unique(nameu,'stable');

% the bins we are looking at are ia(2:end) and one bin before
% but, need to make sure there is data in these bins, if not, extend by 1
% if possible

k=0; rot = [];amp = [];dm = [];in = [];ampd = [];rotd=[];
ig = ~isnan(tbase);
tbase = tbase(ig);
u = u(ig,:);
depu = depu(ig,:);
for iz = 2:length(ia)
    k = k+1;
    data1 = u(:,ia(iz)); % deeper instrument
    data2 = u(:,ia(iz)-1); % shallower instrument
    
    %look at coherence between the adjacent bins for speed and direction:
    r = data2./data1;
    figure(1); clf
    subplot(121); hold on
    ig = find (abs(data1) > 0.05);
    h=histogram(abs(r(ig)),72,'normalization','pdf');
    ylabel('Count')
    xlabel('Speed Coherence, m/s')
    title([char(nameu(ia(iz)-1)),' divided by ',char(nameu(ia(iz)))])
    subplot(122); hold on
    h2 =histogram(angle(r(ig)),72,'normalization','pdf');
    ylabel('Count')
    xlabel('Direction Coherence, degrees')
    
    %get the most common speed
    cc = h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
    [val,ii]=nanmax(h.Values);
    amp(k) =  cc(ii);
    %get the most common direction
    cc = h2.BinEdges(1:end-1)+diff(h2.BinEdges)/2;
    [val,ii]=nanmax(h2.Values);
    rot(k) = cc(ii)*180/pi; %now in radians?
    
    %mean depth
    dm(k) = nanmean(depu(:,ia(iz)));
    %keeping track of index
    in(k) = ia(iz)-1;
    
    % difference in time - tidal step
    ntide = floor( 6/[diff(tbase(1:2))*24]);
    i1 = 1:[length(tbase)-ntide+1];
    i2 = ntide:length(tbase);
    ut = data2(i1,:)-data2(i2,:);
    ub =  data1(i1,:)-data1(i2,:);
    r = ut./ub; %coherence with 2 hour tidal step offset
    
    subplot(121)
    ig = find(abs(ub) > 0.08);
    [bb,cc] = histcounts(abs(r(ig)),72,'normalization','pdf');
    %keep the maximum with the tidal comparison
    [val,ii]=nanmax(bb);
    ampd(k) = cc(ii);
    cc = cc(1:end-1)+diff(cc)/2;
    plot(cc,(bb/sum(bb)),'r')
    grid
    subplot(122)
    [bb,cc] =histcounts(angle(r(ig)),72,'normalization','pdf');
    %keep the maximum with the tidal comparison
    [val,ii]=nanmax(bb);
    rotd(k) =  cc(ii)*180/pi;
    cc = cc(1:end-1)+diff(cc)/2;
    plot(cc,bb/sum(bb),'r')
    grid
    legend('as is','tidal step')
    pause(0.5)
    
    pname = [moorn,'_',char(nameu(iz-1)),'_over_',char(nameu(iz))];
    is = findstr(' ',pname);
    pname(is) = '_';
    orient portrait
    eval( ['print -djpeg ' outputdir 'check_angles_',pname,'.jpg'])
    
end

figure(2),
clf
plot(dm,amp*10,'b*',dm,rot,'r*-')
hold on
h=plot(dm,ampd*10,'bo',dm,rotd,'ro');set(h,'markersize',10)
hold on
xlabel('mean depth (m)')

for a = 1:length(dm)
    text(dm(a),ampd(a)*10,[nameu(in(a)),' speed',num2str(ampd(a))])
    text(dm(a),rotd(a),[' dir ' ,num2str(rotd(a))])
    text(dm(a),amp(a)*10,[nameu(in(a)),' speed',num2str(amp(a))])
    text(dm(a),rot(a),[' dir ',num2str(rot(a))])
end
grid
orient portrait
legend('speed','dir','speed tidal','dir tidal','location','best')
eval( ['print -djpeg ' outputdir 'check_angles_sum',dirn,'.jpg'])

end