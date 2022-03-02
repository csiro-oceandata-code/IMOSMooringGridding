function [tlag,maxcor]=find_timeoffset(timebase,t1,t2);
%  [tlag,maxcor]=find_timeoffset(t1,t2);
% look for time offset between two mooring series t1 and t2
% assume interpolated to same time base, 'timebase'
%
% output:  t2 lags t1 by 'tlag' and correlation maxcor

% SEW, CSIRO MAR  March 2006

dt= nanmean(diff(timebase));
%find number of units equal to a day - Bec Cowley August 2007
% tdiff = diff(timebase);
% it = find(tdiff>0);
% hh = str2num(datestr(nanmean(tdiff(it)),'HH'));
% mm = str2num(datestr(nanmean(tdiff(it)),'MM'));
% ss = str2num(datestr(nanmean(tdiff(it)),'SS'));
% sampleRate = ss/60 + hh*60 + mm; %in minutes
% dy=round(12/(sampleRate/60))-5;
dy = 20; %use 20 lags.

% linearly fill gaps
ib = isnan(t1);
t1(ib) = interp1(timebase(~ib),t1(~ib),timebase(ib));
ib = isnan(t2);
t2(ib) = interp1(timebase(~ib),t2(~ib),timebase(ib));

ig = ~isnan(t1.*t2);
[cc,ll]=xcorr(normal(t1(ig)),normal(t2(ig)),dy,'coeff'); % search 'dy' time intervals for large shifts BecCowley, Aug 07.

plot(ll,cc,'-x'),grid,xlabel('lags')
hold on
grid
[maxcor,im]=max(cc);
% do akima to find max corr for subhour:
in = ([-2:2] + im);
try
    pp = polyfit(ll(in)',cc(in),3);
catch
    pp = polyfit(ll(in),cc(in),3);
end
lp = min(ll(in)):0.1:max(ll(in));
cp=polyval(pp,lp);
plot(lp,polyval(pp,lp),'r')

hold off

[maxcor,im]=max(cp);
tlag = lp(im)*dt;

 return
 
        