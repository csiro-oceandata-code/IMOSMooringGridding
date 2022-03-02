function [tt,toff]= check_t_drift(tbase,t,namet,moorn);
%get rid of empty fields:
inan = isnan(nansum(t,1));
t(:,inan) = [];
namet(inan) = [];

[nt,nd]=size(t);
col = jet(nd);
toff = [];  clear cleg
for id =1:(nd-1);
    tdiff = t(:,id) - t(:,id+1);
    % make legend:
    cleg{id} = [char(namet(id)),' -- ',char(namet(id+1))];
    
    % use 10 day blocks to calculate drift
    k=0;
    for j=1:ceil([5/diff(tbase(1:2))]):length(tbase);
        k=k+1;
        ii = find(abs(tbase - tbase(j)) <= 5 & ~isnan(tdiff));
        toff(id,k)=nanmean( tdiff(ii));   % for pressure verions change to minimum
        tt(k) = tbase(j);
    end
    
end
inanm = find(~isnan(nanmean(toff)));

clf
hold on
for id=1:[nd-1];
    %h(id)=plot(tt,toff(id,:) -  toff(id,1),'-ko','markerfacecolor',col(id,:),'color',col(id,:))
    h(id)=plot(tt,toff(id,:) -  toff(id,inanm(1)),'-ko','markerfacecolor',col(id,:),'color',col(id,:));
end

legend(h,char(cleg),'location','best')
title([moorn ' Block averaged temperature differences between instruments'])
grid
datetick('x','m','keeplimits')

return
