function [ss,soff]= check_s_drift(tbase,s,names,moorn);
%get rid of empty fields:
inan = isnan(nansum(s,1));
s(:,inan) = [];
names(inan) = [];

[nt,nd]=size(s);
col = jet(nd);
soff = [];  clear cleg
for id =1:(nd-1);
    sdiff = s(:,id) - s(:,id+1);
    % make legend:
    cleg{id} = [char(names(id)),' -- ',char(names(id+1))];
    
    % use 10 day blocks to calculate drift
    k=0;
    for j=1:ceil([10/diff(tbase(1:2))]):length(tbase);
        k=k+1;
        ii = find(abs(tbase - tbase(j)) <= 10 & ~isnan(sdiff));
        soff(id,k)=nanmean( sdiff(ii));   % for pressure verions change to minimum
        ss(k) = tbase(j);
    end
    
end
inanm = find(~isnan(nanmean(soff)));

clf
hold on
for id=1:[nd-1];
    %h(id)=plot(tt,toff(id,:) -  toff(id,1),'-ko','markerfacecolor',col(id,:),'color',col(id,:))
    h(id)=plot(ss,soff(id,:) -  soff(id,inanm(1)),'-ko','markerfacecolor',col(id,:),'color',col(id,:));
end

legend(h,char(cleg),'location','best')
title([moorn ' Block averaged salinity differences between instruments'])
grid
datetick('x','m','keeplimits')

return
