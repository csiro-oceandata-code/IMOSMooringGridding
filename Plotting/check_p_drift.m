function [tt,poff]= check_p_drift(tbase,dept,named,tdepth_m,moorn);

[ind_first, ind_secon, npairs ] = get_pairs(tdepth_m);

[nt,nd]=size(dept);

col = jet(npairs);

poff = [];  clear cleg
for ii =1 : npairs;
    id = ind_first(ii);
    id_2 = ind_secon(ii);
    pdiff = dept(:,id) - dept(:,id_2);
    % make legend:
    
    %cleg{id} = [char(named(id)),' -- ',char(named(id_2))];
    cleg{ii} = [char(named(id)),' -- ',char(named(id_2))];
    % use 10 day blocks to calculate drift
    k=0;
    
    for j=1:ceil(5/diff(tbase(1:2))):length(tbase)
        
        k=k+1;
        %ii = find(abs(tbase - tbase(j)) <= 5 & ~isnan(pdiff));
        jj = find(abs(tbase - tbase(j)) <= 5);
        
        poff(id,k) = nanmean( pdiff(jj));
        tt(k) = tbase(j);
    end
    
end
inanm = find(~isnan(nanmean(poff)));

clf
hold on
for ii=1:npairs;
    id = ind_first(ii);
    id_2 = ind_secon(ii);
    h(ii)=plot(tt,poff(id,:) - poff(id,inanm(1)),'-ko','markerfacecolor',col(ii,:),'color',col(ii,:));
%     disp(cleg(ii))
%     pause
end

try
legend(char(cleg),'location','best')
catch
tt=[];poff=[];
end
%legend(char(cleg{2}),char(cleg{3}),char(cleg{5}),char(cleg{7}),char(cleg{9}))
title([moorn ' Block averaged pressure differences between instruments'])
grid
datetick('x','m','keeplimits')

return
