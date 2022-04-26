function [ind_first, ind_secon, npairs] = get_pairs(tdepth_m)

[aa,bb]= size(tdepth_m);

ii = 1;
ind_1 = 0;
ind_2 = 0;

npairs = 0;
ind_first = [];
ind_secon = [];
while ii < bb +1
    
   if(tdepth_m(ii) == 1)
       
       if (ind_1) == 0
           ind_first(1) = ii;
           ind_1 = 1;
       else
           ind_2 = ind_2 + 1;
           ind_secon(ind_2) = ii;
           npairs = npairs + 1;
           ind_1 = ind_1 + 1;
           ind_first(ind_1) = ii;
       end
       
   end
    
   ii = ii + 1;
end




