function [xm,xstd,xmed,nx] = nanmore(x);
% function [xm,xstd,xmed,nx] = nanmore(x);
% finds the mean, etc of vector x excluding the nan's
%
%
% S Wijffels 
dum = x;
[junk,nn] = size(x);
xm=zeros(1,nn);
xstd=xm; xmed=xm;nx=xm;
igood = ~isnan(dum);
for j=1:nn;
ii = find(igood(:,j));
if ~isempty(ii),
xm(j) = mean(dum(ii,j));
xmed(j) = median(dum(ii,j));
nx(j) = sum(igood(:,j));
xstd(j) = std(dum(ii,j));
end
xstd = xstd(:);
xmed = xmed(:);
xm = xm(:);
nx = nx(:);
end
