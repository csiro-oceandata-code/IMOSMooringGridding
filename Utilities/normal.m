function [y,ymean,ystd] = normal(x)
%	function [y,ymean,ystd] = normal(x)
%  demeans and normalizes the colums of matrix x 
% Returns the mean, std 
%  and the demeaned data in y.
%
[N,M] = size(x);
ymean = mean(x);
y = x - ones(N,1)*ymean;
ystd = std(x);
y = y./(ones(N,1)*ystd);

