function y = afiltfilt(b,a,x)
%FILTFILT Zero-phase forward and reverse digital filtering.
%	Y = FILTFILT(B, A, X) filters the data in vector X with the
%	filter described by vectors A and B to create the filtered
%	data Y.  The filter is described by the difference equation:
%
%	  y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%	                   - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%	After filtering in the forward direction, the filtered
%	sequence is then reversed and run back through the filter.
%	The resulting sequence has precisely zero-phase distortion
%	and double the filter order.  Care is taken to minimize
%	startup and ending transients by matching initial conditions.
%
%	The length of the input x must be more than three times
%	the filter order, defined as max(length(b)-1,length(a)-1).
%
%	See also FILTER.

%	Author(s): L. Shure, 5-17-88
%	revised by T. Krauss, 1-21-94
%	initial conditions: Fredrik Gustafsson
%	Copyright (c) 1984-94 by The MathWorks, Inc.
%	$Revision: 1.9 $  $Date: 1994/01/25 17:59:07 $

    error(nargchk(3,3,nargin))
    if (isempty(b)|isempty(a)|isempty(x))
        y = [];
        return
    end

    [m,n] = size(x);
    if (n>1)&(m>1)
        error('Only works for vector input.')
    end
    if m==1
        x = x(:);   % convert row to column
    end
    len = size(x,1);   % length of input
    b = b(:).';
    a = a(:).';
    nb = length(b);
    na = length(a);
    nfilt = max(nb,na);

    nfact = 3*(nfilt-1);  % length of edge transients

    if (len<=nfact),    % input data too short!
        error('Data must have length more than 3 times filter order.');
    end

% set up filter's initial conditions to remove dc offset problems at the 
% beginning and end of the sequence
    if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
    if na < nfilt, a(nfilt)=0; end
  zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
       ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at
% the end points.
% This reduces end effects.
    y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];

% filter, reverse data, filter again, and reverse data again
    y = filter(b,a,y,[zi*y(1)]);
    y = y(length(y):-1:1);
    y = filter(b,a,y,[zi*y(1)]);
    y = y(length(y):-1:1);

% remove extrapolated pieces of y
    y([1:nfact len+nfact+(1:nfact)]) = [];

    if m == 1
        y = y.';   % convert back to row if necessary
    end

