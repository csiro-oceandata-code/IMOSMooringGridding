function s = clean_start_end(s)
%uses s.start_index and s.end_index to set
%values before and after these indexes to NaN
%Bec Cowley, 2 October, 2012


% load start_end_ind %get the start end indices

% ii = strmatch(s.serial,stend.serial);
% if isempty(ii)
%     disp(['No start end times for ' s.serial '!'])
%     return
% end

start_index = find(abs(s.time - datenum(s.time_in,'dd/mm/yyyy HH:MM:SS')) < s.time_int/2/60/60/24);
end_index = find(abs(s.time- datenum(s.time_out,'dd/mm/yyyy HH:MM:SS')) < s.time_int/2/60/60/24);
if isempty(end_index)
    %use last time stamp
    end_index = length(s.time);
end
if isempty(start_index) & isempty(end_index)
    error('Cannot find start/end indicies')
end

[m,n]=size(s.temperature);

fnms = fieldnames(s);

for a=1:length(fnms)
    %skip if it's a qc field or time field
    if strmatch('time',fnms{a})
        continue
    end
    if strfind(fnms{a},'_qc')
        continue
    end
    eval(['dat = s.' char(fnms(a)) ';'])
    
    try
        fn = fieldnames(dat);
        for c = 1:length(dat)
            for b = 1:length(fn)
                eval(['dat2 = dat(c).' char(fn(b)) ';'])
                [o,p] = size(dat2);
                if o~=m
                    continue
                end
                
                dat2(1:start_index-1,:) = NaN;
                dat2(end_index+1:end,:) = NaN;
                
                eval(['s.' char(fnms(a)) '(c).' char(fn(b)) ' = dat2;'])
            end
        end
    catch
        [o,p] = size(dat);
        
        
        if o~=m
            continue
        end
        
        dat(1:start_index-1,:) = NaN;
        dat(end_index+1:end,:) = NaN;
        
        eval(['s.' char(fnms(a)) ' = dat;'])
    end
end

end
