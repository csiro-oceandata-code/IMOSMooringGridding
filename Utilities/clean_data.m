function s = clean_data(s)
%Uses the QC flags available to set bad data to NaN
%Bec Cowley, 20 March, 2013


fnms = fieldnames(s);

for a=1:length(fnms)
    %don't put NaNs in the time field:
%     if ~isempty(strmatch('time',fnms(a)))
%         continue
%     end
    iqc = strfind(fnms{a},'_qc');
    if isempty(iqc)
        %no QC
        continue
    end
    %get the data to go with the QC flags
    str = fnms{a}(1:iqc-1);
    if isempty(strmatch(str,fnms,'exact')) 
        %qc field, but no data field
        continue
    end
    
    eval(['qc = s.' fnms{a} ';'])
 
    eval(['dat = s.' str ';'])
    
    ibad = qc >2 & qc <5;
    
    if ~isempty(strmatch('u',fnms(a)))
        dat(ibad) = NaN + NaN * i;
    else
        dat(ibad) = NaN;
    end
    
    eval(['s.' str ' = dat;'])
    
end
