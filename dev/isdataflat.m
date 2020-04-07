function bool = isdataflat(data,threshold)
if iscolumn(data)
    filteredSignal = stdfilt(data,ones(5,1));
elseif isrow(data)
    filteredSignal = stdfilt(data',ones(5,1));
else
    error('must be vector input')
end
bool = max(filteredSignal)>threshold;
bool = ~bool;

end