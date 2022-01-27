function bool = isdataflat(data,threshold,n_points)
if nargin<3
    n_points = 5;
end
if iscolumn(data)
    filteredSignal = stdfilt(data,ones(n_points,1));
elseif isrow(data)
    filteredSignal = stdfilt(data',ones(n_points,1));
else
    error('must be vector input')
end
bool = (max(filteredSignal)/mean(filteredSignal))<threshold;

end