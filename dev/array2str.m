function cellstr = array2str(cellstr)
if ~iscellstr(cellstr)
    [nrows,ncols] = size(cellstr);
for row = 1:nrows
    for col = 1:ncols
        if isnumeric(cellstr{row,col})
            cellstr{row,col} = num2str(cellstr{row,col});
        elseif isa(cellstr{row,col},'function_handle')
            cellstr{row,col} = func2str(cellstr{row,col});
        end
    end
end
end
end