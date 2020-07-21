read_paths={};

%if we have multiple paths we can steup a looping structure as follows
%paths = {all the different paths}
%for j = 1:length(paths)
%then replace path with paths{j}
if i == 1
    exact_line =zeros(1,numel(paths));
    comment_line = zeros(1,numel(read_paths));
end
for k = 1:numel(paths)
    %bryce: what happens when the path string cant be found? is there an
    %error message?
    
    %select the path we want to change
    path = paths{k};
    
    %we only wish to find the exact path once, as it is computationally
    %expensive
    %if i==1
    if true %hack to see if finding each time helps
        cline = 1;
        for j = 1:length(path)-1
            count = 0;
            toggle = 1; 
            test = 0;
            if iscell(path{j})
                name = path{j}{1};
                num = path{j}{2};
                end_name = strrep(name,'<','</');
            elseif ischar(path{j})
                name = path{j};
                num = 1;
                end_name = '';
            end
            while count<num
                cline = cline + 1;
                if ~isempty(strfind(A{1,cline},name))
                    if test == 0 && count>0
                        toggle = 1 - toggle;
                    end
                    if toggle
                        count = count + 1;
                        toggle = 1 - toggle;
                    end
                    test = test + 1;
                elseif ~isempty(strfind(A{1,cline},end_name))
                    test = test - 1;
                end
            end
        end
        exact_line(1,k) = cline+1;
    end

    %set the value to what we want it to be
    start_ps=strfind(A{exact_line(1,k)},'>');
    end_ps= strfind(A{exact_line(1,k)},'<');
    strt = A{exact_line(1,k)}(1:start_ps(1));
    fin = A{exact_line(1,k)}(end_ps(2):end);
    A{exact_line(1,k)} = strcat(strt,path{end},fin);
end

if i==1
    for k = 1:numel(read_paths)
        %bryce: what happens when the path string cant be found? is there an
        %error message?
        %select the path we want to change
        path = paths{k};
        %we only wish to find the exact path once, as it is computationally
        %expensive

        cline = 1;
        for j = 1:length(path)-1
            count = 0;
            toggle = 1; 
            test = 0;
            if iscell(path{j})
                name = path{j}{1};
                num = path{j}{2};
                end_name = strrep(name,'<','</');
            elseif ischar(path{j})
                name = path{j};
                num = 1;
                end_name = '';
            end
            while count<num
                cline = cline + 1;
                if ~isempty(strfind(A{1,cline},name))
                    if test == 0 && count>0
                        toggle = 1 - toggle;
                    end
                    if toggle
                        count = count + 1;
                        toggle = 1 - toggle;
                    end
                    test = test + 1;
                elseif ~isempty(strfind(A{1,cline},end_name))
                    test = test - 1;
                end
            end
        end
        comment_line(1,k) = cline+1;

        %or use the exact line is already known
        % exact_line = {};

        %print out what the comment says
        j=0;
        end_of_comment = 1;
        while 0
            %display this line somehow
            %A{comment_line(1,k)+j}
            if isempty(strfind(A{comment_line(1,k)+j},'</Val>'))
                j = j+1;
            else
                end_of_comment = 0; 
            end
        end
    end
end
%% write out new file
file_id = fopen(new_path,'w');
for j = 1:numel(A)
    if A{j+1} == -1
        fprintf(file_id,'%s', A{j});
        break
    else
        fprintf(file_id,'%s\n', A{j});
    end
end
fclose(file_id);