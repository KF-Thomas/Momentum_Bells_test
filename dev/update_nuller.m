function update_nuller(n_0,n_1,n_2,n_3,i)
% file = 'C:\Remote\settings202021Jul105753.xml';
file = 'c:\remote\settings202023Jul135844.xml';

%% read in the xml settings file
new_path = 'c:\remote\settings202023Jul120059.xml';%file;
file_id = fopen(file,'r');
j=1;
line = fgetl(file_id);
A{j} = line;
while ischar(line)
    j = j+1;
    line = fgetl(file_id);
    A{j} = line;
end
fclose(file_id);

% delete(file);   % DEBUG - delete xml file
%% make adjustments

nuller_0 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',4},'Final value, Amplitude (exp)',num2str(-n_0)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',5},'Initial Value ',num2str(n_0)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',6},'Initial Value ',num2str(n_0)}};
nuller_1 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',10},'Final value, Amplitude (exp)',num2str(-n_1)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',11},'Initial Value ',num2str(n_1)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',12},'Initial Value ',num2str(n_1)}};
nuller_2 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',16},'Final value, Amplitude (exp)',num2str(-n_2)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',17},'Initial Value ',num2str(n_2)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',18},'Initial Value ',num2str(n_2)}};
nuller_3 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',40},'Final value, Amplitude (exp)',num2str(-n_3)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',41},'Initial Value ',num2str(n_3)},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',42},'Initial Value ',num2str(n_3)}};

paths = {
    nuller_0{:},...
    nuller_1{:},...
    nuller_2{:},...
    nuller_3{:}
    };  


A = xml_edit(A,paths,i);

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

%%clean up
clear('A')
end
% 
% f_log=fopen(path_log,'a');  % append to log-file
% fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : finished. \n']);
% fclose(f_log);

% copyfile(new_path,'C:\Users\BEC Machine\Dropbox\debug\new_path_cp');    % DEBUG