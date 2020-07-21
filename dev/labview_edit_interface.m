%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DOCUMENTATION
%bryce:what needs to be in the path for this to run properly
%explain briefly what this program does and what the user will need to edit

%to do
%could read in number of itterations from settings file in order to be a
%bit intelegent

%INPUT FROM LABVIEW
%mloop,boolean : do you want to interface with mloop or use matlab to scan over variables
%file,string : tells program where to look for setting files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



path_user_config='..\MloopUserConfig.txt';  % this is the user config file regardless of Mloop
user_config_fp=fopen(path_user_config,'r');  % read config file

while true
  this_line = fgetl(user_config_fp);
  if ~ischar(this_line)
      break;
  end  %end of file
  eval(this_line);
end


fclose(user_config_fp);

file = ;

%save(['.\logs\interfacev5',datestr(datetime,'yyyymmdd_HHMMSS'),'.mat'])

%List of all paths of control variables


% %for IT cool opt
% paths = {
%     {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',4},'Initial Value ',''}},...%IT AMP 1
%     {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',5},'Initial Value ',''}},...%IT AMP 2
%     {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',6},'Initial Value ',''}},...%IT AMP 3
%     {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',21},'Initial Value ',''}},...%IT freq 1
%     {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',22},'Initial Value ',''}},...%IT freq 2
%     {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6289M (Ch0...3)',{'<Cluster>',23},'Initial Value ',''}},...%IT freq 3
%     {{{'<Cluster>',1},{'<Cluster>',5},'Channel 0 - Power (dBm), Channel 1 - Frequency (MHz)',{'<Cluster>',14},'Initial Value (dBm/MHz) ',''}},...%evap freq 1
%     {{{'<Cluster>',1},{'<Cluster>',5},'Channel 0 - Power (dBm), Channel 1 - Frequency (MHz)',{'<Cluster>',14},'Final value, Amplitude (exp)',''},...%evap freq 2
%     {{'<Cluster>',1},{'<Cluster>',5},'Channel 0 - Power (dBm), Channel 1 - Frequency (MHz)',{'<Cluster>',15},'Initial Value (dBm/MHz) ',''}},...
%     {{{'<Cluster>',1},{'<Cluster>',5},'Channel 0 - Power (dBm), Channel 1 - Frequency (MHz)',{'<Cluster>',15},'Final value, Amplitude (exp)',''},...%evap freq 3
%     {{'<Cluster>',1},{'<Cluster>',5},'Channel 0 - Power (dBm), Channel 1 - Frequency (MHz)',{'<Cluster>',16},'Initial Value (dBm/MHz) ',''}}
%     };      % 14-param shunting profile
% 

nuller_0 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',4},'Final value, Amplitude (exp)',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',5},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',6},'Initial Value ',''}};
nuller_1 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',10},'Final value, Amplitude (exp)',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',11},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',12},'Initial Value ',''}};
nuller_2 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',16},'Final value, Amplitude (exp)',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',17},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',18},'Initial Value ',''}};
nuller_3 = {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',40},'Final value, Amplitude (exp)',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',41},'Initial Value ',''},...
    {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',42},'Initial Value ',''}};

paths = {
    nuller_0,...
    nuller_1,...
    nuller_2,...
    nuller_3
    };  

%%%%NOTE: QUAD FIRST THEN SHUNT!!!! %%%%
%use a path of names to change variables, with the value you want to change
%the variable too at the end of the path
%eg path = {path to variable, value to change too};
%can also conjoin variables so they are set to the same parameter, simply
%put their paths into a cell array in athe paths master list

% IT cool opt
%volt_lim = [0,10];
%freq_lim=[0.7,22];
%param_limits=[repmat(volt_lim,6,1);repmat(freq_lim,3,1)];    % NOTE: MATLAB's hard-coded parameter limits

%% read in the xml settings file
new_path = file;
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

param = 

conjoined_indx = [];

for idx=1:numel(paths)
    %check if param values are within the hard coded safety range
%     if param(idx)>param_limits(idx,2)
%         param(idx) = param_limits(idx,2);
%     elseif param(idx)<param_limits(idx,1)
%         param(idx) = param_limits(idx,1);
%     end
    
    %place the variable inputs into the paths to be changed
    if iscell(paths{idx}{end})
        % variables are bunched (i.e. conjoined)
        conjoined_indx = [conjoined_indx, idx];
        for j = 1:numel(paths{idx})
            paths{idx}{j}{end} = num2str(param(idx));
        end
    else
        % solitary variable
        paths{idx}{end} = num2str(param(idx));
    end
end

%split up conjoined variables
if numel(conjoined_indx)>0
    temp = paths;
    paths = {};
    r=1;
    for j = 1:numel(temp)
        if any(j==conjoined_indx)
            for k = 1:numel(temp{j})
                paths{r} = temp{j}{k};
                r=r+1;
            end
        else
            paths{r} = temp{j};
            r=r+1;
        end
    end
end

%%clean up
clear('A')

f_log=fopen(path_log,'a');  % append to log-file
fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : finished. \n']);
fclose(f_log);

% copyfile(new_path,'C:\Users\BEC Machine\Dropbox\debug\new_path_cp');    % DEBUG