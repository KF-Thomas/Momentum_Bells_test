%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interface between labview and matlab

%INPUT FROM LABVIEW
%   i          - integer,itteration number
%   auto_enable- boolean
%   file       - string : tells program where to look for setting files
%   file_exact - string
%   mloop      - boolean(legacy) do you want to interface with mloop or use matlab to scan over variables
%   control    - boolean(legacy)

%return
%exact_line- double (legacy just return a zero)
%new_path 

%% USER SETTINGS
exp_type = '';
update_keysight = 1;
% new_path='c:\remote\settings202023Jul083400.xml'; %mag transfer
% new_path='c:\remote\settings202021Jul085215.xml';%momentum transfer
new_path='c:\remote\settings202030Jul080637.xml';
%% Keysight settings
% General settings
sample_rate=1e9;
max_points=double(4e6);
points_min=double(32);
repeats_max=1e6;
f0_AOM=80e6;            % [Hz]     Central AOMs frequency 
Ek=84.9e3;              % beam geometry (90 deg)
srate_all=sample_rate;

T_delay_mix=10e-6;      % Delay between the SRC and MIX pulse
phi1=0;
phi2=0;

% Magnetic Transfer settings
B_trap_bottom=1.0583782e6;
del = 0;%-41.03543351e3; %detuning for second beam 3e3
T_pulse_del = +0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
T_Raman_mix=19.20368817e-6;
Gs_mod_R_mix=2.69400894;
phi1_mix=pi;
K_R_mix=0.31063277;
f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

% B_trap_bottom=1.135e6;
% del = -1e3; %detuning for second beam 3e3
% T_pulse_del = +0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=6.25e-6;
% Gs_mod_R_mix=2.0;
% phi1_mix=pi;
% K_R_mix=0.4;
% f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
% f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

% % double sided bragg settings
dF_Bragg_1=0.020e6;%0.045e6
dF_Bragg_2=0.020e6;%-0.01e5
f1_Bragg_d=f0_AOM-dF_Bragg_1;
f2_Bragg_d=f0_AOM+dF_Bragg_2;
% T_Bragg_d=6e-6;%15e-6;
T_Bragg_d=35e-6;%15e-6;
K_Bragg_d_1=0.155;%28%0.086;
Gs_mod_Bragg_d_1=1.3;%1.8
K_Bragg_d_2=0.155;%28%0.086;
Gs_mod_Bragg_d_2=1.3;%1.8
% 
% % k=+1 diffraction
% dF_Bragg_1=0.0e6;%0.045e6
% dF_Bragg_2=0.0e6;%-0.01e5
% f1_Bragg_p=f0_AOM-dF_Bragg_1;
% f2_Bragg_p=f0_AOM+dF_Bragg_2;
% T_Bragg_p=16.7e-6;%15e-6;
% K_Bragg_p=0.26;%28%0.086;
% Gs_mod_Bragg_p=1.81;%1.8
% 
% % k=-1 diffraction
% dF_Bragg_1=0.037e6;%0.045e6
% dF_Bragg_2=0.037e6;%-0.01e5
% f1_Bragg_m=f0_AOM-dF_Bragg_1;
% f2_Bragg_m=f0_AOM+dF_Bragg_2;
% T_Bragg_m=16.7e-6;%15e-6;
% K_Bragg_m=0.25;%28%0.086;
% Gs_mod_Bragg_m=1.81;%1.8
% 
% % PAL settings
% freq = 1.2e6;
% amp_PAL = sqrt(2)*0.650; %Vrms
% phase_PAL = 0;
% cycles = 6;
% dur_PAL = cycles/freq;

%% Iteration through parameters
addpath('C:\Users\BEC Machine\Documents\MATLAB\Momentum_Bells_test\dev')
params = parameter_reader(i);
%parameters being controlled
% parameter_names = {'dF_Raman','del','T_Raman_mix','Gs_mod_R_mix','K_R_mix'};
% parameter_names = {'dF_Bragg_1','dF_Bragg_2','K_Bragg_d','Gs_mod_Bragg_d'};
% parameter_names = {'dF_Bragg_1','dF_Bragg_2','K_Bragg_d','Gs_mod_Bragg_d','T_Bragg_d'};
% parameter_names = {'dF_Raman','del','Gs_mod_R_mix','K_R_mix'};
% parameter_names = {'dF_Raman','T_pulse_del','Gs_mod_R_mix','K_R_mix'};
% parameter_names = {'n_0','n_1','n_2','n_3'};
parameter_names = {'dF_Bragg','K_Bragg_d_1','K_Bragg_d_2','Gs_mod_Bragg_d_1','Gs_mod_Bragg_d_2'};

for ii = 1:numel(params)
    eval([parameter_names{ii} '=' num2str(params(ii))])
end
%update the nuller
% update_nuller(n_0,n_1,n_2,n_3,i)

% dF_Raman = dF_Raman*1e6;
% del = del*1e3;
% T_Raman_mix = T_Raman_mix*1e-6;
% T_pulse_del = T_pulse_del*1e-6;
% 
%recalculate variables that need to be recalculated
% f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
% f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

% dF_Bragg_1 = dF_Bragg_1*1e3;
% dF_Bragg_2 = dF_Bragg_2*1e3;
% f1_Bragg_d=f0_AOM-dF_Bragg_1;
% f2_Bragg_d=f0_AOM+dF_Bragg_2;
% 
% T_Bragg_d = T_Bragg_d*1e-6;


dF_Bragg_1 = dF_Bragg*1e3;
dF_Bragg_2 = dF_Bragg*1e3;
f1_Bragg_d=f0_AOM-dF_Bragg_1;
f2_Bragg_d=f0_AOM+dF_Bragg_2;

% T_Bragg_d = T_Bragg_d*1e-6;

%% Waveform generation
% null_wf = {{'const',0,sample_rate,1e-6}};
% burst_wf={
%         {'sine',freq,phase_PAL,amp_PAL*sqrt(2),nan,sample_rate,dur_PAL},...
%         {'const',0,sample_rate,1e-6}
%         };
trans_null={'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)};
double_sin_wf={'double_sine',f1_Raman_mix,f2_Raman_mix,phi1_mix,phi2,K_R_mix,K_R_mix,Gs_mod_R_mix,Gs_mod_R_mix,srate_all,T_Raman_mix,T_pulse_del};
mix_delay_wf ={'const',0, srate_all,T_delay_mix};
% %k=+1 bragg
% sin_1_p_wf={'sine',f1_Bragg_p,phi1,K_Bragg_p,Gs_mod_Bragg_p,srate_all,T_Bragg_p};
% sin_2_p_wf={'sine',f2_Bragg_p,phi2,K_Bragg_p,Gs_mod_Bragg_p,srate_all,T_Bragg_p};
% %k=-1 bragg
% sin_1_m_wf={'sine',f1_Bragg_m,phi1,K_Bragg_m,Gs_mod_Bragg_m,srate_all,T_Bragg_m};
% sin_2_m_wf={'sine',f2_Bragg_m,phi2,K_Bragg_m,Gs_mod_Bragg_m,srate_all,T_Bragg_m};
% %doulbe sided bragg
sin_1_d_wf={'sine',f1_Bragg_d,phi1,K_Bragg_d_1,Gs_mod_Bragg_d_1,srate_all,T_Bragg_d};
sin_2_d_wf={'sine',f2_Bragg_d,phi2,K_Bragg_d_2,Gs_mod_Bragg_d_2,srate_all,T_Bragg_d};


waveform.ch1 = {double_sin_wf,mix_delay_wf,sin_1_d_wf};
waveform.ch2 = {trans_null,mix_delay_wf,sin_2_d_wf};
% waveform.ch1 = {double_sin_wf};
% waveform.ch2 = {trans_null};

%% Shot sequence settings
path_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';
path_param_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_KeysightMatlab.txt';

%% convert waveforms to printables for logging
ch1_waveform_str = '';
ch2_waveform_str = '';
for waveforms = 1:numel(waveform.ch1)
    if waveforms>1
        ch1_waveform_str = [ch1_waveform_str,', '];
    end
    ch1_waveform_str = [ch1_waveform_str,cell2str(array2str(waveform.ch1{waveforms}))];
end
for waveforms = 1:numel(waveform.ch2)
    if waveforms>1
        ch2_waveform_str = [ch2_waveform_str,', '];
    end
    ch2_waveform_str = [ch2_waveform_str,cell2str(array2str(waveform.ch2{waveforms}))];
end

ch1_waveform_str = replace(ch1_waveform_str,"'",'');
ch2_waveform_str = replace(ch2_waveform_str,"'",'');

%% Interface

addpath('C:\Users\BEC Machine\cloudstor\MATLAB\keysight-33600a')
% Send waveforms
if update_keysight
    send_waveform(waveform.ch1,waveform.ch2,0);
end
%write to log
f_log=fopen(path_log,'a');  % append to log-file
nowdt=datetime('now');
fprintf(f_log,'shot num:%d, posixtime:%.3f, date:%s, matlab:interfacev8, labview settings:%s\n',...
    i,posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),new_path);
fclose(f_log);
f_log=fopen(path_param_log,'a');  % append to param-log-file
fprintf(f_log,'shot num:%d, posixtime:%.3f, ch1 waveform: %s, ch2 waveform: %s\n',...
    i,posixtime(nowdt),...
    ch1_waveform_str,ch2_waveform_str);
fclose(f_log);
pause(0.1)
% new_path='c:\remote\settings202023Jul120059.xml';