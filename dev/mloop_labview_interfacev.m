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
B_trap_bottom=1.225e6;%0.9e6;%3.47e6; %2.01e6
del = 3e3; %detuning for second beam
T_pulse_del = 0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"
T_Raman_mix=6.29e-6;
Gs_mod_R_mix=2.0;
phi1_mix=pi;
K_R_mix=0.4;

% % double sided bragg settings
% dF_Bragg_1=0.017e6;%0.045e6
% dF_Bragg_2=0.017e6;%-0.01e5
% f1_Bragg_d=f0_AOM-dF_Bragg_1;
% f2_Bragg_d=f0_AOM+dF_Bragg_2;
% T_Bragg_d=6.7e-6;%15e-6;
% K_Bragg_d=0.26;%28%0.086;
% Gs_mod_Bragg_d=1.21;%1.8
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

params = parameter_reader(i);
parameter_names = {'dF_Raman','del','T_Raman_mix','Gs_mod_R_mix','K_R_mix'};

for ii = 1:numel(params)
    eval([parameter_names{ii} '=' num2str(params(ii))])
end

dF_Raman = dF_Raman*1e6;
del = del*1e3;
T_Raman_mix = T_Raman_mix*1e-6;

%recalculate variables that need to be recalculated
f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"


%% Waveform generation
% null_wf = {{'const',0,sample_rate,1e-6}};
% burst_wf={
%         {'sine',freq,phase_PAL,amp_PAL*sqrt(2),nan,sample_rate,dur_PAL},...
%         {'const',0,sample_rate,1e-6}
%         };
trans_null={'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)};
double_sin_wf={'double_sine',f1_Raman_mix,f2_Raman_mix,phi1_mix,phi2,K_R_mix,K_R_mix,Gs_mod_R_mix,Gs_mod_R_mix,srate_all,T_Raman_mix,T_pulse_del};
% mix_delay_wf ={'const',0, srate_all,T_delay_mix};
% %k=+1 bragg
% sin_1_p_wf={'sine',f1_Bragg_p,phi1,K_Bragg_p,Gs_mod_Bragg_p,srate_all,T_Bragg_p};
% sin_2_p_wf={'sine',f2_Bragg_p,phi2,K_Bragg_p,Gs_mod_Bragg_p,srate_all,T_Bragg_p};
% %k=-1 bragg
% sin_1_m_wf={'sine',f1_Bragg_m,phi1,K_Bragg_m,Gs_mod_Bragg_m,srate_all,T_Bragg_m};
% sin_2_m_wf={'sine',f2_Bragg_m,phi2,K_Bragg_m,Gs_mod_Bragg_m,srate_all,T_Bragg_m};
% %doulbe sided bragg
% sin_1_d_wf={'sine',f1_Bragg_d,phi1,K_Bragg_d,Gs_mod_Bragg_d,srate_all,T_Bragg_d};
% sin_2_d_wf={'sine',f2_Bragg_d,phi2,K_Bragg_d,Gs_mod_Bragg_d,srate_all,T_Bragg_d};


waveform.ch1 = double_sin_wf;
waveform.ch2 = trans_null;
%% Shot sequence settings
path_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';
path_param_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_KeysightMatlab.txt';

%% Interface

addpath('C:\Users\BEC Machine\cloudstor\MATLAB\keysight-33600a')
% Send waveforms
if update_keysight
    send_waveform(waveform.ch1,waveform.ch2,0);
end
%write to log
f_log=fopen(path_log,'a');  % append to log-file
nowdt=datetime('now');
fprintf(f_log,'%d,%.3f,%s,interfacev8,%s\n',...
    i,posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'));
fclose(f_log);
f_log=fopen(path_param_log,'a');  % append to param-log-file
fprintf(f_log,'%d,%.3f,%s,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u\n',...
    i,posixtime(nowdt),shot_info.log_entry,...
    f1_Raman_mix,f2_Raman_mix,T_Raman_mix,K_R_mix,Gs_mod_R_mix,phi1_mix,T_pulse_del,del);
fclose(f_log);