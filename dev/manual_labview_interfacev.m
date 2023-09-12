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

update_keysight = 1;

num_points = 100;
shots_per_point = 5;

marker = mod(floor((i-1)/shots_per_point),shots_per_point*num_points)+1; %Counts from 1 to num shots before setpt update


sequence = {'k=+1,0,-1',}; %Construct desired experimental sequenc from sequences above

new_path='c:\remote\settings202308Sep172658.xml';%c:\remote\settings202001Sep155855.xml

%% Keysight settings
% General settings
sample_rate=1e9;
max_points=double(4e6);
points_min=double(32);
repeats_max=1e6;
f0_AOM=80e6;            % [Hz]     Central AOMs frequency 
Ek=84.9e3;              % beam geometry (90 deg)
srate_all=sample_rate;
ampfun = @(b,x) b(1).*x(:,1).^b(2);

%%% MAGNETIC TRANSFER pulse
%--------------------------------------------------------------------------
B_trap_bottom=1.0583782e6;
del = -41.03543351e3; %detuning for second beam 3e3
T_pulse_del = +0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
T_Raman_mix=19.20368817e-6;
Gs_mod_R_mix=2.69400894;
phi1_mix=pi;
K_R_mix=0.31063277;

f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

%%% DELAYs between pulses
%--------------------------------------------------------------------------
% T_delay_mix=3000e-6;      % Delay between the SRC and MIX pulse
T_delay_mix=3e-6;      % Delay between the MAG and Bragg pulse
T_delay_mirror=1500e-6;%500e-6; % Delay between the Bragg pulse and Mirror pulse

%%% phases
%--------------------------------------------------------------------------
% NOTE: separated out since we find zero gives fine results
phi1=0;
phi2=0;

%%% MIRROR pulse
%--------------------------------------------------------------------------
dF_Bragg_1=0.12e6;%0.0849e6;%
dF_Bragg_2=0.12e6;%0.0849e6;%
f1_Bragg_mirror=f0_AOM-dF_Bragg_1;
f2_Bragg_mirror=f0_AOM+dF_Bragg_2;

T_Bragg_mirror=35E-6;%32E-6;
P_Bragg_src = 7.0; %power in mW 7 to 9 works well ~5.9
K_Bragg_mirror_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
K_Bragg_mirror_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;

Gs_mod_Bragg_mirror_1=1.95*T_Bragg_mirror/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01);%~1.83
Gs_mod_Bragg_mirror_2=1.95*T_Bragg_mirror/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01);%~1.83

t0_Bragg_mirror=nan;%3.9895e-6;

sinc_scale_Bragg_mirror=6e-6;%5.5e-6;%5.3e-6;%5.3e-6;
Amp_sinc_Bragg_mirror=sqrt(5);%sqrt(10.0);%0.2;%

wf_mirror_pulse_d = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3).*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_mirror_pulse_1 = @(b,t) ampfun([60.117, 0.5638],abs(wf_mirror_pulse_d(b,t)).^2)/2e3.*sin(2*pi*b(2)*t).*sign(wf_mirror_pulse_d(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%
wf_mirror_pulse_2 = @(b,t) ampfun([132.62, 0.5283],abs(wf_mirror_pulse_d(b,t)).^2)/2e3.*sin(2*pi*b(2)*t).*sign(wf_mirror_pulse_d(b,t));%sinc((t-b(1)/2)./b(4)).*b(3).*sin(2*pi*b(2)*t);%

%%% 50:50 Beam Splitter pulse
%--------------------------------------------------------------------------
dF_Bragg_1=0.1e6;%0.11e6;% or 0.1e6
dF_Bragg_2=0.1e6;%0.11e6;% or 0.1e6
f1_Bragg_splitter=f0_AOM-dF_Bragg_1;
f2_Bragg_splitter=f0_AOM+dF_Bragg_2;

T_Bragg_splitter=32E-6;
P_Bragg_splt = 5.5; %power in mW 7 to 9 works well
K_Bragg_splitter_1=ampfun([60.117, 0.5638],P_Bragg_splt)/2e3;
K_Bragg_splitter_2=ampfun([132.62, 0.5283],P_Bragg_splt)/2e3;
% Gs_mod_Bragg_mirror=0.2;%1;%1.03135*T_Bragg_src/4.2E-6;
Gs_mod_Bragg_splitter_1=1.83*T_Bragg_splitter/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_splitter_2=1.83*T_Bragg_splitter/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01);

t0_Bragg_splitter=nan;%3.9895e-6;

%%% Momentum splitting
%--------------------------------------------------------------------------

%%% Bragg splitting: |k=0> |--> |k=0> + |k=-1K> + |k=-2K>
dF_Bragg_1=0.095e6;%0.099e6;%~0.093;%
dF_Bragg_2=0.095e6;%0.099e6;
f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

T_Bragg_src_f=32E-6;
P_Bragg_f = 7.2;%7.2;%20.5;%~7.8
K_Bragg_src_f_1=1.25*ampfun([60.117, 0.5638],P_Bragg_f)/2e3;%0.08;
K_Bragg_src_f_2=1.33*ampfun([132.62, 0.5283],P_Bragg_f)/2e3;%0.08;%[60.117, 0.5638] multiplier ~1.2 1.4
Gs_mod_Bragg_src_f_1=0.9*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);%~0.7 0.95
Gs_mod_Bragg_src_f_2=0.9*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
t0_Bragg_src_f=nan;

%%% Bragg splitting: |k=0> |--> |k=+1K> + |k=0> + |k=-1K>
dF_Bragg_1= 42.48e3/2;
dF_Bragg_2= 42.48e3/2;
f1_Bragg_sym_f=f0_AOM-dF_Bragg_1;
f2_Bragg_sym_f=f0_AOM+dF_Bragg_2;

T_Bragg_sym_f=500e-6;

K_Bragg_sym_f_1= 0.6;
K_Bragg_sym_f_2= 0.6; 

Gs_mod_Bragg_sym_f_1=40e-6;
Gs_mod_Bragg_sym_f_2=40e-6;
t0_Bragg_sym_f=nan;

phi1_Bragg = 0;
phi2_Bragg = 0;

wf_mirror_pulse = @(b,t) exp(-((t-b(1)/2)./b(4)).^2).*b(3);%

wf_bragg_sym_pulse = @(b,t) sinc((t-b(1)/2)./b(4)).*b(3);%.*cos(pi*(t-b(1)/2)/(b(1))).^2;
wf_bragg_sym_pulse_1 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));
wf_bragg_sym_pulse_2 = @(b,t) wf_mirror_pulse(b,t).*sin(2*pi.*(b(2)-b(6).*t).*t+b(5));


%%% Bragg splitting: |k=0> |--> |k=0> + |k=-1K>
dF_Bragg_1=0.06e6;
dF_Bragg_2=0.06e6;
f1_Bragg_src_t=f0_AOM-dF_Bragg_1;
f2_Bragg_src_t=f0_AOM+dF_Bragg_2;

T_Bragg_src_t=32e-6;
P_Bragg_src = 3.6; %power in mW 7 to 9 works well
K_Bragg_src_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
K_Bragg_src_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;
Gs_mod_Bragg_src_1=1.82*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_src_2=1.81*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01);

t0_Bragg_src_t=nan;

%%% Bragg splitting: |k=0> |--> |k=-1K> + |k=-2K>
dF_Bragg_1=0.085e6;
dF_Bragg_2=0.075e6;
f1_Bragg_src_b=f0_AOM-dF_Bragg_1;
f2_Bragg_src_b=f0_AOM+dF_Bragg_2;

T_Bragg_src_b=10.2e-6;
K_Bragg_src_b=0.24;
Gs_mod_Bragg_src_b=3.0;

%%% PAL settings
%--------------------------------------------------------------------------
freq = 1.2e6;
amp_PAL = sqrt(2)*0.650; %Vrms
phase_PAL = 0;
cycles = 6;
dur_PAL = cycles/freq;

%% Iteration through parameters
Amp_sinc_Bragg_mirror_vec = flip(sqrt([0 0.25 10 19]));
Amp_sinc_Bragg_mirror = Amp_sinc_Bragg_mirror_vec(1);%marker

Gs_mod_Bragg_sym_vec = [0:5:500];

Gs_mod_Bragg_sym_f_1= Gs_mod_Bragg_sym_vec(marker) * 1e-6;
Gs_mod_Bragg_sym_f_2= Gs_mod_Bragg_sym_vec(marker) * 1e-6;


%% Waveform generation
ch1_raw={}; %waveform for chanel 1
ch2_raw={}; %waveform for chanel 2
for ii = 1:length(sequence) %run through each segment
    segment = sequence{ii}; %the current segment
    switch segment %add to the waveforms the desired segment
        case 'mag_transfer'
            %%% Magentic transfer
            ch1_raw=[ch1_raw(:)',...
                {{'double_sine',f1_Raman_mix,f2_Raman_mix,phi1_mix,phi2,K_R_mix,...
                K_R_mix,Gs_mod_R_mix,Gs_mod_R_mix,srate_all,T_Raman_mix,T_pulse_del}}
                ];
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)}}
                ];
        case 'k=0,-1,-2'
            %%% Full Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'sine',    f1_Bragg_src_f       ,phi1,          K_Bragg_src_f_1,       Gs_mod_Bragg_src_f_1, srate_all,   T_Bragg_src_f,   t0_Bragg_src_f}},...
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'sine',    f2_Bragg_src_f       ,phi2,          K_Bragg_src_f_2,       Gs_mod_Bragg_src_f_2, srate_all,   T_Bragg_src_f,    t0_Bragg_src_f}},...
                ];
        case 'k=0,-1'
            %%% Top Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'sine',    f1_Bragg_src_t       ,phi1,          K_Bragg_src_1,       Gs_mod_Bragg_src_1, srate_all,   T_Bragg_src_t, t0_Bragg_src_t}},...
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'sine',    f2_Bragg_src_t       ,phi2,          K_Bragg_src_2,       Gs_mod_Bragg_src_2, srate_all,   T_Bragg_src_t, t0_Bragg_src_t}},...
                ];
        case 'k=-1,-2'
            %%% Bottom Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'sine',    f1_Bragg_src_b       ,phi1,          K_Bragg_src_b,       Gs_mod_Bragg_src_b, srate_all,   T_Bragg_src_b}},...
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_Raman_mix+abs(T_pulse_del)}},...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'sine',    f2_Bragg_src_b       ,phi2,          K_Bragg_src_b,       Gs_mod_Bragg_src_b, srate_all,   T_Bragg_src_b}},...
                ];

        case 'k=+1,0,-1'
            %%% Full Halo
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_sym_pulse_1,  srate_all,  T_Bragg_sym_f,  f1_Bragg_sym_f,  K_Bragg_sym_f_1,  Gs_mod_Bragg_sym_f_1, phi1_Bragg,0}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mix}},...
                {{'arb',   wf_bragg_sym_pulse_2,  srate_all,  T_Bragg_sym_f,  f2_Bragg_sym_f,  K_Bragg_sym_f_2,  Gs_mod_Bragg_sym_f_2, phi2_Bragg,0}}
                ];

        case 'mirror'
            %%% Mirror pulse
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_mirror}},...
                {{'arb',   wf_mirror_pulse_1,  srate_all,  T_Bragg_mirror,  f1_Bragg_mirror,  Amp_sinc_Bragg_mirror,  sinc_scale_Bragg_mirror}}
                ];
%             {'sine',    f1_Bragg_mirror       ,phi1,          K_Bragg_mirror_1,       Gs_mod_Bragg_mirror_1, srate_all,   T_Bragg_mirror, t0_Bragg_mirror}
%             
%             
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_mirror}},...
                {{'arb',   wf_mirror_pulse_2,  srate_all,   T_Bragg_mirror, f2_Bragg_mirror,  Amp_sinc_Bragg_mirror,  sinc_scale_Bragg_mirror}}
                ];
%             {'sine',    f2_Bragg_mirror       ,phi2,          K_Bragg_mirror_2,       Gs_mod_Bragg_mirror_2, srate_all,   T_Bragg_mirror, t0_Bragg_mirror}
%             
%             
        case 'splitter'
            %%% 50:50 Beam splitter pulse
            ch1_raw=[ch1_raw(:)',...
                {{'const',0, srate_all,T_delay_splitter}},...
                {{'sine',    f1_Bragg_splitter       ,phi1,          K_Bragg_splitter_1,       Gs_mod_Bragg_splitter_1, srate_all,   T_Bragg_splitter, t0_Bragg_splitter}}
                ];
            
            ch2_raw=[ch2_raw(:)',...
                {{'const',0, srate_all,T_delay_splitter}},...
                {{'sine',    f2_Bragg_splitter       ,phi2,          K_Bragg_splitter_2,       Gs_mod_Bragg_splitter_2, srate_all,   T_Bragg_splitter, t0_Bragg_splitter}}
                ];
        otherwise
            error('invalid sequence segment');
    end
end

%% Shot sequence settings
path_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';
path_param_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_KeysightMatlab.txt';

%% convert waveforms to printables for logging
ch1_waveform_str = '';
ch2_waveform_str = '';
addpath('C:\Users\BEC Machine\Documents\MATLAB\Momentum_Bells_test\dev')
addpath('C:\Users\BEC Machine\cloudstor\PROJECTS\keysight-33600a\ch_to_waveforms.m')
for waveforms = 1:numel(ch1_raw)
    if waveforms>1
        ch1_waveform_str = [ch1_waveform_str,', '];
    end
    ch1_waveform_str = [ch1_waveform_str,cell2str(array2str(ch1_raw{waveforms}))];
end
for waveforms = 1:numel(ch2_raw)
    if waveforms>1
        ch2_waveform_str = [ch2_waveform_str,', '];
    end
    ch2_waveform_str = [ch2_waveform_str,cell2str(array2str(ch2_raw{waveforms}))];
end

ch1_waveform_str = replace(ch1_waveform_str,"'",'');
ch2_waveform_str = replace(ch2_waveform_str,"'",'');

%% Interface

addpath('C:\Users\BEC Machine\OneDrive - Australian National University\PROJECTS\keysight-33600a')%C:\Users\BEC Machine\OneDrive - Australian National University\PROJECTS\keysight-33600a\WaveformGenMain.m
%C:\Users\BEC Machine\cloudstor\MATLAB\keysight-33600a
% shot_info = shots.(shot_sequence{1});
% new_path=shot_info.LVfile;
% Send waveforms
chanels_dev1={ch_to_waveforms(ch1_raw),ch_to_waveforms(ch2_raw)};
if update_keysight && mod((i-1),shots_per_point) == 0
    send_segments(chanels_dev1,1);
end
%write to log
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