bragg setting log

log for different settings, best settings in chronological order

for reference for some of the settings
%% General configs
f0_AOM=80e6;            % [Hz]     Central AOMs frequency
Ek=84.9e3;              % beam geometry (90 deg)

## magnetic transfer

% B_trap_bottom=1.06310965e6;
% del = -1.95413903e3; %detuning for second beam 3e3
% T_pulse_del = +0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=19.50368817e-6;
% Gs_mod_R_mix=1.06945736;
% phi1_mix=pi;
% K_R_mix=0.14872109;

% B_trap_bottom=1.225e6;%0.9e6;%3.47e6; %2.01e6
% del = 0.05; %detuning for second beam
% T_pulse_del = +0.25e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=6.29e-6;
% Gs_mod_R_mix=2.0;
% phi1_mix=pi;
% K_R_mix=0.4;

% B_trap_bottom=1.135e6;
% del = -1e3; %detuning for second beam 3e3
% T_pulse_del = +0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=6.25e-6;
% Gs_mod_R_mix=2.0;
% phi1_mix=pi;
% K_R_mix=0.4;

% B_trap_bottom=1.06777807e6;
% del = -21.08970536; %detuning for second beam 3e3
% T_pulse_del = +0.0e-6;%delay between pulses
% dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
% T_Raman_mix=13.93735153e-6;
% Gs_mod_R_mix=1.08638709;
% phi1_mix=pi;

## momentum splitting

### transfer to k=0,-1
	(2019-11)
    dF_Bragg_1=0.037e6;%0.045e6
    dF_Bragg_2=0.037e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=16.7e-6;%15e-6;
    K_Bragg_src=0.25;%28%0.086;
    Gs_mod_Bragg_src=1.81;%1.8

    (2020-08-03)
	dF_Bragg_1=0.037e6;
    dF_Bragg_2=0.037e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=16.7e-6;
    K_Bragg_src=0.23;
    Gs_mod_Bragg_src=1.81;
        

### trasnfer to k=+1,0
	(2019-11)
	dF_Bragg_1=0.0e6;%0.045e6
    dF_Bragg_2=0.0e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=16.7e-6;%15e-6;
    K_Bragg_src=0.25;%28%0.086;
    Gs_mod_Bragg_src=1.81;%1.8

### trasnfer to k=+1,0,-1
	(2019-11)
	dF_Bragg_1=0.017e6;%0.045e6
    dF_Bragg_2=0.017e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=6.7e-6;%15e-6;
    K_Bragg_src=0.250;%28%0.086;
    Gs_mod_Bragg_src=1.21;%1.8

### transfer to k=0,-1,-2
	(2020-07-31)
	dF_Bragg_1=0.108e6;%0.037e6;%0.037e6;%0.045e6
    dF_Bragg_2=0.108e6;%0.037e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=4.2E-6;%16.7e-6;%15e-6;
    K_Bragg_src=0.18;%0.25;%28%0.086;
    Gs_mod_Bragg_src=1.0;%1.9;%2.09;%1.81;%1.8

### transfer to k=-1,-2
	(2020-08-05)
	dF_Bragg_1=0.085e6;
    dF_Bragg_2=0.075e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=10.2E-6;
    K_Bragg_src=0.24;
    Gs_mod_Bragg_src=3.0;

    ()
    dF_Bragg_1=0.090e6;
    dF_Bragg_2=0.085e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=10.0E-6;
    K_Bragg_src=0.23;
    Gs_mod_Bragg_src=3.0;

    ()
    dF_Bragg_1=0.090e6;
    dF_Bragg_2=0.085e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=10.0E-6;
    K_Bragg_src=0.195;
    Gs_mod_Bragg_src=2.4;

    ()
    dF_Bragg_1=0.0793e6;
    dF_Bragg_2=0.095e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=11.6E-6;
    K_Bragg_src=0.165;
    Gs_mod_Bragg_src=2.523;

    ()
    dF_Bragg_1=0.089e6;
    dF_Bragg_2=0.092e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=11.5E-6;
    K_Bragg_src=0.168;
    Gs_mod_Bragg_src=2.564;

## transfer to k=-1
%         dF_Bragg_1=0.067e6;
%         dF_Bragg_2=0.067e6;
%         f1_Bragg_src=f0_AOM-dF_Bragg_1;
%         f2_Bragg_src=f0_AOM+dF_Bragg_2;
%         
%         T_Bragg_src=14e-6;
%         K_Bragg_src=0.23;
%         Gs_mod_Bragg_src=1.81;

## mirror

## 50:50 beam spiltter