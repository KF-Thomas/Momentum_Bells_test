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

B_trap_bottom=1.0603782e6;
del = -41.03543351e3; %detuning for second beam 3e3
T_pulse_del = +0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
T_Raman_mix=19.20368817e-6;
Gs_mod_R_mix=2.69400894;
phi1_mix=pi;
K_R_mix=0.33063277;%0.31

f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"


B_trap_bottom=1.0751e6;%1.025e6;%1.0603782e6;
del = -45e3;%-41.03543351e3; %detuning for second beam 3e3
T_pulse_del = +0.0e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
T_Raman_mix=23e-6;%19.20368817e-6;
Gs_mod_R_mix=1.8;%2.69400894;
phi1_mix=pi;
phi2=0;
K_R_mix=0.338;%0.34063277;%0.31

f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

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
        
    (2020-08-18)
	dF_Bragg_1=0.037e6;
    dF_Bragg_2=0.037e6;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
    T_Bragg_src=32e-6;%16.7e-6;
    K_Bragg_src=0.2315;%0.2;%0.23;%0.12;%
    Gs_mod_Bragg_src=1.8*T_Bragg_src/16.7e-6;
        
    t0_Bragg_src=nan;

    (2020-08-19)
    ampfun = @(b,x) b(1).*x(:,1).^b(2);
        dF_Bragg_1=0.055e6;
        dF_Bragg_2=0.055e6;
        f1_Bragg_src=f0_AOM-dF_Bragg_1;
        f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
        T_Bragg_src=32e-6;%16.7e-6;
        P_Bragg_src = 4.5; %power in mW 7 to 9 works well
        K_Bragg_src_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;%0.09;%0.2;%0.23;%0.12;%
        K_Bragg_src_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;%0.177;%0.2;%0.23;%0.12;%
        Gs_mod_Bragg_src_1=1.82*T_Bragg_src/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01);
        Gs_mod_Bragg_src_2=1.81*T_Bragg_src/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01);
        
        t0_Bragg_src=nan;

    (2020-08-19)
	ampfun = @(b,x) b(1).*x(:,1).^b(2);
        dF_Bragg_1=0.06e6;
        dF_Bragg_2=0.06e6;
        f1_Bragg_src=f0_AOM-dF_Bragg_1;
        f2_Bragg_src=f0_AOM+dF_Bragg_2;
        
        T_Bragg_src=32e-6;%16.7e-6;
        P_Bragg_src = 4.5; %power in mW 7 to 9 works well
        K_Bragg_src_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;%0.09;%0.2;%0.23;%0.12;%
        K_Bragg_src_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;%0.177;%0.2;%0.23;%0.12;%
        Gs_mod_Bragg_src_1=1.81*T_Bragg_src/16.7e-6*sqrt(2);
        Gs_mod_Bragg_src_2=1.81*T_Bragg_src/16.7e-6*sqrt(2);
        
        t0_Bragg_src=nan;
    (2020-08-19)
    ampfun = @(b,x) b(1).*x(:,1).^b(2);
dF_Bragg_1=0.06e6;
dF_Bragg_2=0.06e6;
f1_Bragg_src_t=f0_AOM-dF_Bragg_1;
f2_Bragg_src_t=f0_AOM+dF_Bragg_2;

T_Bragg_src_t=32e-6;
P_Bragg_src = 3.5; %power in mW 7 to 9 works well
K_Bragg_src_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
K_Bragg_src_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;
Gs_mod_Bragg_src_1=1.82*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_src_2=1.81*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01);

t0_Bragg_src_t=nan;

(2020-11-20)

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

    (2020-08-21)
    dF_Bragg_1=0.086e6;
dF_Bragg_2=0.086e6;
f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

T_Bragg_src_f=32E-6;
P_Bragg_f = 14;
K_Bragg_src_f_1=ampfun([60.117, 0.5638],P_Bragg_f)/2e3;%0.08;
K_Bragg_src_f_2=ampfun([60.117, 0.5638],P_Bragg_f)/2e3;%0.08;
Gs_mod_Bragg_src_f_1=0.5*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_src_f_2=0.5*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
t0_Bragg_src_f=nan;

(2020-08-31)
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

(2020-09-23)
dF_Bragg_1=0.095e6;%0.099e6;%~0.093;%
dF_Bragg_2=0.095e6;%0.099e6;
f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

T_Bragg_src_f=32E-6;
P_Bragg_f = 8.0;%7.2;%20.5;%~7.8
K_Bragg_src_f_1=1.25*ampfun([60.117, 0.5638],P_Bragg_f)/2e3;%0.08;
K_Bragg_src_f_2=1.33*ampfun([132.62, 0.5283],P_Bragg_f)/2e3;%0.08;%[60.117, 0.5638] multiplier ~1.2 1.4
Gs_mod_Bragg_src_f_1=0.9*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);%~0.7 0.95
Gs_mod_Bragg_src_f_2=0.9*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
t0_Bragg_src_f=nan;

(2020-11-18)
dF_Bragg_1=0.095e6;%0.099e6;%~0.093;%
dF_Bragg_2=0.095e6;%0.099e6;
f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

T_Bragg_src_f=32E-6;
P_Bragg_f = 7.9;%7.2;%20.5;%~7.8
K_Bragg_src_f_1=1.25*ampfun([60.117, 0.5638],P_Bragg_f)/2e3;%0.08;
K_Bragg_src_f_2=1.33*ampfun([132.62, 0.5283],P_Bragg_f)/2e3;%0.08;%[60.117, 0.5638] multiplier ~1.2 1.4
Gs_mod_Bragg_src_f_1=0.75*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);%~0.7 0.95
Gs_mod_Bragg_src_f_2=0.75*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
t0_Bragg_src_f=nan;

(2020-12-04)
dF_Bragg_1=0.095e6;%0.099e6;%~0.093;% %%
dF_Bragg_2=0.095e6;%0.099e6; %%
f1_Bragg_src_f=f0_AOM-dF_Bragg_1;
f2_Bragg_src_f=f0_AOM+dF_Bragg_2;

T_Bragg_src_f=31E-6; %%
P_Bragg_f = 7.5 ;%7.2;%20.5;%~7.8 %%
K_Bragg_src_f_1=1.37*ampfun([60.117, 0.5638],P_Bragg_f)/2e3;%0.08; 1.25 1.35
K_Bragg_src_f_2=1.49*ampfun([132.62, 0.5283],P_Bragg_f)/2e3;%0.08;%[60.117, 0.5638] multiplier ~1.2 1.4 1.33 1.45
Gs_mod_Bragg_src_f_1=0.75*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.63806937736142e-01);%~0.7 0.95
Gs_mod_Bragg_src_f_2=0.75*T_Bragg_src_f/4.2E-6*sqrt(2)*sqrt(5.28341254744861e-01);
t0_Bragg_src_f=nan;

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
	(2020-08-20) rough
	ampfun = @(b,x) b(1).*x(:,1).^b(2);
dF_Bragg_1=0.11e6; or 0.1e6
dF_Bragg_2=0.11e6; or 0.1e6
f1_Bragg_mirror=f0_AOM-dF_Bragg_1;
f2_Bragg_mirror=f0_AOM+dF_Bragg_2;

	T_Bragg_mirror=32E-6;
	P_Bragg_src = 5.5; %power in mW 7 to 9 works well
K_Bragg_mirror_1=ampfun([60.117, 0.5638],P_Bragg_src)/2e3;
K_Bragg_mirror_2=ampfun([132.62, 0.5283],P_Bragg_src)/2e3;
% Gs_mod_Bragg_mirror=0.2;%1;%1.03135*T_Bragg_src/4.2E-6;
Gs_mod_Bragg_mirror_1=1.83*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.63806937736142e-01);
Gs_mod_Bragg_mirror_2=1.83*T_Bragg_src_t/16.7e-6*sqrt(2)*sqrt(5.28341254744861e-01);
	

t0_Bragg_mirror=nan;%3.9895e-6;