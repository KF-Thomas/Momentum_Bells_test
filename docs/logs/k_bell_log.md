s-wave multihalo
c:\remote\settings201901Oct145400.xml
halo 3s shunt, clean seq
c:\remote\settings201901Oct151809.xml

c:\remote\settings201902Oct095311.xml
c:\remote\settings201902Oct114828.xml
shunt=0.6 V, quad=0.5 V
c:\remote\settings201902Oct135534.xml
shunt=0.75 V, quad=0.5 V
c:\remote\settings201902Oct145138.xml

my attempts
c:\remote\settings201903Oct124134.xml
c:\remote\settings201903Oct133928.xml
c:\remote\settings201903Oct143310.xml

c:\remote\settings201904Oct143628.xml

c:\remote\settings201904Oct161536.xml

based of mutilhalo
c:\remote\settings201903Oct164108.xml
c:\remote\settings201903Oct165706.xml
c:\remote\settings201904Oct100509.xml

current
c:\remote\settings201904Oct165528.xml

c:\remote\settings201909Oct135953.xml

david's settings
c:\remote\settings201802Jun153228.xml

my updates to davids: c:\remote\settings201909Oct150832.xml

incorporatyingh bragg beams: c:\remote\settings201910Oct144643.xml

end of day 2019-10-10: c:\remote\settings201910Oct170057.xml

my seetings with davids ramp: c:\remote\settings201911Oct113421.xml
c:\remote\settings201911Oct115118.xml
c:\remote\settings201911Oct132212.xml

end of day 2019-10-11: c:\remote\settings201911Oct164346.xml
end of day 2019-10-14: c:\remote\settings201914Oct151232.xml
c:\remote\settings201915Oct092527.xml

David's ramp settings
Shunt:
sub		wvf		int		fin		dwell		dur		tau
0		con		0		0		100			1900
1		con		0		0		0.1			154.2
2		con		0		0		0.1			0.4
3		con		10		0		10			700
4		lin		10		5.5		10			200
5		con		5.5		-0.2	10			300
6		lin		5.5		-0.05	10			250
7		con 	-0.05	0		500			16000
8		lin		-0.05	0.2971	5			195
9	 	lin		0.2971	0.4773	5			195
10		lin		0.4773	0.5866	5			195
11		lin		0.5866	0.6529	5			195
12		lin		0.6529	0.6931	5			195
13		lin		0.6931	0.7175	5			195
14		lin		0.7175	0.7323	5			195
15		lin		0.7323	0.7413	5			195
16 		lin		0.7413	0.7467	5			195
17		lin		0.7467	0.75	5			200
18		con		0.75	0		5			1290
19		con		0.75	0		0.1			14
20		exp		0.75	0.78	0.1			0.8		0.2

Quad:
sub		wvf		int		fin		dwell		dur		tau
0		con		1.4		0		100			1900
1		con		1.4		0		0.1			154.2
2		con		1.4		0		0.1			0.4
3		con		5.25	5.25	10			700
4		lin		5.25	7.2		10			200
5		con		7.2		0		10			300
6		lin		7.2		5.74	10			250
7		con		5.74	7.71	100			6900
8		lin		5.74	3.4		20			1000
9		con		3.4		3.6		50			7600
10		con		3.4		-0.2	0.1			39
11		con		3.4		-0.2	100			700
12		lin		3.4		1.9752	5			195
13		lin		1.9752	1.1933	5			195
14		lin		1.1933	0.7642	5			195
15		lin		0.7642	0.5286	5			195
16		lin		0.5286	0.3994	5			195
17		lin		0.3994	0.3285	5			195
18		lin		0.3285	0.2895	5			195
19		lin		0.2895	0.2682	5			195
20		lin		0.2682	0.2564	5			195
21		lin		0.2564	0.25	5			200
22		con		0.25	0		100			1200
23		exp		0.25	0.26	0.1			0.8		0.2


## 2019-10-15
shots 795 and 796 regained bragg driffraction I had miss judge trap switch off
currently achieving bragg difraction c:\remote\settings201915Oct112216.xml
keysight trig at 23.455, trap switch off at ~23.459
attempting to see some halos (803-815)
zero detuning (816-982)

# attempting to achieve equal state splitting among k=-1,0,+1
set detuning to be zero, i.e. dF_Bragg_src=0;
short time (3 mus) high power (0.5) produces side orders and has the most pop in k=0
side orders begin to die off around 10 mus

# going to move position of bragg to after trap switch off
doing so causes many more diffraction order to appear (see shot 987-988 for 3ms diffraction after trap switch off, 989-990 for 2ms, 991-992 for 1ms, and 993-994 for at trap switch off [well just slightly before])
moving back to 1 ms before trap switch off (see 997 - 998 for reference)
13mus pulse no longer has side modes though transfer greatly decreased (999-1000)
11 mus pulse have slight side modes but transfer is ok (1001 - 1002)
halfing the power of the pulse (to 0.25) seems to equalize the state transfer (and remove side modes) though k=+/-1 is much smaller than k=0 (1003-4)
reducing pulse length with smaller power (0.25 and 9 mus) increases state transfer while keeping side modes suppresed (1005-1006)
reducing evap (-20kHz) to remove saturation and see how populations compare
still saturating reducing another 20 kHz (again 20 kHz)
10 kHz up (from what I can see fairly equal transfer 0.25 and 7 mus)
back up 30 kHz (back up 20 kHz)
optimal parameters (T=7.8e-6,K=0.26,Gs=2)

will attempt in trap state rotation
with raman activating in trap but no kick yet 1037
with kick activating 1ms after trap switch off (1038 to 1039 nuller was on during kick, 1040 1060 nuller off and kick activating properly, can't see a difference)
kick without in trap transfer 1063 (still no noticble effect)
1069 kick activating (slight shift in z)
increasing raman power to see if there is transfer (1075)
achieved state transfer by removing bragg from keysight (1079)

c:\remote\settings201915Oct164955.xml
c:\remote\settings201915Oct165559.xml

current settings with davids nuller config: c:\remote\settings201915Oct165925.xml

## 2019-10-16
Bragg laser (1550 nm) remained on the whole night with power remaining stable
moving keysight trig 3 ms after trap switch off

## 2019-10-17
current settings: c:\remote\settings201917Oct122440.xml
switching to old settings to understand why i can't move bec with push coils
getting some movement in the bec, will try to get the two peaks back (reducing bragg 1226)
seeing two peak (1226) but not well separated, increasing kick (1227), c:\remote\settings201917Oct124413.xml
removing first push coil (1228) greatly increased separation of the clouds, reference for no push (1229)
(1230) testin shorter kick, (1231) even shrt kick
transporting over to setting with nuller in them (current:c:\remote\settings201917Oct125550.xml, swap:c:\remote\settings201917Oct125907.xml)
another no kick ref (1236)
moving keysight to after trap switch off (1238), almost no transger

current transfer (ch1 double sine)

B_trap_bottom=3.475e6; %2.01e6
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

T_Raman_mix=5.5e-6;
Gs_mod_R_mix=3;
phi1_mix=pi;
phi2=0;
K_R_mix=0.7;
see (1255 for ref) trigeering at 23.458
attemting two beam state transfer
getting ok transfer with two beams (large leakage to mj=-1)

c:\remote\settings201917Oct171006.xml

## 2019-10-18
Attmepting to get transfer with RF antenna (note detuning/ bias in trap is 3.5 MHz)
sees to be about 2.05 +/- 0.05 MHz

keep amplitude constant and varry time and trap bottom

following settings giving good single beam transfer
%%% MIXING pulse
B_trap_bottom=1.08e6;%0.9e6;%3.47e6; %2.01e6
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

T_Raman_mix=5.5e-6;
Gs_mod_R_mix=2;
phi1_mix=pi;

% % MIX 1
K_R_mix=0.4;

with nuller to set to point along the lvis
best so far 
%%% MIXING pulse
B_trap_bottom=1.12e6;%0.9e6;%3.47e6; %2.01e6
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

T_Raman_mix=5.9e-6;
Gs_mod_R_mix=2;
phi1_mix=pi;


% % MIX 1
K_R_mix=0.4;

%%%%
%%% MIXING pulse
B_trap_bottom=1.22e6;%0.9e6;%3.47e6; %2.01e6
del = 0;
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

T_Raman_mix=6.4e-6;
Gs_mod_R_mix=2;
phi1_mix=pi;

% % MIX 1
K_R_mix=0.4;

Ajusted polarisation and alignment (see 1853 to 1880)

Introduced a time delay variable (unsure if implementation is optimal)

adjusted nuller and alignment (again) 380 in the mj=1 state is the number to beat (~79.8% transfer)
going param by param optimise to see what we can get
T_Raman_mix
6.35e-6
2007-2009 shots to beat 
6.3e-6
2013-2015
2022-2024
6.29e-6
2028-2030

testing a tukeywin with 6.291e-6
2037-2039 (no noticeable improvement)
trying tukey again with 6.29e-6
2040-2042 (again seems worse)
trying with chebwin 6.29e-6
2043-2045 (much worse, by a factor of 2 ish)

I think the current parameters are just not optimised for these windows

moving on to del (intial value is 0.05)
0.15
2046-2047
0.3
2048-2051
-0.3
2052-2054
0.0
2055-2057
0.06
2058-
0.04
2061-2064
keeping at 0.05

moving on to phi1_mix (intial value pi)
pi*1.1: 2064-2066
pi*0.9: 2067-2069
pi*0.8: 2070-2072
remaining at (pi)

T_pulse_del (intial -0.2e-6)
-0.8e-6: 2073-2075
0.0e-6: 2076-2078
-0.3e-6: 2079-2081
-0.2e-6: 2082-2084
-0.1e-6: 2085 - 2087
+0.2e-6: 2088-2090
+0.25e-6: 2091-2096 (best so far)

going on to Gs_mod_R_mix (int 2.0)
2.2: 2097 - 2099
1.8: 2100 - 2102
2.5: 2103 - 2105
remaining (2.0)

K_R_mix (0.49)
0.48: 2106-2108
0.5: 2109-2114
0.4: 2115 - 2117, 2124-2126
0.3: 2118 - 2120
0.39: 2121 - 2123
leaving at (0.4)

lastly Im going to adjust B_trap_bottom ( int 1.225e6)
1.235e6: 2127 - 2129
1.265e6: 2130 - 2132
keeping at (1.225e6)

last thing on the transfer list im going to try the nuller again

we have our transfer sequence! 

B_trap_bottom=1.225e6;%0.9e6;%3.47e6; %2.01e6
del = 0.05; %detuning for second beam
T_pulse_del = +0.25e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

T_Raman_mix=6.29e-6;
Gs_mod_R_mix=2.0;
phi1_mix=pi;
% % MIX 1
K_R_mix=0.4;
transfer to mj =-1,0,+1 of about [17%,75%,8%] (note transfer mj = 0 my be higher, it's currently saturated pretty strongly)
knocked synth down (60 kHz)

messing around with keysight trigger time (23.462)

Reintroducing the bragg pulse into the mix (might need to blast magnetically sesative off the detector to see Bragg clearly)
Ajusting the delay between mixing and bragg
intial bragg settings
% % dF_Bragg_src=abs(Ek);
% % % Bragg: |1_0> |--> |1_0> + |1_2K> + |1_-2K>
dF_Bragg_src=0;
f1_Bragg_src=f0_AOM-dF_Bragg_src/2;
 f2_Bragg_src=f0_AOM+dF_Bragg_src/2;
    
    % ~50/50: 24-10-2017
    T_Bragg_src=7.8e-6;%15e-6;
    K_Bragg_src=0.26;%0.086;
    Gs_mod_Bragg_src=2.0;
    T_delay_mix=2000e-6;      % Delay between the SRC and MIX pulse

 Noticed the nuller seems to be pushing the condensate through the evap
 ramp point was 3.56 moved it 3.58
 not enough moving final point to 3.62
 moving the synth doesn't seem to help, will try just turning it off
 that did get rid of it

 going to just leave it in for now , settings : c:\remote\settings201921Oct161725.xml
 nuller was switching at 23.468 now switching off at 23.464

 c:\remote\settings201922Oct093724.xml

## 2019-10-22
Re did the cables and regained Jacobs big push, clears off all the magneticlly sensative atoms, currently connected to 6602 Counter 4

Going to attempt to get bragg again, looks ok (see 2264) have about four or five orders
going to add a detuning for the fact the atoms are moving

Note that the undifracted mj=0's land at 3.8769 seconds (with dld triggered at 20s)
need to account for the doppler shift
getting somewhere (the freq shift due to doppler should be on the order of 10's of kHz
shot 2380 single side diffraction settings 
dF_Bragg_1=0.045e6;
    dF_Bragg_2=0.01e5;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    % ~50/50: 24-10-2017
    T_Bragg_src=8.7e-6;%15e-6;
    K_Bragg_src=0.28;%0.086;
    Gs_mod_Bragg_src=1.4;
    T_delay_mix=10e-6;
slightly better single side (still have some side lobes though)    
dF_Bragg_1=0.045e6;
    dF_Bragg_2=-0.01e5;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=8.7e-6;%15e-6;
    K_Bragg_src=0.28;%0.086;
    Gs_mod_Bragg_src=1.4;

flipping detunings    (2383 vs 2384)


 
 equal(ish) diffraction
 dF_Bragg_1=0.045e6;
    dF_Bragg_2=-0.01e5;
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=9.7e-6;%15e-6;
    K_Bragg_src=0.28;%0.086;
    Gs_mod_Bragg_src=1.0;  

 There's a complicated (not very clear) relation ship between all the variables and the diffraction amount 

 good equal diffraction
     dF_Bragg_1=0.015e6;%0.045e6
    dF_Bragg_2=0.015e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=16.7e-6;%15e-6;
    K_Bragg_src=0.28;%0.086;
    Gs_mod_Bragg_src=1.6;

    halo attemp 2424-2826

    nuller set to 0 in Z (2842-2843)
    nuller 3 0.1 (2844-2845)
    nuller 3 -0.5 (2048-2049)
    nuller 3 -1.5 (2850-2851)
    nuller 3 +1.5 (2852-2853)
    nuller 3 +2.5 (2854-2855)
    nuller 3 -2.5 (2856-2859)

    %%% MIXING pulse
B_trap_bottom=1.225e6;%0.9e6;%3.47e6; %2.01e6
del = 3e3; %detuning for second beam
T_pulse_del = +0.75e-6;%delay between pulses
dF_Raman=-(B_trap_bottom);     %[Hz]    Raman detuning
f1_Raman_mix=f0_AOM-dF_Raman/2-del;     %[Hz]    45(P) RAMAN   "top"                          45(S) RAMAN   "top"
f2_Raman_mix=f0_AOM+dF_Raman/2;     %[Hz]   -45(S) RAMAN   "horizontal"                 -45(P) RAMAN   "horizonatal"

T_Raman_mix=6.29e-6;
Gs_mod_R_mix=2.0;
phi1_mix=pi;

% % MIX 1
K_R_mix=0.4;

nuller 3 +0.1, nuller 2 -0.1 (2868-2869)
nuller 3 +0.3, nuller 2 -0.3 (2870-2871)

c:\remote\settings201928Oct094851.xml

bragg seq (push magnetically sensative atoms away):c:\remote\settings201928Oct154459.xml
transfer seq (keeps magnetically sensative atoms buts separates them):c:\remote\settings201928Oct154724.xml


transfer to k=-1
    dF_Bragg_1=0.037e6;%0.045e6
    dF_Bragg_2=0.037e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=16.7e-6;%15e-6;
    K_Bragg_src=0.25;%28%0.086;
    Gs_mod_Bragg_src=1.81;%1.8
trasnfer to k=+1
dF_Bragg_1=0.0e6;%0.045e6
    dF_Bragg_2=0.0e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=16.7e-6;%15e-6;
    K_Bragg_src=0.25;%28%0.086;
    Gs_mod_Bragg_src=1.81;%1.8
two sided
dF_Bragg_1=0.017e6;%0.045e6
    dF_Bragg_2=0.017e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;n
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=6.7e-6;%15e-6;
    K_Bragg_src=0.250;%28%0.086;
    Gs_mod_Bragg_src=1.21;%1.8


bragg: c:\remote\settings201929Oct153336.xml
trans: c:\remote\settings201929Oct153457.xml 

jacob: c:\remote\settings201913Nov130202.xml 

2020-07-22 c:\remote\settings202022Jul114948.xml  

attempting to get k=0,-1,-2

best I have is

    dF_Bragg_1=0.108e6;%0.037e6;%0.037e6;%0.045e6
    dF_Bragg_2=0.108e6;%0.037e6;%-0.01e5
    f1_Bragg_src=f0_AOM-dF_Bragg_1;
    f2_Bragg_src=f0_AOM+dF_Bragg_2;
    
    T_Bragg_src=4.2E-6;%16.7e-6;%15e-6;
    K_Bragg_src=0.18;%0.25;%28%0.086;
    Gs_mod_Bragg_src=1.0;%1.9;%2.09;%1.81;%1.8
    