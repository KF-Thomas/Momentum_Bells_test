%% monitor the magnetic transfer
clear all
this_folder = fileparts(fileparts(which(mfilename)));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=true;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


% anal_opts.atom_laser.pulsedt=8.000e-3;
% anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
% anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
% anal_opts.atom_laser.pulses=100;
%
% anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
% anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.95;
%
% anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=2.3;
anal_opts.trig_ai_in=20;


% anal_opts.osc_fit.binsx=1000;
% anal_opts.osc_fit.blur=1;
% anal_opts.osc_fit.xlim=[-20,20]*1e-3;
% anal_opts.osc_fit.tlim=[0.86,1.08];
% anal_opts.osc_fit.dimesion=2; %Sel ect coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=600;

hebec_constants
combined_struct = @(S,T) cell2struct(cellfun(@vertcat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);
const.fall_distance = 8.52925545e-01;


% END USER VAR-----------------------------------------------------------
fclose('all')
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

hebec_constants
anal_opts.tdc_import.mat_save=false;
anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;

if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end

%anal_out.dir=sprintf('%sout\\monitor\\',...
%    anal_opts.tdc_import.dir);

anal_out.dir=[fullfile(anal_opts.tdc_import.dir,'out','monitor'),filesep];
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;

frac_opts.num_lim = 1e3;
frac_opts.transfer_state = 'momentum';
frac_opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds


%%
% opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
 opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
log_folder = 'log_Phi.txt';
data_folders = {
%     ''
%     '20210524_bec_MZ_10'
    'full_interferometer\mach-zender\becs\20210520_bec_MZ_4'
%     '20210521_bec_MZ_7'
%     '20210519_bec_MZ_2'
    };
%%

loop_num=0;
opts.num_lim = 1e3;%
phi_vec= [];
out_data = {};
opts.import.force_reimport = true;
for folder_indx = 1:length(data_folders)
    %% import raw data
    data_folder = data_folders{folder_indx};
    opts.import.dir = fullfile(opts.data_root, data_folder);
    opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');
    opts.logfile = fullfile(opts.import.dir,log_folder);
    
    %% import raw data
    [data, ~] = import_mcp_tdc_data(opts.import);
    
    %% remove any ringing
    opts.ring_lim = 0.09e-6;%0.1e-6;%-1;%0;%0.101 %how close can points be in time
    data_masked = ring_removal(data,opts.ring_lim);
    
    %% set up relevant constants
    hebec_constants
    %%
    phi_logs = table2array(readtable(opts.logfile));
    l_t=phi_logs(:,2);
    d_t=data.time_create_write(:,1);
    phi_log_matched = zeros(length(d_t),1).*nan;
    phi_log_check = zeros(length(d_t),1);
    for ii = 1: length(l_t)
        phi_c = phi_logs(ii,3);
        l_c=l_t(ii);
        t_mask=l_c+17<d_t & l_c+29.5>d_t;
        d_indx=find(t_mask);
        phi_log_matched(d_indx,1) = phi_c;
        phi_log_check(d_indx,1) = 1;
        
    end
    num_check = data_masked.num_counts>opts.num_lim;
    data_masked_halo = struct_mask(data_masked,num_check & phi_log_check');
    phi_logs_masked = phi_log_matched(num_check& phi_log_check',1);
    out_frac = fraction_calc(data_masked_halo,frac_opts);
    
    N_mask = out_frac.Ntotal>50;
    
    out_frac = struct_mask(out_frac,N_mask);
    phi_logs_masked = phi_logs_masked(N_mask,:);
    
    unique_phi = unique(phi_logs_masked); %unique phases used in this scan
    for ii = 1:length(unique_phi)
        current_phi = unique_phi(ii,1);
        phi_mask = (phi_logs_masked ==current_phi);
        frac_masked = out_frac.fracs(phi_mask,2);
        if ~ismember(current_phi,phi_vec)
            phi_vec = [phi_vec,current_phi];
        end
        phi_indx = find(phi_vec==current_phi);
        if length(out_data)<phi_indx
            out_data{phi_indx} = frac_masked;
        else
            out_data{phi_indx} = [out_data{phi_indx};frac_masked];
        end
    end
    
    
    %                 cost = momentum_transfer_cost(out_frac.shot_num');
    
    %                 mag_history.all_shots=[mag_history.all_shots,batch_data.mcp_tdc.shot_num];
    %                 mag_history.shot_num=[mag_history.shot_num,out_frac.shot_num'];
    %                 mag_history.trans_frac=[mag_history.trans_frac;out_frac.fracs];
    %                 mag_history.Ns=[mag_history.Ns;out_frac.Ns];
    %                 mag_history.cost =[mag_history.cost;cost.val];
    
    
    
end
%%
for ii = 1:length(out_data)
    k0_mean(ii) = mean(out_data{ii});
    k0_std(ii) = std(out_data{ii});
end

xp = linspace(min(phi_vec),max(phi_vec));
% fit = @(b,x)  b(1).*cos(x.*b(2) + 2*pi/b(6)).*(cos(x.*b(5) + 2*pi/b(3))) + b(4);    % Function to fit
% fit = @(b,x)  b(1).*cos(x.*b(2) + b(3)) + b(4);    % Function to fit [1.229,1,0.8088,0.906]
fit = @(b,x)  b(1).*cos(x + b(2)) + b(3);    % Function to fit
% [1.229,0.8088,0.906]
best_fit = fitnlm(phi_vec,k0_mean,fit,[1,pi,1],'CoefficientNames',{'Amp','Phase','Offset'}); %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
[ysamp_val,ysamp_ci]=predict(best_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'


stfig('Momentum Transfer Fraction against phase');
clf
%                 plot(phi_vec,...
%                     mag_history.trans_frac(:,1:3)',...
%                     'k.')
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
errorbar(phi_vec,k0_mean,k0_std,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
hold on
    plot(xp,ysamp_val,'r','LineWidth',1.5)
    drawnow
    yl=ylim*1.1;
    plot(xp,ysamp_ci,'color',[1,1,1].*0.5)

grid on
h=gca;
grid on    % turn on major grid lines
grid minor % turn on minor grid lines
% Set limits and grid spacing separately for the two directions:
% Must set major grid line properties for both directions simultaneously:
h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
h.GridAlpha=1;  % the default is partially transparent
h.GridColor=[0,0,0]; % here's the color for the major grid lines
% Idem for minor grid line properties:
h.MinorGridLineStyle='-';
h.MinorGridAlpha=0.1;
h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
xlabel('Final Beam Splitter Phase')
ylabel('Tranfer Fraction')
%                 legend('$k=+1$','$k=0$','$k=-1$')
%                 legend('$k=-2$','$k=-1$','$k=0$')

