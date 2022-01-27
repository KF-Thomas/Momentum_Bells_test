% short function to monitor approximate number in halos
this_folder = fileparts(fileparts(which(mfilename)));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20200803_early_k=0,-1,-2_halo_data\';
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
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


% anal_opts.osc_fit.binsx=1000;
% anal_opts.osc_fit.blur=1;
% anal_opts.osc_fit.xlim=[-20,20]*1e-3;
% anal_opts.osc_fit.tlim=[0.86,1.08];
% anal_opts.osc_fit.dimesion=2; %Sel ect coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=200;

hebec_constants
const.fall_distance = 8.52925545e-01;
%% find centers
opts.cent.visual = 2;
opts.cent.threshold = [100,30,30].*1e3; %set in inverse units (Hz for time 1/m for space)
opts.cent.sigma = [8e-5,25e-5,25e-5];
% opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972],[3.8,3.95]}; %time bounds for the different momentum states k=+1,0,-1 respectively
opts.cent.t_bounds = {[1.741,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};;

opts.vel_conv.plot_percentage = 0.2;
opts.vel_conv.visual = 0;
opts.vel_conv.title = 'top halo';
% top halo
trap_switch_off= 0;
opts.vel_conv.const.g0 = const.g0;
opts.vel_conv.const.fall_distance = const.fall_distance;
opts.vel_conv.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.v_mask=[0.8,1.2]; %bounds on radisu as multiple of radius value
opts.vel_conv.z_mask = [-0.04,0.04];

opts.num_lim = 2e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = 10;%0;% %minimum allowed number in halo 10

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



%%

halo_history.corr_length=[];
halo_history.trans_frac=[];
halo_history.shot_num=[];

halo_history.top.halo_N=[];
halo_history.top.m=[];
halo_history.top.g2=[];

halo_history.btm.halo_N=[];
halo_history.btm.m=[];
halo_history.btm.g2=[];

halo_history.shot_num_masked = [];
stfig('Halo Num History');
clf;
stfig('Ratio of Halo Number');
clf;
% stfig('Transfer Fraction History');
% clf
% stfig('Mode Occupancy and Correlation Amplitude');
% clf

%%

loop_num=0;
while true
    pause(0.1)
    batch_data=[];
    batch_data.shot_num=[];
    anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
    
    max_shot_num=max(anal_opts.tdc_import.shot_num);
    anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
        anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));
    %remove processed ones
    
    anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
        ~ismember(anal_opts.tdc_import.shot_num, halo_history.shot_num ) );
    
    if numel(anal_opts.tdc_import.shot_num)==0
        if mod(loop_num,4)==0
            pause(.2)
            fprintf('\b\b\b')
            loop_num=1;
        else
            pause(.1) %little wait animation
            fprintf('.')
            loop_num=loop_num+1;
        end
    else
        try
            batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
            %just to give me a logical vector
            batch_data.mcp_tdc.all_ok=batch_data.mcp_tdc.num_counts>1e3;
%             batch_data.mcp_tdc.all_ok(batch_data.mcp_tdc.all_ok)=...
%                 cellfun(@(x) x(end,1),batch_data.mcp_tdc.counts_txy(batch_data.mcp_tdc.all_ok))>anal_opts.dld_aquire*0.8;
            if sum(batch_data.mcp_tdc.all_ok)==0
                fprintf('waiting for file to be writen\n')
                pause(1.0)
            else
                num_check = batch_data.mcp_tdc.num_counts>opts.num_lim;
                num_shots = size(batch_data.mcp_tdc.num_counts,2);
                masked_data = hotspot_mask(batch_data.mcp_tdc);
                top_halo_num_counts = zeros(1,num_shots);
                btm_halo_num_counts = zeros(1,num_shots);
                for ii = 1:num_shots
                    this_shot = masked_data.counts_txy{ii};
%                     top_halo_num_counts(ii) = size(masktxy_square(this_shot, [3.866, 3.874; -0.03, 0.03; -0.03, 0.03]),1);
%                     btm_halo_num_counts(ii) = size(masktxy_square(this_shot, [3.854, 3.861; -0.03, 0.03; -0.03, 0.03]),1);
                    
                    top_halo_num_counts(ii) = size(masktxy_square(this_shot, [1.76, 1.7665; -0.03, 0.03; -0.03, 0.03]),1);
                    btm_halo_num_counts(ii) = size(masktxy_square(this_shot, [1.747, 1.754; -0.03, 0.03; -0.03, 0.03]),1);
                end
                
                halo_history.shot_num=[halo_history.shot_num,anal_opts.tdc_import.shot_num];
                
                halo_history.btm.halo_N=[halo_history.btm.halo_N,btm_halo_num_counts];
                
                halo_history.top.halo_N=[halo_history.top.halo_N,top_halo_num_counts];
                
                %trim the history vectors
                if numel(halo_history.shot_num)>anal_opts.history.shots
                    %bit sloppy but will assume they are the same length
                    halo_history.shot_num=halo_history.shot_num(end-anal_opts.history.shots:end);
                    halo_history.top.halo_N=halo_history.top.halo_N(end-anal_opts.history.shots:end);
                    halo_history.btm.halo_N=halo_history.btm.halo_N(end-anal_opts.history.shots:end);
                end
                
                stfig('Halo Num History');
                plot(halo_history.shot_num,...
                    halo_history.top.halo_N,...
                    'kx-','LineWidth',1.5)
                hold on
                plot(halo_history.shot_num,...
                    halo_history.btm.halo_N,...
                    'bx-','LineWidth',1.5)
                hold off
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
                xlabel('Shot Number')
                ylabel('N in Halo')
                legend('top halo','bottom halo')
                
                stfig('Ratio of Halo Number')
                plot(halo_history.shot_num,...
                    halo_history.btm.halo_N./halo_history.top.halo_N,...
                    'kx-','LineWidth',1.5)
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
                xlabel('Shot Number')
                ylabel('Ratio')
                pause(1e-6)
                %             saveas(gcf,fullfile(anal_out.dir,'freq_history.png'))
                
            end
        catch
            pause(1.0)
        end
    end
    
end