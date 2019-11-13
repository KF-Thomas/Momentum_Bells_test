
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20191101_brief_halo_data';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=true;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=100;

anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.95;

anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


anal_opts.osc_fit.binsx=1000;
anal_opts.osc_fit.blur=1;
anal_opts.osc_fit.xlim=[-20,20]*1e-3;
anal_opts.osc_fit.tlim=[0.86,1.08];
anal_opts.osc_fit.dimesion=2; %Sel ect coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=50;

hebec_constants
%% find centers
opts.cent.visual = 0;
opts.cent.bin_size = 0.5e-5 * [1, 1, 1];
opts.cent.threshold = 2.5;
opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972]}; %time bounds for the different momentum states k=+1,0,-1 respectively

opts.vel_conv.visual = 0;
opts.vel_conv.title = 'top halo';
% top halo
trap_switch_off= 0;
opts.vel_conv.const.g0 = const.g0;
opts.vel_conv.const.fall_distance = const.fall_distance;
opts.vel_conv.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.v_mask=[0.8,1.2]; %bounds on radisu as multiple of radius value
opts.vel_conv.z_mask = [-0.04,0.04];

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
halo_history.halo_N=[];
halo_history.corr_length=[];
halo_history.trans_frac=[];
halo_history.shot_num=[];
stfig('Halo Num History')
clf;
stfig('Transfer Fraction History')
clf

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
        batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
        %data.mcp_tdc=mcp_tdc_data;
        %just to give me a logical vector
        batch_data.mcp_tdc.all_ok=batch_data.mcp_tdc.num_counts>1e3;
        batch_data.mcp_tdc.all_ok(batch_data.mcp_tdc.all_ok)=...
            cellfun(@(x) x(end,1),batch_data.mcp_tdc.counts_txy(batch_data.mcp_tdc.all_ok))>anal_opts.dld_aquire*0.8;
        if sum(batch_data.mcp_tdc.all_ok)==0
            fprintf('waiting for file to be writen\n')
            pause(0.1)
        else
bec = halo_cent(batch_data.mcp_tdc,opts.cent);
trans_frac = [bec.trans_top';bec.trans_mid';bec.trans_btm'];
%% convert data to velocity
opts.vel_conv.bec_center.north = bec.centre_top;
opts.vel_conv.bec_center.south = bec.centre_mid;
opts.vel_conv.bec_width.north = bec.width_top;
opts.vel_conv.bec_width.south = bec.width_mid;
%%
top_halo = halo_vel_conv(batch_data.mcp_tdc,opts.vel_conv);

            halo_history.shot_num=[halo_history.shot_num,anal_opts.tdc_import.shot_num];
            halo_history.halo_N=[halo_history.halo_N,top_halo.num_counts];
            halo_history.trans_frac=[halo_history.trans_frac,trans_frac];

            %trim the history vectors
            if numel(halo_history.shot_num)>anal_opts.history.shots
                %bit sloppy but will assume they are the same length
                halo_history.shot_num=halo_history.shot_num(end-anal_opts.history.shots:end);
                halo_history.halo_N=halo_history.halo_N(end-anal_opts.history.shots:end);
                halo_history.trans_frac=halo_history.trans_frac(:,end-anal_opts.history.shots:end);
            end

%             sfigure(1);
%             scatter(halo_history.shot_num,...
%                 halo_history.halo_N,...
%                 'kx')
            stfig('Halo Num History')
            plot(halo_history.shot_num,...
                halo_history.halo_N,...
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
            ylabel('N in Halo')
            
            stfig('Transfer Fraction History')
            plot(halo_history.shot_num,...
                halo_history.trans_frac,...
                'LineWidth',1.5)
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
            ylabel('Transfer Fraction')
            legend('top','mid','btm')
            pause(1e-6)
%             saveas(gcf,fullfile(anal_out.dir,'freq_history.png'))

        end
    end

end