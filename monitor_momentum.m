%% monitor the magnetic transfer

this_folder = fileparts(fileparts(which(mfilename)));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
% BEGIN USER VAR-------------------------------------------------
% anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20210727_bragg_amp_scan_new_trap_2';
% 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20210724_bragg_amp_scan_new_trap';
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20210726_bragg_width_scan_new_trap_3';
% 

% 
% anal_opts.tdc_import.dir='Z:\EXPERIMENT-DATA\2020_Momentum_Bells\pulse_characterisation\20210430_bragg_amp_scans\20210511_bragg_amp_scan_8\';
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

anal_opts.history.shots=50;

hebec_constants
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

% if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end
% if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end

%anal_out.dir=sprintf('%sout\\monitor\\',...
%    anal_opts.tdc_import.dir);

anal_out.dir=[fullfile(anal_opts.tdc_import.dir,'out','monitor'),filesep];
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;

frac_opts.num_lim = 1e3;
frac_opts.transfer_state = 'momentum';
frac_opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds

%%
mag_history.trans_frac=[];
mag_history.shot_num=[];
mag_history.all_shots=[];
mag_history.cost=[];
mag_history.Ns=[];
stfig('Momentum Transfer Fraction History');
clf;
% hold on
stfig('Momentum Transfer Cost History');
clf;
stfig('External Fraction');
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
        ~ismember(anal_opts.tdc_import.shot_num, mag_history.all_shots ) );
    
    if numel(anal_opts.tdc_import.shot_num)==0
        if mod(loop_num,4)==0
            pause(.7)
            fprintf('\b\b\b')
            loop_num=1;
        else
            pause(.5) %little wait animation
            fprintf('.')
            loop_num=loop_num+1;
        end
    else
        try
            batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
            %just to give me a logical vector
            batch_data.mcp_tdc.all_ok=batch_data.mcp_tdc.num_counts>1e3;
            batch_data.mcp_tdc.all_ok(batch_data.mcp_tdc.all_ok)=...
                cellfun(@(x) x(end,1),batch_data.mcp_tdc.counts_txy(batch_data.mcp_tdc.all_ok))>anal_opts.dld_aquire*0.8;
            if sum(batch_data.mcp_tdc.all_ok)==0
                fprintf('waiting for file to be writen\n')
                pause(1.0)
            else
                out_frac = fraction_calc(batch_data.mcp_tdc,frac_opts);
%                 cost = momentum_transfer_cost(out_frac.shot_num',anal_opts.tdc_import.dir);
                
                mag_history.all_shots=[mag_history.all_shots,batch_data.mcp_tdc.shot_num];
                mag_history.shot_num=[mag_history.shot_num,out_frac.shot_num'];
                mag_history.trans_frac=[mag_history.trans_frac;out_frac.fracs];
                mag_history.Ns=[mag_history.Ns;out_frac.Ns];
%                 mag_history.cost =[mag_history.cost;cost.val];
                
                %trim the history vectors
                if numel(mag_history.shot_num)>anal_opts.history.shots
                    %bit sloppy but will assume they are the same length
                    mag_history.shot_num=mag_history.shot_num(end-anal_opts.history.shots:end);
                    mag_history.trans_frac=mag_history.trans_frac(end-anal_opts.history.shots:end,:);
%                     mag_history.cost=mag_history.cost(end-anal_opts.history.shots:end,:);
                    mag_history.Ns=mag_history.Ns(end-anal_opts.history.shots:end,:);
                end
                
                stfig('Momentum Transfer Fraction History');
                plot(mag_history.shot_num,...
                    mag_history.trans_frac(:,1:3)',...
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
                ylabel('Tranfer Fraction')
                legend('$k=+1$','$k=0$','$k=-1$')
%                 legend('$k=-2$','$k=-1$','$k=0$')
                
                pause(0.1)
                %             saveas(gcf,fullfile(anal_out.dir,'freq_history.png'))
%                 stfig('Momentum Transfer Cost History');
%                 plot(mag_history.shot_num,...
%                     mag_history.cost',...
%                     'LineWidth',1.5)
%                 hold on
%                 scatter(mag_history.shot_num,...
%                     mag_history.cost',...
%                     'LineWidth',1.5)
%                 hold off 
%                 grid on
%                 h=gca;
%                 grid on    % turn on major grid lines
%                 grid minor % turn on minor grid lines
                % Set limits and grid spacing separately for the two directions:
                % Must set major grid line properties for both directions simultaneously:
%                 h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
%                 h.GridAlpha=1;  % the default is partially transparent
%                 h.GridColor=[0,0,0]; % here's the color for the major grid lines
%                 % Idem for minor grid line properties:
%                 h.MinorGridLineStyle='-';
%                 h.MinorGridAlpha=0.1;
%                 h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
%                 xlabel('Shot Number')
%                 ylabel('Tranfer Cost Function')
%                 legend('cost')
                
                stfig('External Fraction');
                num_mask = ~isnan(mag_history.shot_num)';
% %                 plot(mag_history.shot_num(num_mask),1-1./mag_history.trans_frac(num_mask,4),'LineWidth',1.5)
% %                 hold on
% %                 scatter(mag_history.shot_num(num_mask),1-1./mag_history.trans_frac(num_mask,4),'LineWidth',1.5)
%                 plot(mag_history.shot_num(num_mask),1-1./mag_history.trans_frac(num_mask,4),'LineWidth',1.5)
%                 hold on
%                 scatter(mag_history.shot_num(num_mask),1-1./mag_history.trans_frac(num_mask,4),'LineWidth',1.5)
                %external to k=0,-1 states
                plot(mag_history.shot_num(num_mask),(mag_history.Ns(num_mask,4)-mag_history.Ns(num_mask,2)-mag_history.Ns(num_mask,3))./mag_history.Ns(num_mask,4),'LineWidth',1.5)
                hold on
                scatter(mag_history.shot_num(num_mask),(mag_history.Ns(num_mask,4)-mag_history.Ns(num_mask,2)-mag_history.Ns(num_mask,3))./mag_history.Ns(num_mask,4),'LineWidth',1.5)
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
                ylabel('External Fraction')
                
                pause(0.1)
            end
        catch
            pause(1.0)
        end
    end
    
end