function halo_history = halo_num_check(import_dir,num_shots,shot_mask)
% short function to monitor approximate number in halos
this_folder = fileparts(fileparts(which(mfilename)));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir=import_dir;%'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
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

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=19.5;
anal_opts.dld_aquire=2;


anal_opts.history.shots=num_shots;%50;

% hebec_constants
% const.fall_distance = 8.52925545e-01;
%% find centers
opts.cent.visual = 0;
opts.cent.threshold = [100,30,30].*1e3; %set in inverse units (Hz for time 1/m for space)
opts.cent.sigma = [8e-5,25e-5,25e-5];
opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.884,3.896],[3.75,4]};%
% opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972],[3.8,3.95]}; %time bounds for the different momentum states k=+1,0,-1 respectively
% opts.cent.t_bounds = {[1.741,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};;

opts.vel_conv.plot_percentage = 0.2;
opts.vel_conv.visual = 0;
opts.vel_conv.title = 'top halo';
% top halo
trap_switch_off= 0;
% opts.vel_conv.const.g0 = const.g0;
% opts.vel_conv.const.fall_distance = const.fall_distance;
opts.vel_conv.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.v_mask=[0.8,1.2]; %bounds on radisu as multiple of radius value
opts.vel_conv.z_mask = [-0.04,0.04];

opts.num_lim = 2e3;%0.5e3;% %minimum atom number 1.5e3
% opts.halo_N_lim = 10;%0;% %minimum allowed number in halo 10

% END USER VAR-----------------------------------------------------------
fclose('all')
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% hebec_constants
anal_opts.tdc_import.mat_save=false;
% anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;



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


%%


not_done = 1;
while not_done
    try
        batch_data=[];
        batch_data.shot_num=[];
        anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
        if nargin>2
            if length(shot_mask)<length(anal_opts.tdc_import.shot_num)
                shot_mask = logical([zeros(1,length(anal_opts.tdc_import.shot_num)-length(shot_mask)),shot_mask]);
            end
            anal_opts.tdc_import.shot_num = anal_opts.tdc_import.shot_num(logical(shot_mask));
        end
        
        max_shot_num=max(anal_opts.tdc_import.shot_num);
        anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
            anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));
        
        
        %remove processed ones
        
        anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
            ~ismember(anal_opts.tdc_import.shot_num, halo_history.shot_num ));
        
        batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
        %just to give me a logical vector
        batch_data.mcp_tdc.all_ok=batch_data.mcp_tdc.num_counts>1e3;
        num_check = batch_data.mcp_tdc.num_counts>opts.num_lim;
        num_shots = size(batch_data.mcp_tdc.num_counts,2);
        masked_data = hotspot_mask(batch_data.mcp_tdc);
        top_halo_num_counts = zeros(1,num_shots);
        btm_halo_num_counts = zeros(1,num_shots);
        for ii = 1:num_shots
            this_shot = masked_data.counts_txy{ii};
%             top_halo_num_counts(ii) = size(masktxy_square(this_shot, [1.76, 1.7665; -0.03, 0.03; -0.03, 0.03]),1);
%             btm_halo_num_counts(ii) = size(masktxy_square(this_shot, [1.747, 1.754; -0.03, 0.03; -0.03, 0.03]),1);

            top_halo_num_counts(ii) = size(masktxy_square(this_shot, [3.881, 3.886; -0.03, 0.03; -0.03, 0.03]),1);
            btm_halo_num_counts(ii) = size(masktxy_square(this_shot, [3.867, 3.872; -0.03, 0.03; -0.03, 0.03]),1);
        end
        not_done = 0;
    catch
        pause(1.0)
    end
end


halo_history.shot_num=[halo_history.shot_num,anal_opts.tdc_import.shot_num];

halo_history.btm.halo_N=[halo_history.btm.halo_N,btm_halo_num_counts];

halo_history.top.halo_N=[halo_history.top.halo_N,top_halo_num_counts];

end
