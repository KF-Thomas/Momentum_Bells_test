
%%
% WARNING: use "Run Section" otherwise the output will make no sense... or a mess


%%
clc
close all 
clear all 

%%

clear Nz_results
clear Nzp_results
clear szf_out
clear phase_results

%%  ====    ====    ====    use halo_analysis
%% run halo_analysis everytime on data_folder path change

halo_analysis % run only once  
halo_counts_data = halo{1}.counts_vel';


%% Wrap halo_analysis data 
halo_re_out = halo_reanalysis(halo{1}.counts_vel); 
halo_counts_data = halo_re_out';
halo_reanalysis(halo_re_out);

% halo_reanalysis(halo_re_out)


%%  ====    ====    ====    use simulation
%% Use simulation data
sim_radius = [0.065 0.065*0.05];
sim_shots = 1000;
sim_hits = [15 3];
% sim_detector_noise = [0, 0.0044*0.065];
sim_detector_noise = [0, 0.001];
% sim_detector_noise = [0, 0.01];
% sim_detector_noise = [0, 0];
sim_cut_pi_ang = 0.0*pi;
% sim_cut_pi_ang = 1 * pi;
% sim_cut_pi_ang = 0.25 * pi;
% sim_cut_pi_ang = 0;
sim_cut_pi_per = 0.9;
detector_noise = [0, 0*0.065*0.5];


halo_sim = sim_halo_exp_samples(sim_radius, sim_shots, sim_hits, sim_detector_noise, sim_cut_pi_ang, sim_cut_pi_per);
halo_counts_data = halo_sim.cartesian';


% data_folder_backup = data_folder;
data_folder = '20230213_BEC_mj=0_k=0,-1_helo_new_plate';

halo_reanalysis(halo_sim.cartesian); 


%% 
%% ==== PLOTS ==== 
%% 
%% Calculate stuff

zones_azm = 4;
shift_around = 0.0*pi;
% random_throw_away_perc = 0.1;
random_throw_away_perc = 0.0;
% fprintf(""

%%
if ~exist("Nz_results", "var")
    fprintf("\n    Calculating squeezing_zones ... ");
    Nz_test = [(2:2:50) 60:10:180]'; % faster
    Nz_results = squeezing_zones(halo_counts_data,true, Nz_test,random_throw_away_perc,shift_around); %%%% 
    fprintf("Done --------"+'\n\n');
end 

if ~exist("Nzp_results", "var")
    fprintf("\n    Calculating squeezing_zones ... ");
    % Npz_test = [(2:2:50) 60:10:180]'; % faster
    Nzp_results = squeezing_zones_phav(halo_counts_data, Nz_test, random_throw_away_perc);
    fprintf("Done --------"+'\n\n');
end 

if ~exist("phase_results", "var")
    fprintf("\n\n    Calculating squeezing_phase ... ");
    phase_test = (0:0.01:2)'*pi; % default
    phase_results = squeezing_phase(halo_counts_data, true, phase_test,zones_azm,random_throw_away_perc);
    fprintf("Done --------"+'\n\n');
end 

% if ~exist("szf_out" ,"var")
%     fprintf("\n    Calculating squeezing_zones_filtered ... ");
%     random_throw_away_perc_list = [0, 0.25, 0.50, 0.75];
%     szf_out = squeezing_zones_filtered(halo_counts_data, random_throw_away_perc_list, Nz_test,shift_around); %%%% 
%     fprintf("Done --------"+'\n\n');
% end 

disp("=======    All calculations done    ========    ");

%% Actual plot

squeezing_new(halo_counts_data,true,zones_azm,0,shift_around);
% squeezing_new(halo_counts_data,true,4,0,0*pi);

squeezing_zones_plot(Nz_results) %%%%

squeezing_zones_plot(Nzp_results, 213) %%%%

squeezing_phase_plot(phase_results)

%TODO 
% squeezing_zones_mode(halo_counts_data, true, [(2:2:50) 60:10:180]', [5,20,40], 0);
squeezing_zones_mode(halo_counts_data, true, Nz_test, [5,20,40], 0);

% squeezing_zones_filtered_plot(halo_counts_data, szf_out);

disp("=======    All plotting done    ========    ");

%%  because matlab crashes every 30min ish so... the next few sections to rearrange windows for every relaunch 
%   

fig_name_list = {
    1,  'halo_analysis_radial'
    2,  'halo_analysis_halo'
    3,  'halo_analysis_density'
  200,  'squeezing_new'
  201,  'squeezing_new_halo'
  210,  'squeezing_zones_plot'
  211,  'squeezing_phase'
  213,  'squeezing_zones_phav'
  220,  'squeezing_zones_filtered'
  231,  'halo_reanalysis'
  233,  'halo_reanalysis_fit'
  234,  'halo_reanalysis_azm'
  235,  'halo_reanalysis_azm_elev'
  236,  'halo_reanalysis_wrap_pi'
  237,  'halo_reanalysis_azm_rad'
  238,  'halo_reanalysis_azm_ele'
  239,  'halo_reanalysis_ele_rad'
  250,  'squeezing_zones_mode'
};

% findobj('Type','figure')

%%

for il = 1:size(fig_name_list,1)
    figure(fig_name_list{il,1});
end

%% Laptop screen
for il = 1:size(fig_name_list,1)
    set(figure(fig_name_list{il,1}), 'Position',  [100, 100, 497, 317])
end

movegui(figure(1), [1500, -30]);
movegui(figure(2), [1500,-430]);
movegui(figure(3), [1500,-830]);

movegui(figure(200), [0,-30])
movegui(figure(201), [500,-30]); 
movegui(figure(210), [0,-430]);
movegui(figure(211), [500,-430]);
movegui(figure(213), [1000,-430]);
movegui(figure(220), [0,-830]);
movegui(figure(231), [1000,-30]);
movegui(figure(233), [1000,-430]);
movegui(figure(234), [1000,-830]);
movegui(figure(235), [1300,-930]);
movegui(figure(236), [1300,-930]);
movegui(figure(237), [1300,-930]);
movegui(figure(238), [1300,-930]);
movegui(figure(239), [1300,-930]);


%% bigger monitor config
for il = 1:size(fig_name_list,1)
    set(figure(fig_name_list{il,1}), 'Position',  [100, 100, 597, 417])
end

movegui(figure(1), [1800, -30]);
movegui(figure(2), [1800,-530]);
movegui(figure(3), [1800,-1030]);

movegui(figure(200), [0,-30])
movegui(figure(201), [600,-30]); 
movegui(figure(210), [0,-530]);
movegui(figure(211), [600,-530]);
movegui(figure(220), [0,-1030]);
movegui(figure(231), [1200,-30]);
movegui(figure(233), [1200,-530]);
movegui(figure(234), [1200,-1030]);
movegui(figure(235), [1200,-1530]);

%% vertical monitor config 
for il = 1:size(fig_name_list,1)
    set(figure(fig_name_list{il,1}), 'Position',  [-2000, 2000, 597, 417])
end

movegui(figure(1), [1200,-1030]);
movegui(figure(2), [600,-530]);
movegui(figure(3), [600,-1030]);

movegui(figure(200), [0,-30])
movegui(figure(201), [600,-30]); 
movegui(figure(210), [0,-530]);
movegui(figure(211), [0,-2030]);
movegui(figure(213), [0,-1030]);
movegui(figure(220), [0,-1530]);
movegui(figure(231), [1200,-30]);
movegui(figure(233), [1200,-530]);
movegui(figure(234), [1200,-1530]);
movegui(figure(235), [1200,-2030]);
movegui(figure(236), [600,-1530]);
movegui(figure(237), [600,-2030]);
movegui(figure(238), [600,-2530]);
movegui(figure(239), [1200,-2530]);

% findobj('Type','figure')

%%

% save_figures = true; % this crash my laptop sometimes... I hate matlab 

fig_output_path = "Momentum_Bells_test/Measuring_Quantum_Effeciency/output_plots/";
save_figures = true;

dir_fig = dir(fullfile(fig_output_path, '*.fig'));
files_count = size(dir_fig,1);
found_clash = false;
for id = 1:files_count
    fig_path = fig_output_path + dir_fig(id).name;
    [filepath, name, ~] = fileparts(fig_path);
%     disp(name);
    if strlength(name) > size(data_folder,2)
        name_cut = extractBetween(name, 2, size(data_folder, 2)+1);
%         disp(name_cut)
%         disp(strcmp(name_cut, data_folder))
        if strcmp(name_cut, data_folder)
            found_clash = true; 
        end
    end 
end 
clear id  fig_path name_cut dir_fig files_count name;

% disp("print out section called");
if save_figures==true
% if strcmp('yes',questdlg('Are you sure to export (this crashes the script like 90% of the time)','Continue?','yes','no','no'))
%     questdlg('qstring','title','str1','str2',default)
    pause(0.5); % wait bit more the questdlg to close... otherwise still may crash... I hate matlab
    warning off export_fig:exportgraphics
    fprintf('\n ~~~~ ~~! exporting figures !~~~ (very unstable, I hate matlab) \n');
    fprintf("data_folder = " + data_folder + '\n');

    if found_clash
        fprintf('exist file (prefix) clashed, some files may get overridded, continue? \n');
        pause;
    end
    clear found_clash 

    fprintf('ARE YOU SURE? this thing may crash! do not touch until finish (or matlab crash) \n')
    fprintf('(Press any key to continur or Ctrl-C to excit) \n')
    pause;

    pause(0.5); % wait bit more the questdlg to close... otherwise still may crash... I hate matlab

    last_size_fprintf = 0;
    for il = 1:size(fig_name_list,1)
        fprintf(repmat('\b', 1, last_size_fprintf));
        last_size_fprintf = fprintf("exporting " + num2str(il) + " of " + num2str(size(fig_name_list,1)) + "   ");
        
        fig_id = fig_name_list{il,1};
        fig_name = fig_name_list{il,2};

        figure(fig_id);
        pause(0.1);
        savefig(figure(fig_id),    fig_output_path+'('+data_folder+') '+fig_name+'.fig');     
%         figure(fig_id);
%         pause(0.3);
%         export_fig(figure(fig_id), fig_output_path+'('+data_folder+') '+fig_name+'.pdf');
%         figure(fig_id);
        pause(0.3);
        export_fig(figure(fig_id), fig_output_path+'('+data_folder+') '+fig_name+'.png','-m6');
        
        pause(0.5);
    end 
    clear il fig_id fig_name last_size_fprintf 
    fprintf("\nhorray, all done,  amazing matlab did not crash this time :)))))  \n ");
else 
    disp("skipped");
end 


%     figure(200);
% %     movegui(figure(200), [0,-30]);
%     savefig(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.fig'); 
%     export_fig(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.png');
% 
%     figure(201);
% %     movegui(figure(201), [500,-30]);
%     savefig(figure(201), fig_output_path+'('+data_folder+') '+'squeezing_new_halo.fig'); 
% 
%     figure(210);
% %     movegui(figure(210), [0,-430]);
%     savefig(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.fig'); 
%     
%     figure(211);
% %     movegui(figure(211), [500,-430]);
%     savefig(figure(211), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot2.fig'); 
%     
%     figure(220);
% %     movegui(figure(220), [0,-830]);
%     savefig(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered.fig'); 
%     
%     figure(231);
% %     movegui(figure(231), [1000,-30]);
%     savefig(figure(231), fig_output_path+'('+data_folder+') '+'halo_reanalysis.fig'); 
%     
%     figure(233);
% %     movegui(figure(233), [1000,-430]);
%     savefig(figure(232), fig_output_path+'('+data_folder+') '+'halo_reanalysis2.fig'); 


    
%     exportgraphics(figure(231),'test.tiff'); % worked once, then crashed
%     exportgraphics(figure(220),'test.png'); % ccrash
%     exportgraphics(figure(220),'test.png','Resolution',300); % crash
%     exportgraphics(figure(220),'test.pdf')
%       set(gcf,'renderer','zbuffer');
%     %  print(figure(200), 'test.png', '-dpng', '-r300'); % don't know... this crashes on macOS every time
%     
%     saveas(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.pdf'); 
%     saveas(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.pdf'); 
%     saveas(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered_plot.pdf'); 
%  
% %     saveas(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.png'); 
% %     saveas(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.png'); 
% %     saveas(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered_plot.png'); 
    



%%



%%
























