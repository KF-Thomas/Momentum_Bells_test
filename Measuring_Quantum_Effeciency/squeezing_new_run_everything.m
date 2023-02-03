
% Tip: use "Run Section"
save_figures = false; % this crash my laptop sometimes... I hate matlab 

fig_output_path = "Momentum_Bells_test/Measuring_Quantum_Effeciency/output_plots/";

%%
% clear halo 
clear Nz_results
clear szf_out
clear phase_results
clear halo_counts_data


%% run halo_analysis everytime on data_folder path change
halo_analysis % run only once  

% if ~exist("halo" ,"var")
%     halo_analysis;
%     halo_counts_data = halo{1}.counts_vel';
% 
% end
halo_counts_data = halo{1}.counts_vel';


%%
halo_re_out = halo_reanalysis(halo{1}.counts_vel); 
halo_counts_data = halo_re_out';


%% Use simulation data
sim_radius = [0.065 0.065*0.05];
sim_shots = 100;
sim_hits = [20 0.0001];
halo_sim = sim_halo_exp_samples(sim_radius, sim_shots, sim_hits);
halo_counts_data = halo_sim.cartesian';
detector_noise = [0, 0.065*0.5];

data_folder_backup = data_folder;
data_folder = 'sim_halo_exp_samples_2';

%%
zones_azm = 14;
shift_around = 0.21991;

squeezing_new(halo_counts_data,true,zones_azm,0,shift_around);

movegui(figure(200), [0,-30])
movegui(figure(201), [600,-30]);


if ~exist("Nz_results", "var")
    Nz_test = [(2:2:50) 60:10:180]'; % faster
    Nz_results = squeezing_zones(halo_counts_data,true, Nz_test,0,shift_around); %%%% 
end 
squeezing_zones_plot(Nz_results) %%%%
movegui(figure(210), [0,-430]);


if ~exist("phase_results", "var")
    phase_test = (0:0.01:2)'*pi; % default
    phase_results = squeezing_phase(halo_counts_data, true, phase_test,zones_azm,0);
end 
squeezing_phase_plot(phase_results)
movegui(figure(211), [600,-430]);



if ~exist("szf_out" ,"var")
    random_throw_away_perc_list = [0, 0.25, 0.50, 0.75];
    szf_out = squeezing_zones_filtered(halo_counts_data, [0, 0.25, 0.50, 0.75], Nz_test,shift_around); %%%% 
end 
squeezing_zones_filtered_plot(halo_counts_data, szf_out);
movegui(figure(220), [0,-830]);



%%
if save_figures==true
    savefig(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.fig'); 
    savefig(figure(201), fig_output_path+'('+data_folder+') '+'squeezing_new_halo.fig'); 
    savefig(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.fig'); 
    savefig(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered_plot.fig'); 
    
    
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

close figure(200)








end 