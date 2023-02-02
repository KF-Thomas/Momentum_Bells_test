
% Tip: use "Run Section"
save_figures = false; % this crash my laptop sometimes... I hate matlab 

fig_output_path = "Momentum_Bells_test/Measuring_Quantum_Effeciency/output_plots/";


%% run halo_analysis everytime on data_folder path change
% halo_analysis % run only once  

if ~exist("halo" ,"var")
    halo_analysis
end


%%
squeezing_new(halo{1}.counts_vel',true,10);
movegui(figure(200), [0,-30]);
movegui(figure(210), [600,-30]);



if ~exist("Nz_results", "var")
    Nz_test = [(2:2:50) 60:10:180]'; % faster
    Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test); %%%% 
end 
squeezing_zones_plot(Nz_results) %%%%
movegui(figure(210), [0,-530]);




if ~exist("szf_out" ,"var")
    random_throw_away_perc_list = 0:0.3:0.9; %%%% 
    szf_out = squeezing_zones_filtered(halo{1}.counts_vel', [0, 0.25, 0.50, 0.75], Nz_test); %%%% 
end 
squeezing_zones_filtered_plot(halo{1}.counts_vel', szf_out);
movegui(figure(220), [0,-1030]);



%%
if save_figures
    savefig(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.fig'); 
    savefig(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.fig'); 
    savefig(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered_plot.fig'); 
    
    
    %  print(figure(200), 'test.png', '-dpng', '-r300'); % don't know... this crashes on macOS every time
    
    saveas(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.pdf'); 
    saveas(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.pdf'); 
    saveas(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered_plot.pdf'); 
 
%     saveas(figure(200), fig_output_path+'('+data_folder+') '+'squeezing_new.png'); 
%     saveas(figure(210), fig_output_path+'('+data_folder+') '+'squeezing_zones_plot.png'); 
%     saveas(figure(220), fig_output_path+'('+data_folder+') '+'squeezing_zones_filtered_plot.png'); 
    
end 