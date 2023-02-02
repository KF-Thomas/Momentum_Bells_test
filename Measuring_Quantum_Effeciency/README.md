# Steps to measure quantum efficiency



## 1. `halo_analysis.m`

File is in `/Momentum_Bells_test/halo_analysis.m`.

Need to set `opts.data_root` and `data_folder`, then run 

* might takes a while to run, need to wait until warning(?) and plots to show up

## 2. `squeezing_new`

This generates Fig.2 as https://link.aps.org/doi/10.1103/PhysRevLett.105.190402 [1]

```matlab
squeezing_new(halo{1}.counts_vel',true,10);
```

## 3. `squeezing_zones`

Fig 3 as [1]

```matlab
Nz_test = [(2:2:50) 60:10:360]'; % default
Nz_test = [(2:2:50) 60:10:180]'; % faster
Nz_results = squeezing_zones(halo{1}.counts_vel',true);
Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test); %%%% 
```

![squeezing_zones_plot (on data 20221212_new_plates_halo_test_34_combined_notHe34_just_test_34)](/Users/tonyyan/Documents/_ANU/_He_BEC_Group/Measuring_Quantum_Effeciency/output_plots/squeezing_zones_plot (on data 20221212_new_plates_halo_test_34_combined_notHe34_just_test_34).png)



###  `squeezing_zones_plot`

```matlab
squeezing_zones_plot(Nz_results) %%%%
squeezing_zones_plot(squeezing_zones(halo{1}.counts_vel',false))
squeezing_zones_plot(squeezing_zones(halo{1}.counts_vel',false, Nz_test)) 
squeezing_zones_plot(squeezing_zones(halo{1}.counts_vel',false, Nz_test)) 

```

### `squeezing_zones_filtered`

```
random_throw_away_perc_list = 0:0.3:0.9; %%%% 
szf_out = squeezing_zones_filtered(halo{1}.counts_vel');
szf_out = squeezing_zones_filtered(halo{1}.counts_vel', random_throw_away_perc_list, Nz_test);
szf_out = squeezing_zones_filtered(halo{1}.counts_vel', [0, 0.25, 0.50, 0.75], Nz_test); %%%% 

squeezing_zones_filtered_plot(halo{1}.counts_vel', szf_out);
```



![squeezing_zones_filtered_plot (on data 20221212_new_plates_halo_test_34_combined_notHe34_just_test_34)](/Users/tonyyan/Documents/_ANU/_He_BEC_Group/Measuring_Quantum_Effeciency/output_plots/squeezing_zones_filtered_plot (on data 20221212_new_plates_halo_test_34_combined_notHe34_just_test_34).png)











