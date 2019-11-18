function m = halo_mode_occupancy(data,opts)
atom_num = data.num_counts./opts.qe; %total number of the acounts in the halo
vr = data.rad; %radius of the halo
bec_width= data.bec_vel_width;%average width of the BEC in velocity space
m = atom_num.*bec_width.^3./(4*pi.*vr.^2.*bec_width);
end