function bec = halo_cent(data,opts_cent)
%find the center of the top middle and bottom bec
%also calculates the transfer fraction
opts_cent.crop = [opts_cent.t_bounds{1}; -0.03, 0.03; -0.03, 0.03];
[bec.centre_top,bec.width_top,bec.counts_top,bec.centre_OK_top] =  find_dist_cen(data,opts_cent);
opts_cent.crop = [opts_cent.t_bounds{2}; -0.03, 0.03; -0.03, 0.03];
[bec.centre_mid,bec.width_mid,bec.counts_mid,bec.centre_OK_mid] =  find_dist_cen(data,opts_cent);
opts_cent.crop = [opts_cent.t_bounds{3}; -0.03, 0.03; -0.03, 0.03];
[bec.centre_btm,bec.width_btm,bec.counts_btm,bec.centre_OK_btm] =  find_dist_cen(data,opts_cent);
bec.trans_top = bec.counts_top./(bec.counts_top+bec.counts_mid+bec.counts_btm);
bec.trans_mid = bec.counts_mid./(bec.counts_top+bec.counts_mid+bec.counts_btm);
bec.trans_btm = bec.counts_btm./(bec.counts_top+bec.counts_mid+bec.counts_btm);
end