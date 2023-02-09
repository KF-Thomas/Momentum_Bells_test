function out=create_bins(bin_info)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%creates bins for squeezing
%inputs 
    %bin_info ; struct containing specifications needed for binning

%outputs
    %bin_centers ; contains center of each bin
    %bin_pairs ; contains start and end of each bin
%--------------------------------------------------------------------------
range_azm=bin_info.range_azm;
range_elev=bin_info.range_elev;
steps_azm = bin_info.steps_azm;
steps_elev= bin_info.steps_elev;
bin_width_azm = bin_info.bin_width_azm;
bin_width_elev= bin_info.bin_width_elev;
mirror_azm=bin_info.mirror_azm; %use for the mag sensitive halo which you must cut out a region ??

% create the bins
bin_centers_azm=wrapTo2Pi(linspace(range_azm(1),range_azm(2),steps_azm+1)-range(range_azm)/(2*steps_azm));
bin_centers_azm=bin_centers_azm(2:end); %first and last elements are separated by 2pi and we only need one of them

%pairs containing start and end of each bin
bin_pairs_azm=transpose([wrapTo2Pi(bin_centers_azm-bin_width_azm/2) ; wrapTo2Pi(bin_centers_azm+bin_width_azm/2)]);


if mirror_azm % ? not sure what this does ?
    bin_centers_azm=[bin_centers_azm, bin_centers_azm+pi];
    bin_pairs_azm=[bin_pairs_azm;bin_pairs_azm+pi];
end

bin_centers_elev=wrapTo2Pi(linspace(range_elev(1),range_elev(2),steps_elev+1)-range(range_elev)/(2*steps_elev));
bin_centers_elev=bin_centers_elev(2:end);
bin_pairs_elev=transpose([wrapTo2Pi(bin_centers_elev-bin_width_elev/2) ; wrapTo2Pi(bin_centers_elev+bin_width_elev/2)]);

out={bin_centers_azm,bin_pairs_azm,bin_centers_elev,bin_pairs_elev};

end