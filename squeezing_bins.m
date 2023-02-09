function out=squeezing_bins(halo_centered_cells,bin_info,isverbose)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%calulates squeezing, with bin inputs
%inputs 
    %cell array containing Z,X,Y counts for centered halo
    %bin_info ; struct containing specifications needed for binning
    %isverbose = 1 displays some info, =0 does not. 
%outputs
    %bin_pairs ; array containing counts for each bin combination for each shot
    %total_counts ; array containing total counts in each shot
    %angle_pairs ; array containing angle between each bin pair
%--------------------------------------------------------------------------

% define bins
range_azm=bin_info.range_azm;
range_elev=bin_info.range_elev;
steps_azm = bin_info.steps_azm;
steps_elev= bin_info.steps_elev;
bin_width_azm = bin_info.bin_width_azm;
bin_width_elev= bin_info.bin_width_elev;
cuts_azm = bin_info.cuts_azm;
cuts_elev= bin_info.cuts_elev;
mirror_azm=bin_info.mirror_azm; %use for the mag sensitive halo which you must cut out a region ??

%remove empty bins
halo_centered_cells=halo_centered_cells(~cellfun('isempty',halo_centered_cells));

% create the bins
bins = create_bins(bin_info);
bin_centers_azm = bins{1};
bin_pairs_azm = bins{2};
bin_centers_elev = bins{3};
bin_pairs_elev = bins{4}; 

if isverbose
disp('binning')
end

angle_counts={}; %stores number of counts in each bin for each shot. Each element is an array containing the shot number, azm/elev bin centres and counts
shot_num=size(halo_centered_cells,2);
total_counts=zeros(shot_num,1);

for n=1:shot_num
    angle_counts_file={};
    single_halo=halo_centered_cells{n};
    
    total_counts(n)=size(single_halo,1);
    
    [halo_radial,halo_azm,halo_elev]=ConvToSph(single_halo);

    %launch into the nested loop of binning
    for m=[1:size(bin_pairs_azm,1)]
        for p=[1:size(bin_pairs_elev,1)]
            %count the number of atoms in bin m,p
        
            %separate clause for bins wrapping around 0,2pi to prevent lost counts
            if bin_pairs_azm(m,1)>bin_pairs_azm(m,2)                
                %break bin into two parts (x1,2pi) and (0,x2) and count separately
                counts = sum( halo_azm>bin_pairs_azm(m,1) & halo_azm<2*pi ...
                    & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2))...
                    + sum( halo_azm>0 & halo_azm<bin_pairs_azm(m,2) ...
                    & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2));           
            else 
            
            counts=sum(halo_azm>bin_pairs_azm(m,1) & halo_azm<bin_pairs_azm(m,2)...
                & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2)); %mask_rad
            end
                   
            angle_counts_file{m,p}=[n,bin_centers_azm(m),bin_centers_elev(p),counts];        
        end
    end
    angle_counts{n}=vertcat(angle_counts_file{:});
end


%now we have the bins we can find the normalized number vairance
%between them as defined in http://arxiv.org/pdf/1008.0845.pdf

num_bins = size(angle_counts{1},1);
%all possible combinations of two bins is given by
index_combs=nchoosek(1:num_bins,2);
if isverbose
    disp(['number of combinations of bins ',int2str(size(index_combs,1))]);
end

%exclude combinations if bins overlap (these will have artificially low variance)
exclude_overlapping_bins = 1;
if exclude_overlapping_bins
    if bin_width_azm > range(range_azm)/steps_azm || bin_width_elev > range(range_elev)/steps_elev 
        if isverbose
            disp(['deleting overlapping bins'])
        end
        index_combs = delete_overlapping_bins(index_combs, steps_azm, steps_elev, bin_width_azm, bin_width_elev, range_azm, range_elev,isverbose);
    end
end

%exclude combinations if any bin overlaps with removed parts of halo
%find bins from bin_pair_azm that overlap with cuts, work out what number
%this corresponds to and get rid of those
%only works for azm cuts (doesn't distinguish elevation)
exclude_cut_bins=1;
if exclude_cut_bins
    cut_bins_azm=[];
    for m=1:size(bin_pairs_azm,1)
        if bin_pairs_azm(m,1)>bin_pairs_azm(m,2)
            interval1=fixed.Interval(bin_pairs_azm(m,1),2*pi,'()');
            interval2=fixed.Interval(0,bin_pairs_azm(m,2),'()');
            bin_range=union(interval1,interval2);
        else
        bin_range=fixed.Interval(bin_pairs_azm(m,1),bin_pairs_azm(m,2),'()'); %work out what to do with 2pi
        end
        %check if it overlaps with any cuts (including opposite to a cut!)
        for i=1:size(cuts_azm,1)
            cut_range=fixed.Interval(cuts_azm(i,1),cuts_azm(i,2),'()');
            cut_range_opp = fixed.Interval(wrapTo2Pi(cuts_azm(i,1)+pi),wrapTo2Pi(cuts_azm(i,2)+pi),'()');
            if any(overlaps(bin_range,cut_range)) || any(overlaps(bin_range,cut_range_opp))
                cut_bins_azm = [cut_bins_azm m];

                %if it does overlap with a cut, delete bin combinations involving this bin
                if m==steps_azm
                    cut_mask = any(mod(index_combs(:,:),steps_azm)==0,2) ;
                    index_combs = index_combs(~cut_mask,:);
                    %disp(num2str(size(index_combs,1)))
                else
                    cut_mask = any(mod(index_combs(:,:),steps_azm)==m,2) ;
                    index_combs = index_combs(~cut_mask,:);
                    %disp(num2str(size(index_combs,1)))
                end

            end
        end

    end 
end

if isverbose
    disp(['number of combinations of bins (cuts removed) ',num2str(size(index_combs,1))])
end

if size(index_combs,1)==0
    disp(['Error: 0 valid bin combinations. Increase bin #, reduce overlaps or reduce cuts'])
end

% number of halos , number of bin combinations, 2
bin_pairs=zeros(size(halo_centered_cells,1),size(index_combs,1),2);


%first we rearange the cell output into an array
%there may be a better way to do this other than a loop but i cant find it
%predefine array ? is the cell actually used for anything or can we just
%put stuff straight into angle_counts_mat as they are calculated ?
if isverbose
disp('Finding Norm Var ')
end

num_shots = size(angle_counts,2); 
angle_counts_mat=zeros(num_bins,size(angle_counts{1},2),num_shots);
for n=1:num_shots %size(angle-counts,2) is the number of shots
    angle_counts_mat(:,:,n)=angle_counts{n};
end
%angle_counts_mat [bin number, angle_counts_info, shot]
%to get counts in bin x, shot n, angle_counts_mat(x,4,n)

%this selects the appropriate bins so we have a matix of filex*index in
%comb list*(bin1,bin2)
bin_pairs=[angle_counts_mat(index_combs(:,1),4,:), angle_counts_mat(index_combs(:,2),4,:)]; %[counts in bin1, counts in bin2]
%permute these get desired format: first index = shot, second = bin comb.,
%third index = which of the two bins in the combination
bin_pairs=permute(bin_pairs,[3 1 2]);

%want to calculate the angle between the bin_centers 
%https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
%https://math.stackexchange.com/questions/243142/what-is-the-general-formula-for-calculating-dot-and-cross-products-in-spherical
%r1r2(sin?1sin?2cos(?1??2)+cos?1cos?2)
%acos(sin(elev1)sin(elev2)cos(azm1?azm2)+cos(elev1)cos(elev2))
angle_pairs=zeros(1,size(index_combs,1));
for n=[1:size(index_combs,1)]
    azm1=angle_counts{1}(index_combs(n,1),2);
    elev1=angle_counts{1}(index_combs(n,1),3);
    azm2=angle_counts{1}(index_combs(n,2),2);
    elev2=angle_counts{1}(index_combs(n,2),3);
    
    angle_pairs(n)= real(acos(sin(elev1)*sin(elev2)*cos(azm1-azm2)+cos(elev1)*cos(elev2)));
    %only take real part to avoid error where acos returns a tiny imaginary component
end

out={bin_pairs,total_counts,angle_pairs};

end