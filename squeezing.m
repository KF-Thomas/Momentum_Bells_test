function out=squeezing(halo_centered_cells,isverbose)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
% new file 1/2/23 - M.W

%calulates squezing
%inputs 
    %cell array for eahc file containing Z,X,Y counts
    %isverbose = 1 displays some info, =0 does not. 
%outputs
    %a cell that contains two matrices, the first contains what is needed
    %for the sqz angle plot (angle,variance,variance_se)
    %the second which contains [mean_opst_bin,unc_opst_bin,mean_other_bin,unc_other_bin]
%this takes the centered halo data as cells for each shot of z,x,y 
%converts it to radial cord and then finds the normalized number vairance
%the traditional way (http://arxiv.org/pdf/1008.0845.pdf)of 
%showing this squezing data is to plot the disjointed pairs and the jointed
%pairs, this seems crude and arbitrary( the x axis has no information)
%instead i plot using the angle between the bin centers when this is pi
%then you expect decreased num var above it not


% ref system
% BEC's are at the poles
% azimuthal angle is arround the equator (zero to 2 pi)
% inclination is angle from the poles (which is zero inclination)(zero to pi)
% this was chosen to reduce the wrapping problems arround zero/2 pi 


%TO FIX
%the double loop just after disp('Fiding Norm Var (rearanging)') takes the
%majority of the time, by doing some smarter matrix stuff it may be sped up

%bins that wrap across 0,2pi lose counts


%================================ START USER INPUT=======================================

plot_sqz_bins=1;%plot sqz with bin #
plot_sqz_angle=1;%plot sqz with angle


% define bins
mirror_azm=0; %use for the mag sensitive halo which you must cut out a region ??

range_azm=[0 2]*pi;%azimuthal range for bins. to look at whole halo, [0 2]*pi
range_elev= [0 1]*pi; %[0.0 1]*pi;%Elevation range. whole halo [0.0 1]*pi

steps_azm=4; %number of bins in azimuthal direction
steps_elev=2; %number of bins in elev direction

%defining the width seprately allows for under or over sampling bins
bin_width_azm=2*range(range_azm)/steps_azm; 
bin_width_elev=range(range_elev)/steps_elev;

%binning notes to work out:
%if you make steps_azm small and steps_elev big something gets upset later
%bin width has a significant effect on results
%note from to fix: apparently bins wrapping across 0,2pi lose counts


%===================================END USER INPUT========================================

%remove empty bins
halo_centered_cells=halo_centered_cells(~cellfun('isempty',halo_centered_cells));

% create the bins
%bit hacky to create bins that dont wrap across 0,2pi
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


% mask=mask_rad & halo_azm>window_azm(1) & halo_azm<window_azm(2);
% num_origin=sum(mask);

if isverbose
disp('binning') %indicate that binning is complete
end

% convert halo to spherical coordinates and count atoms in each bin
% halo_counts=[]; %not used
% halo_centered_spherical={}; %not used

angle_counts={}; %stores number of counts in each bin for each shot. Each element is an array containing the shot number, azm/elev bin centres and counts

for n=1:size(halo_centered_cells,2)
    angle_counts_file={};
    single_halo=halo_centered_cells{n};
    
    [halo_radial,halo_azm,halo_elev]=ConvToSph(single_halo);
    
    %halo_centered_spherical{n} =[halo_radial,halo_azm,halo_elev];
    %[halo_radial,halo_azm/pi,halo_elev/pi]
    %mask_rad=halo_radial>radial_min & halo_radial<radial_max;
    %halo_counts(n)=sum(mask_rad);
    
    %launch into the nested loop of bining
    for m=[1:size(bin_pairs_azm,1)]
        for p=[1:size(bin_pairs_elev,1)]
            %count the number of atoms in bin m,p
        counts=sum(halo_azm>bin_pairs_azm(m,1) & halo_azm<bin_pairs_azm(m,2)  & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2)); %mask_rad
        %store this number in 
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

%exclude combinations if the two bins overlap (these will have artificially
%low variance)
if bin_width_azm > range(range_azm)/steps_azm || bin_width_elev > range(range_elev)/steps_elev 
    if isverbose
        disp(['deleting overlapping bins'])
    end
    index_combs = delete_overlapping_bins(index_combs, steps_azm, steps_elev, bin_width_azm, bin_width_elev, range_azm, range_elev);
    index_combs
end

% number of halos , number of bin combinations, 2
bin_pairs=zeros(size(halo_centered_cells,1),size(index_combs,1),2);

% old way of rearanging the bins
% disp('Fiding Norm Var (rearanging)')
% tic
% for n=[1:size(halo_centered,1)]
%     for m=[1:size(index_combs,1)]
%         %finds the counds in the bin defined by the index_combs
%         temp_angle_counts=angle_counts{n};
%         bin_pairs(n,m,:)=[temp_angle_counts(index_combs(m,1),4),temp_angle_counts(index_combs(m,2),4)];
%     end
% end
% clear temp_angle_counts;
% toc
% size(bin_pairs)

%first we rearange the cell output into an array
%there may be a better way to do this other than a loop but i cant find it
%predefine array ? is the cell actually used for anything or can we just
%put stuff straight into angle_counts_mat as they are calculated ?
if isverbose
disp('Finding Norm Var ')
end

% angle_counts_mat=[];
% for n=1:size(angle_counts,2) %size(angle-counts,2) is the number of shots
%     angle_counts_mat(:,:,n)=angle_counts{n};
% end

%this should be faster
num_shots = size(angle_counts,2); %check what happens if you have multiple halos
angle_counts_mat=zeros(num_bins,size(angle_counts{1},2),num_shots);
for n=1:size(angle_counts,2) %size(angle-counts,2) is the number of shots
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


% %the above is little more that rearangeing values, below is where the norm
% %var is actualy calculated
% diffs=(bin_pairs(:,:,1)-bin_pairs(:,:,2));
% diffs2=(diffs).^2;
% norm_var=(mean(diffs2,1)-mean(diffs,1).^2)./(mean(bin_pairs(:,:,1))+mean(bin_pairs(:,:,2)));

%the above is little more that rearangeing values, below is where the norm
%var is actualy calculated
diffs=(bin_pairs(:,:,1)-bin_pairs(:,:,2));
diffs2=(diffs).^2;
norm_var=(mean(diffs2,1,'omitnan')-mean(diffs,1, 'omitnan').^2)./(mean(bin_pairs(:,:,1),'omitnan')+mean(bin_pairs(:,:,2), 'omitnan'));



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
    angle_pairs(n)= acos(sin(elev1)*sin(elev2)*cos(azm1-azm2)+cos(elev1)*cos(elev2));
end

%need to exclude pairs with adjacent bins
%most efficient way to do this in general? Check if bins overlap, if so don't include them in index_combs 

%now i want to sort through the angle and norm var values, average the
%norm_var values that have the same angle_pairs
angle_tol= 0.0001;% 0.0001;
uniq_angles=uniquetol(angle_pairs,angle_tol);
angle_var_sd=zeros(size(uniq_angles,2),3);
%here i also remove zeros
for n=1:size(uniq_angles,2)
    matching=norm_var((abs(angle_pairs-uniq_angles(n))<angle_tol) & norm_var~=0);
    angle_var_sd(n,:)=[uniq_angles(n)/pi, mean(matching),std(matching)/sqrt(size(matching,2))]; %omitnan
end

if plot_sqz_angle
    
    figure(10)
    errorbar(angle_var_sd(:,1),angle_var_sd(:,2),angle_var_sd(:,3),'x')
    %title(['Squeezing (',num2str(steps_azm),' azm ',num2str(steps_elev),' elev )']);
    set(gcf,'Color',[1 1 1]);
    xlabel('Angle Between Bins /Pi')
    ylabel('Normalised Variance')
    line([-2 2], [1 1],'Color','red');
    set(gca,'Xlim',[-0.1 1.1])

end

%also plot the more traditional type of plot where the x axis is the bin
%comb number, two types of points one for the opst and the other for all
%others

norm_var_opst=norm_var((abs(angle_pairs-pi)<angle_tol) & norm_var~=0);
norm_var_rest=norm_var(~(abs(angle_pairs-pi)<angle_tol) & norm_var~=0);

norm_var_avg=angle_var_sd(:,1);

mean_opst_bin=angle_var_sd(angle_var_sd(:,1)==1,2);
unc_opst_bin=angle_var_sd(angle_var_sd(:,1)==1,3);
non_opst_bins=angle_var_sd(angle_var_sd(:,1)~=1,2);
mean_other_bin=mean(non_opst_bins); %omitnan
min_other_bin=min(non_opst_bins);
unc_other_bin=std(non_opst_bins)/sqrt(length(mean_other_bin));
min_opst_bin=min(norm_var(angle_pairs==pi));

if  isverbose
    disp(['mean opst bin' ,num2str(mean_opst_bin),'±',num2str(unc_opst_bin)])
    disp(['mean other bins',num2str(mean_other_bin),'±',num2str(unc_other_bin)])
    disp(['min opst bin ',num2str(min_opst_bin)])
end
if plot_sqz_bins
    figure(11)
    plot(1:size(norm_var_opst,2),norm_var_opst,'x',1:size(norm_var_rest,2),norm_var_rest,'+')
    xlabel('Bin Pair Number')
    ylabel('Normalised Variance')
    line([0 size(norm_var_rest,2)], [1 1],'Color','red');
end

out={[angle_var_sd(:,1),angle_var_sd(:,2),angle_var_sd(:,3)],[[mean_opst_bin,unc_opst_bin,min_opst_bin];[mean_other_bin,unc_other_bin,min_other_bin]]};

end