function index_combs = delete_overlapping_bins(index_combs, steps_azm, steps_elev, bin_width_azm, bin_width_elev, range_azm, range_elev,isverbose)
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%M.W.
%this function removes combinations containing two overlapping bins
%
%it uses the bin_width to determine is adjacent bins overlap, 
%then removes overlapping combinations from index_combs
%
%It works if each bin only overlaps with bins adjacent to it.
%If the bin widths are large enough so that it overlaps with more bins,
%it displays a warning. choose better bins. 
%
%Bin indexing
%index_combs indexes the bins by number, as azm increases then elev (from
%the for loop that produces angle_counts in the squeezing.m code)
% if we have m azm sections and p elev sections, then the nth bin is
% the one with the (n mod m) azm section and the ceiling(n/m) elev section
%
%bins 1 to m have the lowest elevation, bins m+1 to 2m have the next lowest etc...


%number of overlapping bins on either side of a given bin
num_overlapping_azm = ceil(bin_width_azm/(range(range_azm)/steps_azm)) - 1;
num_overlapping_elev = ceil(bin_width_elev/(range(range_elev)/steps_elev)) - 1;

if num_overlapping_azm>1 || num_overlapping_elev>1
    disp(['warning: bins overlap too much']) 
end

%delete overlapping bins
overlapping_bins=[];

if num_overlapping_azm>0 || num_overlapping_elev>0 %find bins overlapping in azm direction
    for n = 1:size(index_combs,1)
        %hacky way of determining if two bins are overlapping: 
        if (mod(index_combs(n,2),steps_azm)==mod(index_combs(n,1),steps_azm) || mod(abs( mod(index_combs(n,2),steps_azm)-mod(index_combs(n,1),steps_azm) ),steps_azm-2) == num_overlapping_azm )... %check if two bins are within num_overlapping_azm of each other
                && (ceil(index_combs(n,2)/steps_azm)==ceil(index_combs(n,1)/steps_azm) || abs( ceil(index_combs(n,2)/steps_azm)-ceil(index_combs(n,1)/steps_azm) ) == num_overlapping_elev) %check if two bins are within num_overlapping_elev of each other
            %disp([index_combs(n,1) index_combs(n,2) (mod(index_combs(n,2),steps_azm)==mod(index_combs(n,1),steps_azm) || mod(abs( mod(index_combs(n,2),steps_azm)-mod(index_combs(n,1),steps_azm) ),steps_azm-2) == num_overlapping_azm ) (ceil(index_combs(n,2)/steps_azm)==ceil(index_combs(n,1)/steps_azm) || abs( ceil(index_combs(n,2)/steps_azm)-ceil(index_combs(n,1)/steps_azm) ) == num_overlapping_elev)])
            overlapping_bins(end+1)= n;
        end
    end
end

index_combs = removerows(index_combs,'ind',overlapping_bins); %check that this works
if isverbose
    disp(['number of combinations of bins (overlaps removed)',int2str(size(index_combs,1))]);
end


% %slow method that checks if bins overlap using coordinates
% %requires more inputs from squeezing.m (angle_counts)
% bin_list = angle_counts{1}(:,[2,3]); %list of centre coordinates of all bins
% if bin_width_azm > range(range_azm)/steps_azm || bin_width_elev > range(range_elev)/steps_elev %first check if bins widths are big enough to overlap
%     if isverbose
%         disp('accounting for overlapping bins')
%     end
%     index_combs_coords = zeros(size(index_combs,1),6);
%     for n = 1:size(index_combs,1)
%         index_combs_coords(n,:)= [index_combs(n,:) bin_list(index_combs(n,1),:) bin_list(index_combs(n,2),:)];
%     end
%     %go through and check if they overlap
%     overlapping_bins = [];
%     for n = 1:size(index_combs,1)
%         i_azm_center = index_combs_coords(n,3);
%         j_azm_center = index_combs_coords(n,5);
%         i_elev_center = index_combs_coords(n,4);
%         j_elev_center = index_combs_coords(n,6);
% 
%         i_azm_range = fixed.Interval(i_azm_center-bin_width_azm/2, i_azm_center+bin_width_azm/2,'()');
%         j_azm_range = fixed.Interval(j_azm_center-bin_width_azm/2, j_azm_center+bin_width_azm/2,'()');
%         i_elev_range = fixed.Interval(i_elev_center-bin_width_elev/2, i_elev_center+bin_width_elev/2,'()');
%         j_elev_range = fixed.Interval(j_elev_center-bin_width_elev/2, j_elev_center+bin_width_elev/2,'()');
%         %open intervals so that they are only flagged as overlapped if
%         %they share more than the endpoints in common
% 
% 
%         if overlaps(i_azm_range, j_azm_range) && overlaps(i_elev_range, j_elev_range)
%             % if so, delete the nth entry of index_combs_coords
%             overlapping_bins(end+1) = n;
%             % disp(['overlapping bins: ',num2str(index_combs(n,:))]) %check which bins are flagged as overlapping
%         end    
%     end
%     %delete overlapping bins
%     index_combs = removerows(index_combs,'ind',overlapping_bins); %check that this works
%     if isverbose
%         disp(['number of combinations of bins (overlaps removed)',int2str(size(index_combs,1))]);
%     end
% end



end