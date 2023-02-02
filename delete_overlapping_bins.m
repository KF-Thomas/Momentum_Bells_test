function index_combs = delete_overlapping_bins(index_combs, steps_azm, steps_elev, bin_width_azm, bin_width_elev, range_azm, range_elev)
%exclude combinations if the two bins overlap (these will have artificially
%low variance)
%index_combs indexes the bins by number, as azm increases then elev (from
%the for loop that produces angle_counts)
% if we have m azimuthal sections and p elev sections, then the nth bin is
% the one with the (n mod m) azm section and the ceiling(n/m) elev section

%if there are m=steps_azm bins in azm direction, then bin
                %numbering swaps at multiples of m. So bins 1:m are adjacent in azm direction, m and m+1 are not

                
%use bin_width to tell how many adjacent bins overlap, then remove these
%pairs from index_combs

%number of overlapping bins on either side of a given bin
num_overlapping_azm = ceil(bin_width_azm/(range(range_azm)/steps_azm)) - 1
num_overlapping_elev = ceil(bin_width_elev/(range(range_elev)/steps_elev)) - 1

%delete overlapping bins
%need to adapt to take into account num_overlapping
overlapping_bins=[];
if num_overlapping_azm>0 %find bins overlapping in azm direction
    disp(['finding overlapping azm'])
    for n = 1:size(index_combs,1)
        if ( abs( mod(index_combs(n,2)-index_combs(n,1),steps_azm) ) <= num_overlapping_azm )... %check if two bins are within num_overlappin_azm of each other
                && abs( ceil(index_combs(n,2)/steps_azm)-ceil(index_combs(n,1)/steps_azm) ) <= num_overlapping_elev %check if two bins are within num_overlappin_elev of each other
            overlapping_bins(end+1)= n;
        end
    end
end

index_combs = removerows(index_combs,'ind',overlapping_bins); %check that this works
%disp(index_combs)
disp(['number of combinations of bins (overlaps removed)',int2str(size(index_combs,1))]);


% %method with lots of holes
% if bin_width_azm > range(range_azm)/steps_azm || bin_width_elev > range(range_elev)/steps_elev %first check if bins widths are big enough to overlap
%     if isverbose
%         disp('accounting for overlapping bins')
%     end
%     
%     %number of overlapping bins on either side of a given bin
%     num_overlapping_azm = ceil(bin_width_azm/(range(range_azm)/steps_azm)) - 1
%     num_overlapping_elev = ceil(bin_width_elev/(range(range_elev)/steps_elev)) - 1
%     
%     %delete overlapping bins
%     %need to adapt to take into account num_overlapping
%     overlapping_bins=[];
%     if num_overlapping_azm>0 %find bins overlapping in azm direction
%         disp(['finding overlapping azm'])
%         for n = 1:size(index_combs,1)
%             if index_combs(n,2)-index_combs(n,1)<=num_overlapping_azm && mod(index_combs(n,1),steps_azm)~=0 ...
%                     || index_combs(n,2)-index_combs(n,1)==steps_azm-1 && mod(index_combs(n,1),steps_azm)==1 && mod(index_combs(n,2),steps_azm)==0
% 
%                 %if there are m=steps_azm bins in azm direction, then bin
%                 %numbering swaps at multiples of m. So bins 1:m are adjacent in azm direction, m and m+1 are not
%                 overlapping_bins(end+1)= n;
%             end
%         end
%     end
%     if num_overlapping_elev>0 %find bins overlapping in elev direction
%         disp(['finding overlapping elev'])
%         for n = 1:size(index_combs,1)
%             if index_combs(n,2)-index_combs(n,1)==steps_azm %make this more sophisticated to account for m and p
%                 %if there are m=steps_azm bins in azm direction, then bin
%                 %numbering swaps at multiples of m. So bins 1:m are adjacent in azm direction, m and m+1 are not
%                 overlapping_bins(end+1)= n;
%             end
%         end
%     end
%     if num_overlapping_elev>0 && num_overlapping_azm>0 %find bins overlapping in both directions
%         disp(['finding overlapping elev and azm'])
%         for n = 1:size(index_combs,1)
%             if index_combs(n,2)-index_combs(n,1)==steps_azm+num_overlapping_azm... %these all have holes
%                     ||index_combs(n,2)-index_combs(n,1)==steps_azm-num_overlapping_azm...
%                     ||index_combs(n,2)-index_combs(n,1)==steps_azm-1 && mod(index_combs(n,1),steps_azm)==1 && mod(index_combs(n,2),steps_azm)==0
%  %make this more sophisticated to account for m and p
%                 %if there are m=steps_azm bins in azm direction, then bin
%                 %numbering swaps at multiples of m. So bins 1:m are adjacent in azm direction, m and m+1 are not
%                 overlapping_bins(end+1)= n;
%             end
%         end
%     end
%     index_combs = removerows(index_combs,'ind',overlapping_bins); %check that this works
%     disp(index_combs)
%     if isverbose
%         disp(['number of combinations of bins (overlaps removed)',int2str(size(index_combs,1))]);
%     end
% end

% %slow method
% %don't need all 6 in index_combs_coords
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