%testing out a correlation idea
counts = out_data{6}.bottom_halo.counts_txy;
top_bound = [1.762,1.7655];%[1.7603,1.767];
btm_bound = [1.7485,1.753];%[1.7476,1.754];
one_d_edges = linspace(0,0.3,300);
one_d_cens = (one_d_edges(1:end-1)+one_d_edges(2:end))./2;
one_d_bins = zeros(size(one_d_edges,2)-1,1);

for ii=1:length(counts);
this_txy = counts{ii};

time_mask = (this_txy(:,1)>top_bound(1) & this_txy(:,1)<top_bound(2)) |...
    (this_txy(:,1)>btm_bound(1) & this_txy(:,1)<btm_bound(2));
masked_counts = this_txy(time_mask,:);

masked_counts(:,1) = masked_counts(:,1).*9.81.*0.417;
delta = pdist(masked_counts);

one_d_bins=one_d_bins+hist_adaptive_method(delta.',one_d_edges.');
end

