function counts_chunked=chunk_data_complete(counts,sample_proportion,sort_dir,num_chunks)
%data chunking function for normalization
%see calc_any_g2_type for usage

%input
%norm_samp_factor-should be about the correlation amp for equal noise contibution from in-shot & between shot
%data.txy
%data.data.num_counts
%to do
%   -documentation

all_counts = cat(1, counts{:});

num_counts=cellfun(@(x)size(x,1),counts);

num_shots = size(counts,2);

total_counts=sum(num_counts);

sum_counts=cumsum(num_counts);

if sample_proportion>1 || sample_proportion<=0
   error('Sample proportion must be between 0 and 1') 
end

if nargin<4
    num_chunks = total_counts-num_counts(1);
end

if num_chunks>total_counts-num_counts(1)
    num_chunks = total_counts-num_counts(1);
end

counts_chunked=cell(1,num_chunks);
% for ii=1:num_shots
%     this_counts = counts{ii};
%     chunk_dist = randperm(num_chunks);
%     for jj = 1:num_counts(ii)
%         box = chunk_dist(jj);
%         counts_chunked{box}=[counts_chunked{box};this_counts(jj,:)];
%     end
% end
% for ii=1:num_shots
%     this_counts = counts{ii};
%     chunk_dist = randperm(num_chunks);
%     for jj = 1:num_chunks%num_counts(ii)
%         box = chunk_dist(jj);
%         counts_chunked{box}=[counts_chunked{box};this_counts(mod(jj,num_counts(ii))+1,:)];
%     end
% end
total_indx = 1;
for ii = 2:num_shots%sum_counts(1)+1:total_counts
    for jj = 1:num_counts(ii)
    counts_chunked{total_indx}=[counts_chunked{total_indx};counts{ii}(jj,:)];
    indxs = 1:(total_indx+sum_counts(1)-1);
    indxs = setdiff(indxs,sum_counts(ii-1)+1:sum_counts(ii));
    sampel_factor = round(length(indxs)*sample_proportion);
    sampled_indxs = randsample(indxs,sampel_factor);
    counts_chunked{total_indx}=[counts_chunked{total_indx};all_counts(sampled_indxs,:)];
    total_indx =  total_indx + 1;
    if total_indx>num_chunks
        break
    end
    end
    if total_indx>num_chunks
        break
    end
end

if ~isnan(sort_dir)
    for ii = 1:num_chunks
        tmp_data = counts_chunked{ii};
        [~,order]=sort(tmp_data(:,sort_dir));
        tmp_data=tmp_data(order,:);
        counts_chunked{ii} = tmp_data;
    end
end
% if size(vertcat(counts_chunked{:}),1)~=total_counts
%     warning('lost counts')
% end
end