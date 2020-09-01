function counts_chunked=chunk_data_complete(counts,sample_proportion,sort_dir)
%data chunking function for normalization
%see calc_any_g2_type for usage

%input
%norm_samp_factor-should be about the correlation amp for equal noise contibution from in-shot & between shot
%data.txy
%data.data.num_counts
%to do
%   -documentation

% all_counts = cat(1, counts{:});

% num_counts=cellfun(@(x)size(x,1),counts);

num_shots = size(counts,2);

% total_counts=sum(num_counts);

% sum_counts=cumsum(num_counts);

if sample_proportion>1 || sample_proportion<=0
    error('Sample proportion must be between 0 and 1')
end
num_chunks= round(num_shots^2.*sample_proportion);
sample_mask = [zeros(1,num_shots^2-num_chunks) ones(1,num_chunks)];  %# Fill the vector with 0 and 1
sample_mask = sample_mask(randperm(num_shots^2));  %# Randomly reorder it
counts_chunked=cell(2,num_chunks);
total_indx = 1;
chunk_indx = 1;

if size(counts,1)>1
    for ii = 1:num_shots%sum_counts(1)+1:total_counts
        for jj = 1:num_shots
            if sample_mask(total_indx)
                counts_chunked{1,chunk_indx} = counts{1,ii};
                counts_chunked{2,chunk_indx} = counts{2,jj};
                chunk_indx = chunk_indx + 1;
            end
            total_indx =  total_indx + 1;
        end
    end
else
    for ii = 1:num_shots%sum_counts(1)+1:total_counts
        for jj = 1:num_shots
            if sample_mask(total_indx)
                counts_chunked{1,chunk_indx} = counts{ii};
                counts_chunked{2,chunk_indx} = counts{jj};
                chunk_indx = chunk_indx + 1;
            end
            total_indx =  total_indx + 1;
        end
    end
end

num_chunks = size(counts_chunked,2);
if ~isnan(sort_dir)
    for ii = 1:num_chunks
        for jj = 1:2
            tmp_data = counts_chunked{jj,ii};
            [~,order]=sort(tmp_data(:,sort_dir));
            tmp_data=tmp_data(order,:);
            counts_chunked{jj,ii} = tmp_data;
        end
    end
end

end