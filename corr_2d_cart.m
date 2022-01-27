function out=corr_2d_cart(corr_opts,counts)
% corr_1d_cart - %a low/high memory implementation of the one D correlation
% Notes:
% This g2 correlator is very fast and can be thought as taking a tube along the full 3d cartesian g2 but without having
% to deal with how slow that is. It does this by windowing the pair differences in two dimensions then histograming the
% differences in the remaining dimension(one_d_dimension). Because the correlation is effectively integerated in those
% two window dimensions the window width (one_d_window) must be a small fraction of the correlation length (in those dimensions)
%to avoid suppressing the measured corelation amplitude.
%
% Usage:
% pass in a cell array of counts and then calculate the correlations using a
% bin on the fly method (low mem) or a bin all diferences (high mem)
% differences can be prewindowed using fast_sorted_mask (optional) in both cases to reduce memory & comp requirements.
% differences are then windowed in two dimensions and the correlation calculated in the remaining axis.
% conforms to output formats of import_mcp_tdc_data
%
% Syntax:  out=corr_1d_cart(corr_opts,counts_txy)
%
% Inputs:
%	counts_txy                  - as a cell array of zxy counts size=[num_counts,3]
%               must be ordered in z for pre window optimzation
%   corr_opts                   - structure of input options
%           .one_d_window       - [[tmin,tmax];[xmin,xmax];[ymin,ymax]]; %limiits in the one_d_dimension are ignored
%           .cl_or_bb           - colinear(false) or bb (true)
%           .one_d_edges        - edges of histogram bins for the correlation
%           .one_d_dimension    - dimension to calculate correlation
%           .one_d_smoothing    - [float] gaussian smoothing for the  coincidence density output use nan or 0 for no
%                                 smoothing
%           .attenuate_counts   - value 0-1 the ammount of data to keep, like simulating a lower QE
%           .do_pre_mask        - logical, use fast sorted mask to window to a range of
%                                 one_d_window in the sorted axis arround each count
%               .sorted_dir     - must be passed when corr_opts.do_pre_mask used, usualy 1 for time
%           .low_mem            - use lowmem(true) or highmem(false) options, use nan to set dynamicaly using the
%                                 remaining memory on your computer.
%
% Outputs:
%   out                         - output structure
%       .one_d_corr_density     - (optional smoothed) coincidence density
%       .one_d_corr_density_raw - raw coincidence density
%       .x_centers              - centers of the histogram vector
%       .pairs                  - total number of pairs
%
% Example:
%    [TO DO]
%
% Other m-files required: fast_sorted_mask
% Subfunctions: none
% MAT-files required: none
% See also: corr_unit_testing, calc_any_g2_type, import_mcp_tdc_data
%
% Known BUGS/ Possible Improvements
%   - lowmem option seems to be much faster at the moment (with less cpu usage)
%   - see if better to prevent prevent diff vector initalization when using high meme with prewindowing to save memory
%   - auto memory option does not work on linux, should defalult to low_mem
%   - input checking
%	- document output better
%	- more outputs
%   - for the high memory option I dont have to store the prewindowed dimension
%     - might be a little messy to implement
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-12-18

%------------- BEGIN CODE --------------
if size(counts,1)==2
    corr_opts.between_sets = true;
else
    corr_opts.between_sets = false;
end


if ~isfield(corr_opts,'attenuate_counts')
    corr_opts.attenuate_counts=1;
else
    if corr_opts.attenuate_counts>1 || corr_opts.attenuate_counts<0
        error('invalid attenuation value, pass a value between 0 and 1\n')
    end
end

if ~isfield(corr_opts,'cl_or_bb')
    error('corr_opts.cl_or_bb not specified!do you want co-linear(false) or back-to-back (true)? ')
end

if isfield(corr_opts,'do_pre_mask') && ~isfield(corr_opts,'sorted_dir')
    error('you must pass corr_opts.sorted_dir with corr_opts.do_pre_mask')
elseif ~isfield(corr_opts,'do_pre_mask')
    corr_opts.do_pre_mask=false;
end

if ~isfield(corr_opts,'progress_updates')
    corr_opts.progress_updates=50;
end

num_counts=cellfun(@(x)size(x,1),counts);

% convert edges to col vec
% corr_opts.two_d_edges=col_vec(corr_opts.two_d_edges);

%for the dimension of interest we set the window to be the region specified in corr_opts.one_d_edges
%this simplifies the code and can speed up the historgraming
corr_opts.one_d_window(corr_opts.two_d_dimensions,:)=[min(corr_opts.two_d_edges{1}),max(corr_opts.two_d_edges{1});...
    min(corr_opts.two_d_edges{2}),max(corr_opts.two_d_edges{2})];

updates=corr_opts.progress_updates; %number of updates in the progress bar to give the user, can slow thigs down if too high


shots =size(counts,2);
update_interval=ceil(shots/updates);
N = ceil(shots/update_interval);
% dvs = divisors(shots);
% [~,idx]=min(abs(dvs-N));
% N=dvs(idx);
if ~isfield(corr_opts,'print_update')
    corr_opts.print_update = true;
end
if corr_opts.print_update
    parfor_progress_imp(N);
end

pairs_count=zeros(1,shots);
delta_multiplier=[1,-1];
delta_multiplier=delta_multiplier(1+corr_opts.cl_or_bb); %gives 1 when cl and -1 when bb

%build a queue of the dimensions that should be masked in before the histogram
%should not mask in corr_opts.one_d_dimension & corr_opts.sorted_dir (if corr_opts.do_pre_mask )
mask_dimensions=[1,2,3];
mask_dimensions_logic=[true,true,true]; %logic to avoid clashes when corr_opts.one_d_dimension==corr_opts.sorted_dir
mask_dimensions_logic(corr_opts.two_d_dimensions)=false;
if corr_opts.do_pre_mask && ~corr_opts.between_sets, mask_dimensions_logic(corr_opts.sorted_dir)=false; end
mask_dimensions=mask_dimensions(mask_dimensions_logic);

[two_d_bins,~]=histcn([Inf,Inf;Inf,Inf],corr_opts.two_d_edges{1},corr_opts.two_d_edges{2});

    % the low memory mode is serial and is a bit easier on mem requirements
%     two_d_bins=col_vec(two_d_bins);
    for shotnum=1:shots
        if corr_opts.between_sets
            shot_txy_1=counts{1,shotnum};
            shot_txy_2=counts{2,shotnum};
            num_counts_shot_1=num_counts(1,shotnum);
            num_counts_shot_2=num_counts(2,shotnum);
        else
            shot_txy=counts{shotnum};
            num_counts_shot=num_counts(shotnum);
        end
        if corr_opts.attenuate_counts~=1 %randomly keep corr_opts.attenuate_counts fraction of the data
            if corr_opts.between_sets
                mask1=rand(num_counts_shot_1,1)<corr_opts.attenuate_counts;
                mask2=rand(num_counts_shot_2,1)<corr_opts.attenuate_counts;
                shot_txy_1=shot_txy_1(mask1,:);
                shot_txy_2=shot_txy_2(mask2,:);
                num_counts_shot_1=size(shot_txy_1,1);
                num_counts_shot_2=size(shot_txy_2,1);
            else
                mask=rand(num_counts_shot,1)<corr_opts.attenuate_counts;
                shot_txy=shot_txy(mask,:);
                num_counts_shot=size(shot_txy,1); %recalulate the number of counts
            end
        end
        if corr_opts.between_sets
            pairs_count(shotnum)=num_counts_shot_1*num_counts_shot_2;
            num_counts_shot = num_counts_shot_1+1;
        else
            % full number of pairs in the shot
            pairs_count(shotnum)=num_counts_shot^2 -num_counts_shot;
        end
        
        %         if num_counts_shot<2
        %             warning('%u counts input\n',num_counts_shot)
        %         end
        for ii=1:num_counts_shot-1
            %if each shot is sorted in time then this will only produce counts that are one direction in time
            % to deal with this we need to add both the delta and the -ve of the delta to the difference hist
            if corr_opts.do_pre_mask && ~corr_opts.between_sets  %pre mask optimzation using sortd masking
                % gives a small speedup ~%20 depending on the numbe of counts outsize the correlation range of interest
                % it works by using a sorted dimension to pre-window using binary search
                temp_1d_diff=shot_txy(ii+1:num_counts_shot,corr_opts.sorted_dir)...
                    -delta_multiplier*shot_txy(ii,corr_opts.sorted_dir);
                mask_idx=fast_sorted_mask(...
                    temp_1d_diff,...
                    corr_opts.one_d_window(corr_opts.sorted_dir,1),...
                    corr_opts.one_d_window(corr_opts.sorted_dir,2));
                delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+mask_idx(1):ii+mask_idx(2),:);
            else
                if corr_opts.between_sets
                    if isequal(shot_txy_1,shot_txy_2)
                        if ii == num_counts_shot-1
                            continue
                        end
                        delta=shot_txy_1(ii,:)-delta_multiplier*shot_txy_2(ii+1:(num_counts_shot-1),:);
                    else
                        delta=shot_txy_1(ii,:)-delta_multiplier*shot_txy_2(:,:);
                    end
                else
                    delta=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:);
                end
            end
            one_d_mask_pos=true(size(delta,1),1); %initalize the mask
            one_d_mask_neg=one_d_mask_pos;
            %run the mask in all the queued dimensions
            for mask_dim=mask_dimensions
                one_d_mask_pos=one_d_mask_pos & delta(:,mask_dim)>corr_opts.one_d_window(mask_dim,1) ...
                    & delta(:,mask_dim)<corr_opts.one_d_window(mask_dim,2);
                one_d_mask_neg=one_d_mask_neg & -delta(:,mask_dim)>corr_opts.one_d_window(mask_dim,1) ...
                    & -delta(:,mask_dim)<corr_opts.one_d_window(mask_dim,2);
            end
            
            %to be strictly accurate we must calaulate things symetricaly
            % using fast histograming (gives speedup for sparse histogram)
            % gives >4x speedup on corr unit testing
            
            
            if corr_opts.between_sets && ~isequal(shot_txy_1,shot_txy_2)
                two_d_bins=two_d_bins+histcn(delta(one_d_mask_pos,corr_opts.two_d_dimensions),corr_opts.two_d_edges{1},corr_opts.two_d_edges{2});
                two_d_bins=two_d_bins+histcn(-delta(one_d_mask_neg,corr_opts.two_d_dimensions),corr_opts.two_d_edges{1},corr_opts.two_d_edges{2});
            else
                two_d_bins=two_d_bins+histcn(delta(one_d_mask_pos,corr_opts.two_d_dimensions),corr_opts.two_d_edges{1},corr_opts.two_d_edges{2});
                two_d_bins=two_d_bins+histcn(-delta(one_d_mask_neg,corr_opts.two_d_dimensions),corr_opts.two_d_edges{1},corr_opts.two_d_edges{2});
            end
            
            % old brute histogram approach
            %one_d_bins=one_d_bins+histcounts(delta(one_d_mask_pos,corr_opts.one_d_dimension),corr_opts.one_d_edges)';
            %one_d_bins=one_d_bins+histcounts(-delta(one_d_mask_neg,corr_opts.one_d_dimension),corr_opts.one_d_edges)';
        end
        if mod(shotnum,update_interval)==0 && corr_opts.print_update
            parfor_progress_imp;
        end
    end%loop over shots

out.centers={(corr_opts.two_d_edges{1}(2:end)+corr_opts.two_d_edges{1}(1:end-1))/2,...
    (corr_opts.two_d_edges{2}(2:end)+corr_opts.two_d_edges{2}(1:end-1))/2};
out.pairs=sum(pairs_count);
% sub_index=[1,2,3];
% sub_index(corr_opts.one_d_dimension)=[];
% one_d_volume=corr_opts.one_d_window(sub_index,1)-corr_opts.one_d_window(sub_index,2);
% one_d_volume=prod(one_d_volume)*(corr_opts.one_d_edges(2:end)-corr_opts.one_d_edges(1:end-1));
if corr_opts.between_sets
    two_d_bins = two_d_bins./2;
end
if isfield(corr_opts,'normalisation_factor')
    norm_factor = corr_opts.normalisation_factor;
else
    norm_factor = 1;
end

out.two_d_corr_density_raw=two_d_bins./(shots*norm_factor);
%smooth the correlation function for better normalization
if ~(isnan(corr_opts.one_d_smoothing) || corr_opts.one_d_smoothing==0)
    out.two_d_corr_density=gaussfilt(out.centers,out.two_d_corr_density_raw,corr_opts.one_d_smoothing);
else
    out.two_d_corr_density=out.two_d_corr_density_raw;
end


end


% %%you can check that the sorted mask this is equivelent to the below brute approach
% delta_brute=shot_txy(ii,:)-delta_multiplier*shot_txy(ii+1:num_counts_shot,:);
% min_lim=corr_opts.one_d_window(corr_opts.sorted_dir,1);
% max_lim=corr_opts.one_d_window(corr_opts.sorted_dir,2);
% mask=delta_brute(:,corr_opts.sorted_dir)>min_lim & delta_brute(:,corr_opts.sorted_dir)<max_lim;
% delta_brute =delta_brute(mask,:);
% if ~isequal(delta_brute,delta)
%     fprintf('warning mask is not giving correct answer')
% end