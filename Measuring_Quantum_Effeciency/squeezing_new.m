function out=squeezing_new(halo_centered_cells,plot_on,zones_azm, random_throw_away_perc)
% [1] paper Sub-Poissonian Number Diff in 4 wave mixing 
% debug = true;


%  squeezing_new(halo{1}.counts_vel',true,10);

% print = @(str) if debug fprintf(str) end;

% mirror_azm = true;

range_azm = [0 2]*pi;
range_elev = [0 1]*pi;
% range_elev = [0.25 0.75]*pi;

if ~exist('zones_azm', 'var')
    zones_azm = 4;
end

if ~exist('random_throw_away_perc', 'var')
    random_throw_away_perc = 0;
end 

% zones_azm = 4;
zones_elev = 2;
assert(mod(zones_azm,2)==0, "Need need to be even for the correlated zones to work!");
assert(mod(zones_elev,2)==0, "Need need to be even (I think?) for the correlated zones to work!");
zones_azm_half = zones_azm/2;
zones_total = zones_azm * zones_elev;

halo_centered_cells=halo_centered_cells(~cellfun('isempty',halo_centered_cells));
shots_total = size(halo_centered_cells,2); % Ns in [1]

bin_width_azm  = range(range_azm)/zones_azm;
bin_width_elev = range(range_elev)/zones_elev; 
 
% bin_centers_azm = linspace(0, 2*pi, steps+1) - 0.5*(2*pi)/steps
% bin_centers_azm=wrapTo2Pi(linspace(range_azm(1),range_azm(2),zones_azm+1)-0.5*range(range_azm)/zones_azm);
% bin_centers_azm=bin_centers_azm(2:end);
% bin_boundary_azm=transpose([wrapTo2Pi(bin_centers_azm - 0.5*bin_width_azm) ; wrapTo2Pi(bin_centers_azm + 0.5*bin_width_azm)]);

% if mirror_azm %TODO: what is this doing?
%     bin_centers_azm=[bin_centers_azm, bin_centers_azm+pi];
%     bin_pairs_azm=[bin_pairs_azm;bin_pairs_azm+pi];
% end

% bin_centers_elev=wrapTo2Pi(linspace(range_elev(1),range_elev(2),zones_elev+1)-0.5*range(range_elev)/zones_elev);
% bin_centers_elev=bin_centers_elev(2:end);
% bin_boundary_elev=transpose([wrapTo2Pi(bin_centers_elev - 0.5*bin_width_elev) ; wrapTo2Pi(bin_centers_elev + 0.5*bin_width_elev)]);

zones_ae_index = 1:(zones_total); % so that zones_ae_index(j) = j for all j = 1,2,...,zones_total
zones_ae_index_partner = zeros(1,zones_total);
for iz = 1:zones_total
    elev = fix((iz-1)/zones_azm)+1;      % because matlab start counting from 1  
    azm = mod((iz-1), zones_azm)+1; % because matlab start counting from 1
    azm_partner = mod(azm + zones_azm_half-1, zones_azm)+1; % because matlab start counting from 1
    elev_partner = zones_elev + 1 - elev; % because matlab start counting from 1
    zones_ae_index_partner(iz) = azm_partner + (elev_partner-1)*zones_azm;  % because matlab start counting from 1
    % therefore I hate matlab
end 



map_each_shot_to_bin_ae_id = cell(shots_total,1);
binned_hits_for_all_shot = cell(zones_total,1);
for is = 1:shots_total
    this_halo_shot = halo_centered_cells{is};
    random_throw_away_ind = rand(size(this_halo_shot,1),1) > random_throw_away_perc;
    this_halo_shot(random_throw_away_ind, :);

    [halo_radial,halo_azm,halo_elev]=ConvToSph(this_halo_shot);

    bin_a_s = fix(halo_azm/bin_width_azm)+1;
    bin_e_s = fix(halo_elev/bin_width_elev);

    bin_ae_id = bin_a_s + bin_e_s*zones_azm;

    map_each_shot_to_bin_ae_id{is} = bin_ae_id;
    
    for ihit = 1:size(this_halo_shot,1)
        bin_id_hit = bin_ae_id(ihit);
        binned_hits_for_all_shot{bin_id_hit} = [this_halo_shot(ihit,:); binned_hits_for_all_shot{bin_id_hit}];
    end 
end 

map_each_shot_N = zeros(shots_total, zones_total);
for is = 1:shots_total
    this_halo_shot_bin_ae_id = map_each_shot_to_bin_ae_id{is};
    for in = 1:zones_total
        map_each_shot_N(is,in) = sum(this_halo_shot_bin_ae_id == in);
    end 
end 


index_ij = nchoosek(1:1:zones_total, 2);
index_ij_total = size(index_ij,1);

V_ij = zeros(index_ij_total,1);
V_ij_std = zeros(index_ij_total,1);
V_ij_ind = 1:index_ij_total;

V_ij_corr = false(index_ij_total,1);

for ij = 1:index_ij_total
    
    i = index_ij(ij,1);
    j = index_ij(ij,2);

    if j == zones_ae_index_partner(i) 
        V_ij_corr(ij) = true;
    end 
    if i == zones_ae_index_partner(j)
        V_ij_corr(ij) = true;
    end 
    % assert(i == zones_ae_index_partner(j) & j ~= zones_ae_index_partner(i), "I didn't think is possible to land here?")
    
    Ni = map_each_shot_N(:,i);
    Nj = map_each_shot_N(:,j);

    V_ij(ij) = var(Ni-Nj) / (mean(Ni) + mean(Nj));    % Eq(1) of [1]
    V_ij_std(ij) = (mean((Ni-Nj).^4) - mean((Ni-Nj).^2).^2) / mean(Ni+Nj) / shots_total;  
end 

% debug_ij = [index_ij V_ij V_ij_corr];


if plot_on

    mean_corr   = mean(V_ij(V_ij_corr));
    mean_corr_var = mean(V_ij_std(V_ij_corr).^2);
    mean_corr_std = sqrt(mean_corr_var);
    
    mean_uncorr = mean(V_ij(~V_ij_corr));
    mean_uncorr_var = mean(V_ij_std(~V_ij_corr).^2);
    mean_uncorr_std = sqrt(mean_uncorr_var);

    cb = [0 0 1];
    cr = [1 0 0];
    figure(200);clf(200);
    figure(200);
    errorbar(V_ij_ind(~V_ij_corr), V_ij(~V_ij_corr), V_ij_std(~V_ij_corr), ...
        '.', 'Color',[cb,.9],'CapSize',0); hold on;
    errorbar(V_ij_ind( V_ij_corr), V_ij( V_ij_corr), V_ij_std( V_ij_corr), ...
        '.', 'Color',[cr,.9],'CapSize',0); hold on;
    yline(mean_uncorr, '-',num2str(mean_uncorr,4)+"$\pm$"+num2str(mean_uncorr_std,4), ...
        'LineWidth',0.5,'Color',[cb,.9],'Interpreter','latex'); hold on; 
    yline(mean_corr,   '-',num2str(mean_corr,  4)+"$\pm$"+num2str(mean_corr_std,  4), ...
        'LineWidth',0.5,'Color',[cr,.9], 'Interpreter','latex'); hold on; 
    rectangle('Position',[0,mean_uncorr-mean_uncorr_std index_ij_total*1.05 2*mean_uncorr_std], ...
        'LineStyle','none','FaceColor',[cb,.1]); hold on;
    rectangle('Position',[0,mean_corr-mean_corr_std   index_ij_total*1.05 2*mean_corr_std], ...
        'LineStyle','none','FaceColor',[cr,.1]); hold on;
    xlabel('Zone pair $ij$');
    ylabel('Normalised variance $V_{ij}$');
    xlim([0, index_ij_total*1.02]);
    legend('uncorrelated', 'correlated');
    hold off;
    shg;

    figure(201);clf(201);
    figure(201);
    assert(zones_elev==2, "below code only works with two elev zones! ");
    colors = [hsv(zones_total/2); zeros(zones_azm, 3)];
    for ie = zones_azm:zones_total
        colors(ie,:) = colors(zones_ae_index_partner(ie),:);
    end 
    
    for ibin = 1:zones_total
        bhs = binned_hits_for_all_shot{ibin}; % binned hits 
        c = colors(ibin,:);
        plot3(bhs(:,1), bhs(:,2), bhs(:,3), '.', 'Color',[c 0.1], 'LineStyle','none'); hold on;
        
    end
    xlabel("x");
    xlabel("y");
    xlabel("z");
    title("Correlated zones");
    set(gca, 'Projection','perspective');
    pbaspect([1 1 1]);
    view(300, 100);
%     set(gca, 'rat')
    shg;

end 

% V_ij_lin = V_ij(:); 
out = [V_ij, V_ij_std, V_ij_corr];

end









