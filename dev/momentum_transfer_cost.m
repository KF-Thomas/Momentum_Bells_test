%cost function
function cost = momentum_transfer_cost(shot_num)

frac_opts.num_lim = 0.05e3;
frac_opts.transfer_state = 'momentum';
frac_opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds

anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=true;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];
anal_opts.tdc_import.shot_num=shot_num;
data_read_check = true;
while data_read_check
    try
        batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
        if ~isnan(batch_data.mcp_tdc.shot_num)
            data_read_check = false;
        else
            pause(1.0)
        end
    catch
        pause(1.0)
    end
end
out_frac = fraction_calc(batch_data.mcp_tdc,frac_opts);
%k=0,-1,-2 cost
% cost.val = abs(out_frac.Ns(:,1)-out_frac.Ns(:,2))./(out_frac.Ns(:,1)+out_frac.Ns(:,2))...
%     + abs(out_frac.Ns(:,3)-out_frac.Ns(:,2))./(out_frac.Ns(:,3)+out_frac.Ns(:,2))...
%     + abs(out_frac.Ns(:,1)-out_frac.Ns(:,3))./(out_frac.Ns(:,1)+out_frac.Ns(:,3));%...
    %+ abs(out_frac.Ns(:,4)-(out_frac.Ns(:,1)+out_frac.Ns(:,2)+out_frac.Ns(:,3)))./out_frac.Ns(:,4);

%k=-1,-2 cost
% cost.val = abs(out_frac.Ns(:,1)-out_frac.Ns(:,2))./(out_frac.Ns(:,1)+out_frac.Ns(:,2))...
%     + abs(out_frac.Ns(:,4)-(out_frac.Ns(:,1)+out_frac.Ns(:,2)))./out_frac.Ns(:,4);
% 
% %k=0,-1 cost
% cost.val = abs(out_frac.Ns(:,3)-out_frac.Ns(:,2))./(out_frac.Ns(:,3)+out_frac.Ns(:,2))...
%     + abs(out_frac.Ns(:,4)-(out_frac.Ns(:,3)+out_frac.Ns(:,2)))./out_frac.Ns(:,4);

%k=0,-2 cost
cost.val = abs(out_frac.Ns(:,1)-out_frac.Ns(:,3))./(out_frac.Ns(:,1)+out_frac.Ns(:,3))...
    + abs(out_frac.Ns(:,4)-(out_frac.Ns(:,1)+out_frac.Ns(:,3)))./out_frac.Ns(:,4);


%cost.unc = sqrt((1+out_frac.fracs(:,2))./out_frac.Ns(:,2));

%make sure theres enough number in the states
if ~isempty(out_frac.Ns(:,1:3))
    cost.val(sum(out_frac.Ns(:,1:3),2)<frac_opts.num_lim) = nan;
end
end