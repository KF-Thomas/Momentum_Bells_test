%optimisation algorithim

shot_num = 88;
max_num = 151;

path_input = 'C:\Users\BEC Machine\Documents\MATLAB\exp_input.txt';
path_output = 'C:\Users\BEC Machine\Documents\MATLAB\exp_output.txt';
path_param_log = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_OptimiserParams.txt';
formatSpec = '%s';

previous_params = [0];
params = [0];

while shot_num<max_num
    
    %read input file
    while all(previous_params == params)
        try
            f_input=fopen(path_input,'r');
            A = fscanf(f_input,formatSpec);
            param_str = A(15:end-2);
            params = str2num(param_str);
            fclose(f_input);
        end
        pause(1.0)
    end
    
    %write to log
    f_log=fopen(path_param_log,'a');  % append to param-log-file
    nowdt=datetime('now');
    fprintf(f_log,'%d,%.3f,%s\n',shot_num,posixtime(nowdt),param_str);
    fclose(f_log);
    
    %calculate cost of last run
    %     cost = magnetic_transfer_cost(shot_num);
    cost = momentum_transfer_cost(shot_num);
    
    %write the cost of the next file
    cost_val = cost.val;
    bad = 'False';
    if isempty(cost) || isempty(cost.val) || isnan(cost.val) || isinf(cost.val)
        cost_val = 10000000;
        unc = 100000;
        bad = 'True';
    end
    %cost write
    fileID = fopen(path_output,'w');
    fprintf(fileID,'cost = %4.4f\n',cost_val);
    if isfield(cost,'unc')
        unc = cost.unc;
        fprintf(fileID,'uncer = %4.4f\n',unc);
    end
    fprintf(fileID,'bad = %s',bad);
    fclose(fileID);
    
    shot_num = shot_num+1;
    previous_params = params;
end