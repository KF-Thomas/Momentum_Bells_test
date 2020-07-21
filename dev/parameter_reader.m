function params = parameter_reader(line_num)
%parameter read
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '';
opts.import.dir = fullfile(opts.data_root, data_folder);
opts.paramfile = fullfile(opts.import.dir, 'log_OptimiserParams.txt');

trip = true;
while trip
    %check for the right line
    try
        % read in log
        logs = readmatrix(opts.paramfile);
        params = logs(line_num,3:end);
        trip = false;
    catch
        %if fail wait
        pause(2.0)
    end
end
end