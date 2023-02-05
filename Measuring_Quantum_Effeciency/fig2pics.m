%
%
% 
% Another reason I hate matlab
%
%
% Workflow: edit in macOS, generate figure in macOS, 
%   then to generate png: remote desktop into intel mac, then turn on windows simulation, run matlab there
%
%



assert(exist("fig_output_path", "var"), "Run squeezing_new_run_everything.m first!");





% file_list = ls(fig_output_path);

dir_fig = dir(fullfile(fig_output_path, '*.fig'));
dir_png = dir(fullfile(fig_output_path, '*.png'));
dir_pdf = dir(fullfile(fig_output_path, '*.pdf'));

files_count = size(dir_fig,1);

% dir_fig_names = [dir_fig.name];
% dir_pdf_names = [dir_pdf.name];
% dir_png_names = [dir_png.name];

last_size_fprintf = 0;
for id = 1:files_count
    fprintf(repmat('\b', 1, last_size_fprintf));
    last_size_fprintf = fprintf("calculating " + num2str(id) + " out of " + num2str(files_count));
    
    fig_path = fig_output_path + dir_fig(id).name;
    [filepath, name, ext] = fileparts(fig_path);
     
    fig = openfig(fig_path);

    if ~exist(fig_output_path+name+".pdf", "file")
%         fprintf(name+'.pdf'+'\n');
        saveas(fig,fig_output_path+name+'.pdf');
    end

    if ~exist(fig_output_path+name+".png", "file")
%         fprintf(name+'.pnd'+'\n');
        print(fig, fig_output_path+name+'.png', '-dpng', '-r300');
    end
    close(fig);
    pause(0.5);
end




