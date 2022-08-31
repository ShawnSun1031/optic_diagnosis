%{
Run the lookup table simulation

Benjamin Kao
Last update: 2020/10/26
%}

clc;clear;close all;

%% clear the stop flag
to_save=0;
save('stop_flag.txt','to_save','-ascii','-tabs');
mode = ["ground", "train"];
subject = "ZJ";
for i=1:size(mode,2)
    folderpath = fullfile(subject,mode(i))
    fun_MCX_run_lookup(folderpath,subject);
end
% fun_MCX_run_lookup('WW');