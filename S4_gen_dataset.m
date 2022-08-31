%{
get the ANN training data and ground truth data
CS Sun
Last update: 2022/08/31
%}

clear all
counts = zeros(2808,26);
subject = "ZJ";
mode = ["ground","train"];
sim_index_set=load('thisPC_sim_wl_index.txt');
sds = 6;

do_normalize = 0;

for m=1:size(mode,2)
    folder_name = fullfile(subject,mode(m));
    mus_table = load(fullfile(subject,mode(m),'mus_table.txt'));
    if do_normalize==1
        mean_mus_table = mean(mus_table/10);
        std_mus_table = std(mus_table/10);
        for i=sim_index_set(1):sim_index_set(2)
            filename = fullfile(folder_name,['sim_' int2str(i)],'cfg_1.mat');
            load(filename)
            filename = fullfile(folder_name,['sim_' int2str(i)],'PL_1.mat');
            load(filename)
            detp.ppath = 10*SDS_detpt_arr{sds};
            photon_weight = each_photon_weight_arr(sds);
            tof=mcxdettime(detp,cfg.prop);
            [tempcounts, idx]=histc(tof,0:cfg.tstep:cfg.tend);
            tempcounts = tempcounts';
            norm_prop = (cfg.prop(2:5,2)' - mean_mus_table(1:4))./std_mus_table(1:4);
            counts(i,:) = [norm_prop,-log((tempcounts+1)/photon_weight)];
        end
    else
        for i=sim_index_set(1):sim_index_set(2)
            filename = fullfile(folder_name,['sim_' int2str(i)],'cfg_1.mat');
            load(filename)
            filename = fullfile(folder_name,['sim_' int2str(i)],'PL_1.mat');
            load(filename)
            detp.ppath = 10*SDS_detpt_arr{sds};
            photon_weight = each_photon_weight_arr(sds);
            tof=mcxdettime(detp,cfg.prop);
            [tempcounts, idx]=histc(tof,0:cfg.tstep:cfg.tend);
            tempcounts = tempcounts';
            figure('Name',mode(m),'NumberTitle','off');
            plot(0:cfg.tstep:cfg.tend,tempcounts(i,:))
            title(mode(m))
            counts(i,:) = -log((tempcounts+1)/photon_weight);

        end
    end
    mkdir('ANN_'+ mode(m))
    save(fullfile('ANN_'+ mode(m),mode(m)+'.txt'),'counts','-ascii','-tabs')
    % plot(0:cfg.tstep:cfg.tend,counts(i,:))
end

