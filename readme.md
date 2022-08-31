# optic diagnosis final project
>purpose : naive idea to using high photons Monte Carlo result as target, and low photons Monte Carlo target as train set.
---
## Simulation Steps
1. Put the segmented model and position/direction of the probes in the `models` folder.  
    ![](https://i.imgur.com/ol0sWyu.png)
2. Use the `S1_make_the_sim_setting.m` to make the setting, including:  
    1. how many SDS are there to simulate.
    2. which mus combination to simulation.
    3. choose ground / train mode to match each simulation photon number

    We only need to set the mus to simulate, because the forward is using WMC, so the mua can be any combination.  
    For example, the mus for each layer are setting as below:  
    ![](https://i.imgur.com/OdeYvgx.png)  
    Then there will be 13X9X4X6=2808 sets of mus combinations. And the program will auto generate these combinations for you. 
    
    ```matlab
    %% param
    model_folder='models'; % the folder of the models
    target_name_arr={'ZJ'}; % the name of the model to simulate
    num_SDS=7; % number of detectors
    layer_mus={[50:25:350],[50:37.5:350],[10 19 28 37],[50:60:350]}; % mus for each layer, 1/cm
    mode = ["ground","train"];  % [ground/train]
    num_photon=[1E9,1E8]; % the number of photon
    n=1.4;
    g=0.9;
    sim_version=4.41;
    ```
3. Use 'thisPC_sim_wl_index.txt' to set the beginning and endding index for wavelength to simulate by this PC
    ```
    % The beginning and endding index for wavelength to simulate by this PC
    1 2808
    ```

4. Use `S2_run_script.m` to run the simulation.  
    You can set it to run many simulation one-by-one.  
    ![](https://i.imgur.com/Tdp4faW.png)  
    If you hany more than one GPU on your computer, you should set the `GPU_setting.txt` to determine which GPU is used and the load for each GPU.  
    ![](https://i.imgur.com/YEHCKdz.png)  
5. The result folder for each subject will containe many [sim_ + index] folders  
    ![](https://i.imgur.com/VeaoCSL.png)  
    In each folder is the WMC simulation result of the given mus combination. It records the pathlength in each layer for each detected photon.  
    There will also be some files containing the information of the lookup table, e.g., the `mus_table.txt` containing the mus for each set of combination.  
    ![](https://i.imgur.com/EvFX3tw.png)
6. After the simulation, use `S3_exam_the_simed_table.m` to check if there is any error in the simulation result. 
  
7. After checking simulation is done, use `S4_gen_dataset.m` to generate both `ANN_ground` data and `ANN_train` data.
    ![](https://i.imgur.com/cZWPqKL.jpg)

8. Run `S5_ANN_train.m` to train our ANN model.
    ![](https://i.imgur.com/PzvTYey.jpg)
    ![](https://i.imgur.com/2jfbyfs.jpg)
---
TODO
- [ ] Using more mua set to get make model more robustness. 
