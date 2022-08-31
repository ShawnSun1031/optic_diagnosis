clear all
data_input = load(fullfile('ANN_train/train.txt'));
training_data_input = [data_input(1:500,:);data_input(701:2808,:)];
data_output = load(fullfile('ANN_ground/ground.txt'));
training_data_output = [data_output(1:500,:);data_output(701:2808,:)];
validation_data=cell(1,2);
validation_data{1} = data_input(501:600,:)';
validation_data{2} = data_output(501:600,:)';
testing_data = cell(1,2);
testing_data{1} = data_input(601:700,:)';
testing_data{2} = data_output(601:700,:)';

options = trainingOptions('adam', ...
        'Shuffle','every-epoch', ...
        'ExecutionEnvironment','gpu', ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.95, ...
        'LearnRateDropPeriod',10, ...
        'ValidationFrequency',1, ...
        'ValidationData',validation_data, ...
        'MaxEpochs',1000, ...
        'ValidationPatience',20, ...
        'MiniBatchSize',128, ...
        'VerboseFrequency',1, ...
        'Plots','training-progress');
    
% numFeatures = 26;
% numHiddenUnits = 125;
% numResponses = 1;
% 
% layers = [ ...
%     sequenceInputLayer(numFeatures)
%     lstmLayer(numHiddenUnits,'OutputMode','sequence')
%     fullyConnectedLayer(numResponses)
%     regressionLayer];

    layers=[sequenceInputLayer(30,'Name','input')
%         fullyConnectedLayer(650)
%         fullyConnectedLayer(1024)
%         leakyReluLayer
        fullyConnectedLayer(1024,'Name','FC_1')
        reluLayer('Name','Relu_1')
    %     reluLayer
        
        fullyConnectedLayer(512,'Name','FC_2')
        reluLayer('Name','Relu_2')
    %     reluLayer
        fullyConnectedLayer(256,'Name','FC_3')
        reluLayer('Name','Relu_3')
    %     reluLayer
        fullyConnectedLayer(128,'Name','FC_4')
        reluLayer('Name','Relu_4')
    %     reluLayer
        fullyConnectedLayer(26,'Name','Output')
        regressionLayer('Name','Regression')];

%     gpuDevice(1);
    [net,info]=trainNetwork(training_data_input',training_data_output',layers,options);
    mkdir('ANN_model')
    output_dir = 'ANN_model';
    save(fullfile(output_dir,['ANN_model.mat']),'net');
    
    lgraph = layerGraph(layers);
    plot(lgraph)
    
    number = 30
    output = predict(net,testing_data{1,1});
    
    RMSD =  sqrt(sum((exp(-output) - exp(-testing_data{1,2})).^2,1)/size(output,1));
    y_mean = mean(exp(-testing_data{1,2}));
    ANN_output_CV = RMSD./y_mean;
    
    RMSD =  sqrt(sum((exp(-testing_data{1,1}(5:30,:)) - exp(-testing_data{1,2})).^2,1)/size(output,1));
    y_mean = mean(exp(-testing_data{1,2}));
    ANN_input_CV = RMSD./y_mean;
    
    input_mean_CV = mean(ANN_input_CV);
    input_std_CV = std(ANN_input_CV);
    output_mean_CV = mean(ANN_output_CV);
    output_std_CV = std(ANN_output_CV);
%     CV_result = [input_mean_CV,input_std_CV,output_mean_CV,output_std_CV];
%     save('CV_result.txt','CV_result','-ascii','-tabs')
    x = 1:1000;
    figure('Renderer', 'painters', 'Position', [10 10 1600 900])
    hold on
    plot(x,info.TrainingLoss,'b')
    plot(x,info.ValidationLoss,'r')
    ylabel('loss')
    xlabel('Epoch')
    legend('Train Loss' ,'Validation Loss')
    
    x = 1:1000;
    figure('Renderer', 'painters', 'Position', [10 10 1600 900])
    hold on
    plot(x,info.TrainingRMSE,'b')
    plot(x,info.ValidationRMSE,'r')
    ylabel('RMSE')
    xlabel('Epoch')
    legend('Train RMSE' ,'Validation RMSE')

    x = 1:100;
    figure('Renderer', 'painters', 'Position', [10 10 1600 900])
    hold on
    title('Each training set CV')
    plot(x,100*ANN_output_CV,'r')
    plot(x,100*ANN_input_CV,'b')
    ylabel('CV(%)')
    xlabel('Train data number #')
    legend('ANN pred CV' ,'Low photon CV')
    
    figure('Renderer', 'painters', 'Position', [10 10 1600 900])
    subplot(1,3,1);plot(0:2*10^-10:5*10^-9,exp(-testing_data{1,1}(5:30,number)));title('ANN Input(low photon)');xlabel('time(s)');ylabel('reflactance');
    subplot(1,3,2);plot(0:2*10^-10:5*10^-9,exp(-output(:,number)));title('ANN Output(pred)');xlabel('time(s)');ylabel('reflactance');
    subplot(1,3,3);plot(0:2*10^-10:5*10^-9,exp(-testing_data{1,2}(:,number)));title('Ground Truth(high photon)');xlabel('time(s)');ylabel('reflactance');
    
    x = 1:2;
    y = [30 95];
    subplot(1,2,1);bar(x,y);ylabel('Time(hours)');
    set(gca, 'xticklabel', {'Low photon','High photon'});
    
    x = 1:2;
    y = [27.45 12.43];
    subplot(1,2,2);bar(x,y);ylabel('CV(%)');
    set(gca, 'xticklabel', {'Low photon','Denoise ANN'});
    
  
    