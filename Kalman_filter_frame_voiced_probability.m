%% 
% PhD Project, 2015-2019 
% Dionelis Nikolaos, CID: 00690438 

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Clear all the previous data
clear all; clc; 
close all;

% Add the datapath of the voicebox
addpath ./voicebox

% Set up voicebox
%y_voicebox = voicebox('speech_voicebox');
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

addpath ./clean

addpath ./babble/0dB
addpath ./babble/5dB
addpath ./babble/10dB
addpath ./babble/15dB

addpath ./airport/0dB
addpath ./airport/5dB
addpath ./airport/10dB
addpath ./airport/15dB

addpath ./car/0dB
addpath ./car/5dB
addpath ./car/10dB
addpath ./car/15dB
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% % Add the datapaths of the "stprtool" toolbox
% addpath ./External_stprtool
% % Use "stprpath" to add the datapaths of the Statistical Pattern Recognition Toolbox
% run('External_stprtool/stprpath.m'); 

% % Add the datapaths of the noise databases
% addpath ./Noise_Database_44kHz
% addpath ./Noise_Database_8kHz
% addpath ./Noise_Database_from_RTDSP_course
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Add the datapaths of the TIMIT database
% addpath ./TIMIT\TRAIN\DR1\FCJF0
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------
% 
% % Add the datapath of the Nato Noise database
% addpath ./NatoNoise0
% 
% % Add the datapath for binaural noise files
% % this is from the noise P501 database
% addpath ./Binaural
% % Add the datapath for monaural noise files
% % this is from the noise P501 database
% addpath ./Monaural
% %-------------------------------------------------------------------------------------------------------------------------------
% %-------------------------------------------------------------------------------------------------------------------------------

% Add the datapaths of the functions for Noise Estimation
addpath ./Functions_for_Noise_Estimation

% Add the datapaths of the functions for the Enhancement Systems
addpath ./Functions_for_Enhancement_System
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

addpath ./FCJF0
addpath ./NatoNoise0
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

%%
% KF for pv

% Things to do: 
% 1) EM algorithm VS LS/MMSE for training the KF. We use only LS/MMSE here. We need to use the EM algorithm.
% 2) Use future frames. Use KF for smoothing/hindsight and not for filtering. In general, filtering VS prediction VS smoothing.
% 3) Now, pv is the frame voiced probability. We need to try pv(k,l), we need to try fr.bin and frame voiced probability.

%-------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------------

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we store the true values
true_pv = pv;

% we use noise 

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% Read the noise file
[s2, fs2, bits] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

% we find error
error_pv = pv - true_pv;

% fit zero mean gaussian
[m,v,w] = gaussmix(error_pv, [], [], 1);
[m2,v2,w2] = gaussmix(error_pv-m, [], [], 1);

% we use v2
v2

% clear the other temp variables
clear m;
clear v;
clear w;
clear m2;
clear w2;

% this was for one sentence
% we need to repeat with many sentences

% we initialize the counter
counter = 0;

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    total_error_pv = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 5;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        % we find error
        error_pv = pv - true_pv;
        
        total_error_pv = [total_error_pv; error_pv'];
        
        total_total_error_pv{counter} = total_error_pv;
        
    end
end

% this was for 5 dB SNR 
% we now use 0 dB SNR 

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    total_error_pv = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 0;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        % we find error
        error_pv = pv - true_pv;
        
        total_error_pv = [total_error_pv; error_pv'];
        
        total_total_error_pv{counter} = total_error_pv;
        
    end
end

% this was for 0 dB SNR 
% we now use 10 dB SNR 

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    total_error_pv = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 10;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        % we find error
        error_pv = pv - true_pv;
        
        total_error_pv = [total_error_pv; error_pv'];
        
        total_total_error_pv{counter} = total_error_pv;
        
    end
end

% this was for 0 dB SNR 

% we now need one array with all the values

% we initialize the new array
errors_in_pv = [];

for i = 1 : counter
    temp_array = total_total_error_pv{i};
    
    for k = 1 : size(temp_array,1)
        error_in_pv = [errors_in_pv temp_array(k,:)];
    
    end
end

% we use: error_in_pv

error_in_pv = error_in_pv';

% fit zero mean gaussian
[m,v,w] = gaussmix(error_in_pv, [], [], 1);
[m2,v2,w2] = gaussmix(error_in_pv-m, [], [], 1);

% we use v2
v2

% v2 is 0.0942
% we use this value so as to not wait for training
%v2 = 0.0942;

% clear the other temp variables
clear m;
clear v;
clear w;
clear m2;
clear w2;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now fit the transition model 

tr_model = [];

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    tr_model = [tr_model; diff(pv)];

end

% fit zero mean gaussian
[m,v,w] = gaussmix(tr_model, [], [], 1);
[m2,v2_x,w2] = gaussmix(tr_model-m, [], [], 1);

% we use v2_x
v2_x

% v2_x is 0.0198
% we use this value so as to not wait for training
%v2_x = 0.0198;

% clear the other temp variables
clear m;
clear v;
clear w;
clear m2;
clear w2;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now use Kalman filter (KF)

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

reference_true_pv = pv;

% Read the noise file
[s2, fs2] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

predicted_x_m = [];
predicted_P_m = [];

for loop = 1 : length(pv)
    % we initialize x_m and P_m
    if (loop == 1)
        x_m = 0;
        P_m = 10.2;
        
    end
    
    x_m = (((P_m+v2_x)*pv(loop)) + (v2*x_m)) / (P_m+v2_x+v2);
    
    P_m = ((P_m+v2_x) * v2) / (P_m+v2_x+v2);
    
    % store the results
    predicted_x_m = [predicted_x_m x_m];
    predicted_P_m = [predicted_P_m P_m];

end

predicted_x_m = predicted_x_m';
predicted_P_m = predicted_P_m';

figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);

plot(predicted_x_m);
hold on;
plot(pv, 'g');
hold on;
plot(reference_true_pv, 'r');
hold off;

set(gca,'FontSize',10);
title('pv. Probability of frame being voiced.');
xlabel('Time in frames');
ylabel('Probability pv');

legend(' new pv', ' PEFAC pv', ' true pv');

% ------------------------------------------------------------------------

% treat the problem as a classification problem
% we use 0 and 1 only

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) > 0.5)
       reference_true_pv(i) = 1;
   else
       reference_true_pv(i) = 0;
   end
end

for i = 1 : length(pv)
   if (pv(i) > 0.5)
       pv(i) = 1;
   else
       pv(i) = 0;
   end
end

for i = 1 : length(predicted_x_m)
   if (predicted_x_m(i) > 0.5)
       predicted_x_m(i) = 1;
   else
       predicted_x_m(i) = 0;
   end
end

% we find accuracy for pv
pv_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == pv(i))
       pv_accuracy = pv_accuracy + 1;
   end
end

pv_accuracy = pv_accuracy / length(reference_true_pv);

% we find accuracy for x_m
x_m_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == predicted_x_m(i))
       x_m_accuracy = x_m_accuracy + 1;
   end
end

x_m_accuracy = x_m_accuracy / length(reference_true_pv);

% we now use pv
% we compare pv and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

pv_true_positives = 0;
pv_true_negatives = 0;

pv_false_positives = 0;
pv_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && pv(i) == 1)
        pv_true_positives = pv_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 0)
        pv_true_negatives = pv_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && pv(i) == 0)
        pv_false_positives = pv_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 1)
        pv_false_negatives = pv_false_negatives + 1;
    end
end
    
% find precision
pv_precision = (pv_true_positives) / (pv_true_positives+pv_false_positives);

% find recall
pv_recall = (pv_true_positives) / (pv_true_positives+pv_false_negatives);

% we now do the same for x_m

% we now use x_m
% we compare x_m and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

x_m_true_positives = 0;
x_m_true_negatives = 0;

x_m_false_positives = 0;
x_m_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && predicted_x_m(i) == 1)
        x_m_true_positives = x_m_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 0)
        x_m_true_negatives = x_m_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && predicted_x_m(i) == 0)
        x_m_false_positives = x_m_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 1)
        x_m_false_negatives = x_m_false_negatives + 1;
    end
end

% find precision
x_m_precision = (x_m_true_positives) / (x_m_true_positives+x_m_false_positives);

% find recall
x_m_recall = (x_m_true_positives) / (x_m_true_positives+x_m_false_negatives);

% we thus have 3 metrics for the two methods

% we have: pv_accuracy, pv_precision, pv_recall 
% we have: x_m_accuracy, x_m_precision, x_m_recall

% print the metrics
disp(['The pv PEFAC has accuracy: ', num2str(pv_accuracy)]);
disp(['The pv PEFAC has precision: ', num2str(pv_precision)]);
disp(['The pv PEFAC has recall: ', num2str(pv_recall)]);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The new pv method has accuracy: ', num2str(x_m_accuracy)]);
disp(['The new pv method has precision: ', num2str(x_m_precision)]);
disp(['The new pv method has recall: ', num2str(x_m_recall)]);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% we now repeat the KF

% we use a different training method
% we use LS or MMSE training for the KF

% we use non-iterative approach 

% we need to define total_total_observations 
% we need to define the states: total_states

% time increases in the column direction 

counter = 0;

states1 = [];
observations2 = [];

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    states2 = [];
    observations2 = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 5;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        states2 = [states2; true_pv'];
        observations2 = [observations2; pv'];
        
        states1{counter} = states2;
        observations1{counter} = observations2;
        
    end
end

% this was for SNR 5 dB
% we now use SNR 0 dB

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    states2 = [];
    observations2 = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 0;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        states2 = [states2; true_pv'];
        observations2 = [observations2; pv'];
        
        states1{counter} = states2;
        observations1{counter} = observations2;
        
    end
end

% this was for SNR 0 dB
% we now use SNR 10 dB

% Outter for loop
for filename2 = {'SA1.WAV', 'SA2.WAV', 'SI648.WAV', ...
        'SI1027.WAV', 'SI1657.WAV', 'SX37.WAV', ...
        'SX127.WAV', 'SX217.WAV', 'SX307.WAV', ...
        'SX397.WAV'}
    % Read the PHN TIMIT file
    [s,fs,wrd,phn] = readsph(filename2{1},'wt');

    states2 = [];
    observations2 = [];
    
    % Use counter 
    counter = counter + 1;
    
    % we use PEFAC to find voiced speech
    [fx,tx,pv,fv] = fxpefac(s,fs);

    % we store the true values
    true_pv = pv;

    % Inner for loop
    for filename = {'babble.wav', 'buccaneer1.wav', 'buccaneer2.wav', ...
            'destroyerengine.wav', 'destroyerops.wav', 'f16.wav', 'factory1.wav',  ...
            'factory2.wav', 'hfchannel.wav', 'leopard.wav', 'm109.wav', ...
            'machinegun.wav', 'pink.wav', 'volvo.wav', 'white.wav'}
        % Read the noise file
        [s2, fs2] = audioread(filename{1});

        % Define the SNR
        snr = 10;

        % Add noise from column vector n() at sample frequency fn with random start
        % Sample and wrapping around as required with a vorbis cross-fade
        z = v_addnoise(s,fs,snr,'',s2,fs2); 

        s = z;

        % we use PEFAC to find voiced speech
        [fx,tx,pv,fv] = fxpefac(s,fs);

        states2 = [states2; true_pv'];
        observations2 = [observations2; pv'];
        
        states1{counter} = states2;
        observations1{counter} = observations2;
        
    end
end

% SNR 5, 0 and 10 were used
% the training with SNR 5,0 and 10 has finished

% we need to define total_total_observations 
% we need to define the states: total_states

total_states = [];
total_total_observations = [];

temp_array = [];

for i = 1 : counter
    temp_array = states1{i};
    temp_array2 = observations1{i};
    
    for k = 1 : size(temp_array,1)
        for kk = 1 : size(temp_array,2)
            total_states = [total_states temp_array(k,kk)];
            total_total_observations = [total_total_observations temp_array2(k,kk)];
        
        end
    end
end

% KF -- train parameters -- assume constant

% In KF, there is no access to the correct data when the algorithm is running in real time.
% In KF, the error is between the measurement and the Kalman prediction.
% It can be stated that the KF resembles a low pass filter.
% but its transfer characteristic is nonlinear, and the cut-off frequency shifts.

% In KF, there are 2 main things to do in each time step. 
% The first thing is to calculate the apriori state estimate, and the measurement update.
% The second step is to calculate the gain, corresponding covariance, and then the posterior state.

% we use M time steps
M = size(total_states,2);

for loop = 1 : M
    temp = total_states(:,loop);
    
    X(:,loop) = temp';
end

X1 = X(:, 1:(end-1));
X2 = X(:, 2:end);

clear temp;
C = 1;

for loop = 1 : M
    temp = zeros(1,C);

    % temp = observations;
    temp = total_total_observations(:,loop);

    Z(:,loop) = temp';
end

%A = X2 * X1' * (X1 * X1')^(-1);
A = X2 * X1' * pinv(X1 * X1');

%H = Z * X' * (X * X')^(-1);
H = Z * X' * pinv(X * X');

W = (X2 - (A * X1)) * (X2 - (A * X1))' / (M-1);

Q = (Z - (H * X)) * (Z - (H * X))' / M; 

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% do the KF

% A is the same as the F matrix
% W is the same as the Q matrix

F = A;
Q1 = Q;
Q = W;

% H is the same as H
% Q is the same as the R matrix

R = Q1;

% we use noisy speech signal
% we implement the KF filter

% Read the PHN TIMIT file
[s,fs,wrd,phn] = readsph('SA1.WAV','wt');

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

reference_true_pv = pv;

% Read the noise file
[s2, fs2] = readwav('babble.wav');

% Define the SNR
snr = 5;

% Add noise from column vector n() at sample frequency fn with random start
% Sample and wrapping around as required with a vorbis cross-fade
z = v_addnoise(s,fs,snr,'',s2,fs2); 

s = z;

% we use PEFAC to find voiced speech
[fx,tx,pv,fv] = fxpefac(s,fs);

predicted_x_m = [];
predicted_P_m = [];

for loop = 1 : length(pv)
    % we initialize x_m and delta_delta_phi and P_p
    if (loop == 1)
        x_m = 0.5;
        P_m = 1000.2;
    end

    % define l 
    % l is the observation
    l = pv(loop);
    
    K = ((F*P_m*F') + Q) * H' * ((H * ((F*P_m*F') + Q) * H') + R)^(-1);
    
    x_p = (F*x_m) + (K * (l - (H * F * x_m)));
    
    P_p = (1-K) * ((F * P_m * F') + Q);
    
    % constrain x_p
    if (x_p < 0)
       x_p = 0.01; 
    end
    
    % constrain x_p
    if (x_p > 1)
       x_p = 0.99; 
    end
    
    % store the results
    predicted_x_m = [predicted_x_m x_p];
    predicted_P_m = [predicted_P_m P_p];

    % prepare the next iteration
    x_m = x_p;
    P_m = P_p;

end

predicted_x_m = predicted_x_m';
predicted_P_m = predicted_P_m';

figure; 
set(gcf,'Color','w');
set(gca,'FontSize',10);

plot(predicted_x_m);
hold on;
plot(pv, 'g');
hold on;
plot(reference_true_pv, 'r');
hold off;

set(gca,'FontSize',10);
title('pv. Probability of frame being voiced.');
xlabel('Time in frames');
ylabel('Probability pv');

legend(' new pv', ' PEFAC pv', ' true pv');

% ------------------------------------------------------------------------

% treat the problem as a classification problem
% we use 0 and 1 only

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) > 0.5)
       reference_true_pv(i) = 1;
   else
       reference_true_pv(i) = 0;
   end
end

for i = 1 : length(pv)
   if (pv(i) > 0.5)
       pv(i) = 1;
   else
       pv(i) = 0;
   end
end

for i = 1 : length(predicted_x_m)
   if (predicted_x_m(i) > 0.5)
       predicted_x_m(i) = 1;
   else
       predicted_x_m(i) = 0;
   end
end

% we find accuracy for pv
pv_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == pv(i))
       pv_accuracy = pv_accuracy + 1;
   end
end

pv_accuracy = pv_accuracy / length(reference_true_pv);

% we find accuracy for x_m
x_m_accuracy = 0;

for i = 1 : length(reference_true_pv)
   if (reference_true_pv(i) == predicted_x_m(i))
       x_m_accuracy = x_m_accuracy + 1;
   end
end

x_m_accuracy = x_m_accuracy / length(reference_true_pv);

% we now use pv
% we compare pv and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

pv_true_positives = 0;
pv_true_negatives = 0;

pv_false_positives = 0;
pv_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && pv(i) == 1)
        pv_true_positives = pv_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 0)
        pv_true_negatives = pv_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && pv(i) == 0)
        pv_false_positives = pv_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && pv(i) == 1)
        pv_false_negatives = pv_false_negatives + 1;
    end
end
    
% find precision
pv_precision = (pv_true_positives) / (pv_true_positives+pv_false_positives);

% find recall
pv_recall = (pv_true_positives) / (pv_true_positives+pv_false_negatives);

% we now do the same for x_m

% we now use x_m
% we compare x_m and reference_true_pv

% define true positives, true negatives
% define false positives, false negatives

x_m_true_positives = 0;
x_m_true_negatives = 0;

x_m_false_positives = 0;
x_m_false_negatives = 0;

for i = 1 : length(reference_true_pv)
    if (reference_true_pv(i) == 1 && predicted_x_m(i) == 1)
        x_m_true_positives = x_m_true_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 0)
        x_m_true_negatives = x_m_true_negatives + 1;
    elseif (reference_true_pv(i) == 1 && predicted_x_m(i) == 0)
        x_m_false_positives = x_m_false_positives + 1;
    elseif (reference_true_pv(i) == 0 && predicted_x_m(i) == 1)
        x_m_false_negatives = x_m_false_negatives + 1;
    end
end

% find precision
x_m_precision = (x_m_true_positives) / (x_m_true_positives+x_m_false_positives);

% find recall
x_m_recall = (x_m_true_positives) / (x_m_true_positives+x_m_false_negatives);

% we thus have 3 metrics for the two methods

% we have: pv_accuracy, pv_precision, pv_recall 
% we have: x_m_accuracy, x_m_precision, x_m_recall

disp([' ']);
disp([' ']);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The pv PEFAC has accuracy: ', num2str(pv_accuracy)]);
disp(['The pv PEFAC has precision: ', num2str(pv_precision)]);
disp(['The pv PEFAC has recall: ', num2str(pv_recall)]);

disp([' ']);
disp([' ']);

% print the metrics
disp(['The new pv method has accuracy: ', num2str(x_m_accuracy)]);
disp(['The new pv method has precision: ', num2str(x_m_precision)]);
disp(['The new pv method has recall: ', num2str(x_m_recall)]);

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
