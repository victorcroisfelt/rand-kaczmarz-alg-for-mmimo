% This Matlab script can be used to generate the Figure 3 in the article:
%
% Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão,
% "Randomized Kaczmarz Algorithm for Massive MIMO Systems with Channel
% Estimation and Spatial Correlation", to appear, 2019.
%
% Download article:
%
%                   https://arxiv.org/abs/1904.04376
%
% This is version 5.0 (Last edited: 2019-15-07)
%
% License: This code is licensed under the MIT license. If you in any way
% use this code for research that results in publications, please reference
% our original article as shown above.
%
% References:
% [1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive
% MIMO Networks: Spectral, Energy, and Hardware Efficiency", Foundations
% and Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI:
% 10.1561/2000000093
%

% Initialization
close all;
clearvars;

%% Simulation scenarios

% Number of BS antennas
M = 100;

% Range of the number of UEs
Krange = [10 30 50];

% Extract the max number of UEs
Kmax = max(Krange);

% Correlation factor of the exponential correlation model
correlationFactor = 0.5;

% Standard deviation of large-scale fading (shadowing) [dB]
stdLSF = 4;

% Number of rKA iterations
numIterationsRange = [1 100:100:6000];

% Number of setups
numSetups = 10;

% Number of stats ergodic realizations
numRealizations = 100;

%% Propagation parameters

% Communication bandwidth [Hz]
B = 20e6;

% Total uplink transmit power per UE [mW]
p = 100;

% Define noise figure at BS [dB]
noiseFigure = 10;

% Compute total noise power [dBm]
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

% Select length of coherence block [samples]
tau_c = 200;

%% Simulation

% Prepare to save simulation results
meanSE_RZF  = zeros(3,2,length(Krange),numSetups);
meanSE_hyb  = zeros(length(numIterationsRange),3,2,length(Krange),numSetups);

% Go through all setups
for s = 1:numSetups
    
    % Setting a timer
    setup = tic;
    
    % Output setups progress
    disp([num2str(s) ' setups out of ' num2str(numSetups)]);
    
    % Compute the setup characteristics
    [Runcorr,Rcorr,channelGaindB] = functionSetup(M,Kmax,correlationFactor,stdLSF);
    
    % Compute the normalized channel gains
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    % Go through all K values
    for k = 1:length(Krange)
        
        % Extract the current value of UEs
        K = Krange(k);
        
        % Select uncorrelated or correlated fading channels
        for chn = 1:2
            
            if chn == 1
                
                % Extract the current covariance matrices
                RuncorrMatrix = Runcorr(:,:,1:K);
                
                % Compute the channel estimates for LS/MMSE channel estimators
                [H,Hhat_LS,Hhat_MMSE,tau_p] = functionChannelEstimates(M,K,p,RuncorrMatrix,channelGainOverNoise,numRealizations);
                
            else
                
                % Extract the current covariance matrices
                RcorrMatrix = Rcorr(:,:,1:K);
                
                % Compute the channel estimates for LS/MMSE channel estimators
                [H,Hhat_LS,Hhat_MMSE,tau_p] = functionChannelEstimates(M,K,p,RcorrMatrix,channelGainOverNoise,numRealizations);
                
            end
            
            % Compute the UL SE of the RZF canonical scheme
            [SE_RZF,SE_RZF_LS,SE_RZF_MMSE] = functionComputeSE_UL_Canonical(M,K,p,tau_c,tau_p,numRealizations,H,Hhat_LS,Hhat_MMSE);

            % Save simulation results
            meanSE_RZF(1,chn,k,s) = mean(SE_RZF);
            meanSE_RZF(2,chn,k,s) = mean(SE_RZF_LS);
            meanSE_RZF(3,chn,k,s) = mean(SE_RZF_MMSE);
            
            % Run PARL rKA-based schemes with the hybrid initialization 
            % approach
            [V_hyb] = functionPARLrKA(M,K,p,numRealizations,numIterationsRange,H,Hhat_LS,Hhat_MMSE);
        
            % Compute the UL SEs for the hybrid approach
            [SE_hyb_true,SE_hyb_LS,SE_hyb_MMSE] = functionComputeSE_UL_PARLrKA(M,K,p,tau_c,tau_p,numRealizations,numIterationsRange,H,V_hyb);

            % Save simulation results
            meanSE_hyb(:,1,chn,k,s) = mean(SE_hyb_true,1);
            meanSE_hyb(:,2,chn,k,s) = mean(SE_hyb_LS,1);
            meanSE_hyb(:,3,chn,k,s) = mean(SE_hyb_MMSE,1);
            
        end
        
    end
   
    toc(setup)
    
end

%% Plot simulation results

figure;
hold on; box on; ax = gca;

% Extracting average points
y(:,1,:,:) = 100*(abs(mean(meanSE_hyb(:,1,:,:,:),4)-mean(meanSE_RZF(1,:,:,:),4))./mean(meanSE_RZF(1,:,:,:),4));
y(:,2,:,:) = 100*(abs(mean(meanSE_hyb(:,2,:,:,:),4)-mean(meanSE_RZF(2,:,:,:),4))./mean(meanSE_RZF(2,:,:,:),4));
y(:,3,:,:) = 100*(abs(mean(meanSE_hyb(:,3,:,:,:),4)-mean(meanSE_RZF(3,:,:,:),4))./mean(meanSE_RZF(3,:,:,:),4));

% Legend plots
plot(numIterationsRange(1),y(1,1,1,1),'k--'); % Uncorrelated
plot(numIterationsRange(1),y(1,1,2,1),'k-.'); % Correlated
plot(numIterationsRange(1),y(1,1,1,1),'ks');  % True
plot(numIterationsRange(1),y(1,2,1,1),'ko');  % LS
plot(numIterationsRange(1),y(1,3,1,1),'kd');  % MMSE

% Defining a vector with the positions where the markers will be placed
markers = round(linspace(1,length(numIterationsRange),100));

% Define a string containing the markers for true/LS/MMSE channels
markersStr = ['s','o','d'];

% Go through all K values
for k = 1:length(Krange)
    
    % Go through all channel estimate types
    for est = 1:3
        
        % Select uncorrelated or correlated fading channels
        for chn = 1:2
            
            if chn == 1
                
                str = '--';
                
            elseif chn == 2
                
                str = '-.';
                
            end
           
            % Plotting process
            ax.ColorOrderIndex = k;
            plot(numIterationsRange,y(:,est,chn,k),str,'LineWidth',1)
            ax.ColorOrderIndex = k;
            plot(numIterationsRange(markers),y(markers,est,chn,k),markersStr(2),'LineWidth',1)
            
        end
        
    end
    
end

% Plot settings
xlabel('Average number of rKA iterations ($\bar{T}_{\scriptscriptstyle\mathrm{rKA}}$)')
ylabel('Gap to the average UL RZF SE per UE bound $[\%]$');

ylim([10^(-1) 100])

set(gca,'YScale','log');
set(gca,'XScale','log');

legend('Uncorrelated','Correlated','True channel','LS estimator','MMSE estimator','Location','NorthEast');