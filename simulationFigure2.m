% This Matlab script can be used to generate the Figure 2 in the articles:
%
% Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão,
% "Kaczmarz Precoding and Detection for Massive MIMO Systems", IEEE
% Wireless Communications and Networking Conference, Marrakech, Morroco,
% 2019: 1-6.
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
% [2] M. N. Boroujerdi, S. Haghighatshoar, and G. Caire, "Low-Complexity
% Statistically Robust Precoder/Detector  for Massive MIMO
% Systems", IEEE Transactions on Wireless Communications, vol. 17, no. 10,
% pp. 6516–6530, 2018.
%

% Initialization
close all;
clearvars;

%% Simulation parameters

% Number of BS antennas
M = 100;

% Number of UEs per BS
K = 10;

% Correlation factor of the exponential correlation model
correlationFactor = 0.5;

% Standard deviation of the large-scale fading (shadowing) [dB]
stdLSF = 4;

% Number of rKA iterations
numIterationsRange = [1 100:100:500];

% Set the number of setups
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
meanSE_RZF  = zeros(3,2,numSetups);
meanSE_hyb  = zeros(length(numIterationsRange),3,2,numSetups);
meanSE_neat = zeros(length(numIterationsRange),3,2,numSetups);

% Go through all setups
for s = 1:numSetups
    
    % Setting a timer
    setup = tic;
    
    % Output setups progress
    disp([num2str(s) ' setup out of ' num2str(numSetups)]);
    
    % Compute the setup characteristics
    [Runcorr,Rcorr,channelGaindB] = functionSetup(M,K,correlationFactor,stdLSF);
    
    % Compute the normalized channel gains
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    % Select uncorrelated or correlated fading channels
    for chn = 1:2
        
        if chn == 1
            
            % Compute the channel estimates for LS/MMSE channel estimators
            [H,Hhat_LS,Hhat_MMSE,tau_p] = functionChannelEstimates(M,K,p,Runcorr,channelGainOverNoise,numRealizations);
            
        else
            
            % Compute the channel estimates for LS/MMSE channel estimators
            [H,Hhat_LS,Hhat_MMSE,tau_p] = functionChannelEstimates(M,K,p,Rcorr,channelGainOverNoise,numRealizations);
            
        end
        
        % Compute the UL SE of the RZF canonical scheme
        [SE_RZF,SE_RZF_LS,SE_RZF_MMSE] = functionComputeSE_UL_Canonical(M,K,p,tau_c,tau_p,numRealizations,H,Hhat_LS,Hhat_MMSE);
        
        % Save simulation results
        meanSE_RZF(1,chn,s) = mean(SE_RZF);
        meanSE_RZF(2,chn,s) = mean(SE_RZF_LS);
        meanSE_RZF(3,chn,s) = mean(SE_RZF_MMSE);
        
        % Run PARL rKA-based schemes with and w/o the hybrid initialization 
        %approach
        [V_hyb,V_neat] = functionPARLrKA(M,K,p,numRealizations,numIterationsRange,H,Hhat_LS,Hhat_MMSE,1);
        
        % Compute the UL SEs for the hybrid approach
        [SE_hyb_true,SE_hyb_LS,SE_hyb_MMSE] = functionComputeSE_UL_PARLrKA(M,K,p,tau_c,tau_p,numRealizations,numIterationsRange,H,V_hyb);
        
        % Save simulation results
        meanSE_hyb(:,1,chn,s) = mean(SE_hyb_true,1);
        meanSE_hyb(:,2,chn,s) = mean(SE_hyb_LS,1);
        meanSE_hyb(:,3,chn,s) = mean(SE_hyb_MMSE,1);
        
        % Compute the UL SEs for the neat approach
        [SE_neat_true,SE_neat_LS,SE_neat_MMSE] = functionComputeSE_UL_PARLrKA(M,K,p,tau_c,tau_p,numRealizations,numIterationsRange,H,V_neat);
        
        % Save simulation results
        meanSE_neat(:,1,chn,s) = mean(SE_neat_true,1);
        meanSE_neat(:,2,chn,s) = mean(SE_neat_LS,1);
        meanSE_neat(:,3,chn,s) = mean(SE_neat_MMSE,1);
        
    end
    
    toc(setup)
    
end

%% Plot simulation results

% Extracting data
SE_RZF  = mean(meanSE_RZF,3);
SE_hyb  = mean(meanSE_hyb,4);
SE_neat = mean(meanSE_neat,4);

% Set the markers' places
markers = round(linspace(1,length(numIterationsRange),10));

figure;
subplot(1,3,1) % True channel
hold on; box on; ax = gca;

% Legend plots
plot(numIterationsRange(1),SE_neat(1,1,1),'--k','LineWidth',1)
plot(numIterationsRange(1),SE_hyb(1,2,1),'-.k','LineWidth',1) 
plot(numIterationsRange(1),SE_neat(1,1,1),'*k','LineWidth',1) 
plot(numIterationsRange(1),SE_hyb(1,1,1),'pk','LineWidth',1) 

% Select uncorrelated or correlated fading channels
for chn = 1:2
    
    if chn == 1
        
        str = '--';
        
    elseif chn == 2
        
        str = '-.';
        
    end
    
    ax.ColorOrderIndex = 1;
    
    % Plotting process
    plot(numIterationsRange,SE_neat(:,1,chn),'--','LineWidth',1)
    plot(numIterationsRange,SE_hyb(:,1,chn),'--','LineWidth',1) 
   
    ax.ColorOrderIndex = 1;
    
    plot(numIterationsRange(markers),SE_neat(markers,1,chn),'*','LineWidth',1) 
    plot(numIterationsRange(markers),SE_hyb(markers,1,chn),'p','LineWidth',1)
    
end

% Bound plots
plot(numIterationsRange,SE_RZF(1,1)*ones(length(numIterationsRange),1),'k--','LineWidth',1.5)
plot(numIterationsRange,SE_RZF(1,2)*ones(length(numIterationsRange),1),'k-.','LineWidth',1.5)

% Plot settings
xticks([0 100 200 300 400 500])

ylim([0 8])
xlim([0 500])

ylabel('Average UL SE per UE [bit/s/Hz/user]');
xlabel('Number of rKA iterations ($T_{\scriptscriptstyle\mathrm{rKA}}$)');

legend('Uncorrelated','Correlated','PARL rKA-based w/o init. [6]','PARL rKA-based w/ hyb. init.','Location','SouthEast');

title('True channel')

subplot(1,3,2) % LS estimator
hold on; box on; ax = gca;

% Legend plots
plot(numIterationsRange(1),SE_neat(1,2,1),'--k','LineWidth',1)
plot(numIterationsRange(1),SE_hyb(1,2,2),'-.k','LineWidth',1) 
plot(numIterationsRange(1),SE_neat(1,2,1),'*k','LineWidth',1) 
plot(numIterationsRange(1),SE_hyb(1,2,1),'pk','LineWidth',1) 

% Select uncorrelated or correlated fading channels
for chn = 1:2
    
    if chn == 1
        
        str = '--';
        
    elseif chn == 2
        
        str = '-.';
        
    end
    
    ax.ColorOrderIndex = 1;
    
    % Plotting process
    plot(numIterationsRange,SE_neat(:,2,chn),'--','LineWidth',1)
    plot(numIterationsRange,SE_hyb(:,2,chn),'--','LineWidth',1) 
   
    ax.ColorOrderIndex = 1;
    
    plot(numIterationsRange(markers),SE_neat(markers,2,chn),'*','LineWidth',1) 
    plot(numIterationsRange(markers),SE_hyb(markers,2,chn),'p','LineWidth',1)
    
end

% Bound plots
plot(numIterationsRange,SE_RZF(2,1)*ones(length(numIterationsRange),1),'k--','LineWidth',1.5)
plot(numIterationsRange,SE_RZF(2,2)*ones(length(numIterationsRange),1),'k-.','LineWidth',1.5)

% Plot settings
xticks([0 100 200 300 400 500])

ylim([0 8])
xlim([0 500])

xlabel('Number of rKA iterations ($T_{\scriptscriptstyle\mathrm{rKA}}$)')

legend('Uncorrelated','Correlated','PARL rKA-based w/o init. [6]','PARL rKA-based w/ hyb. init.','Location','SouthEast');

title('LS estimator')

subplot(1,3,3) % MMSE estimator
hold on; box on; ax = gca;

% Legend plots
plot(numIterationsRange(1),SE_neat(1,3,1),'--k','LineWidth',1)
plot(numIterationsRange(1),SE_hyb(1,3,2),'-.k','LineWidth',1) 
plot(numIterationsRange(1),SE_neat(1,3,1),'*k','LineWidth',1) 
plot(numIterationsRange(1),SE_hyb(1,3,1),'pk','LineWidth',1) 

% Select uncorrelated or correlated fading channels
for chn = 1:2
    
    if chn == 1
        
        str = '--';
        
    elseif chn == 2
        
        str = '-.';
        
    end
    
    ax.ColorOrderIndex = 1;
    
    % Plotting process
    plot(numIterationsRange,SE_neat(:,3,chn),'--','LineWidth',1)
    plot(numIterationsRange,SE_hyb(:,3,chn),'--','LineWidth',1) 
   
    ax.ColorOrderIndex = 1;
    
    plot(numIterationsRange(markers),SE_neat(markers,3,chn),'*','LineWidth',1) 
    plot(numIterationsRange(markers),SE_hyb(markers,3,chn),'p','LineWidth',1)
    
end

% Bound plots
plot(numIterationsRange,SE_RZF(3,1)*ones(length(numIterationsRange),1),'k--','LineWidth',1.5)
plot(numIterationsRange,SE_RZF(3,2)*ones(length(numIterationsRange),1),'k-.','LineWidth',1.5)

% Plot settings
xticks([0 100 200 300 400 500])

ylim([0 8])
xlim([0 500])

xlabel('Number of rKA iterations ($T_{\scriptscriptstyle\mathrm{rKA}}$)');

legend('Uncorrelated','Correlated','PARL rKA-based w/o init. [6]','PARL rKA-based w/ hyb. init.','Location','SouthEast');

title('MMSE estimator')