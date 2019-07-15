% This Matlab script can be used to generate the Figure 1 in the articles:
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
% (github:https://github.com/emilbjornson/massivemimobook).
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

% Standard deviation of the large scale fading (shadowing) [dB]
stdLSF = 4;

% Define the range of the pathloss exponents
alphaRange = [2 4];

% Set the number of setups
numSetups = 1000;

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

%% Simulation

% Prepare to save simulation results
logSampleProb = zeros(K,numSetups,3,2,length(alphaRange));

% Go through all simulations points
for alpha = 1:length(alphaRange)
    
    % Output simulation progress
    disp([num2str(alpha) ' alpha out of ' num2str(length(alphaRange))]);
    
    for s = 1:numSetups
        
        % Output simulation progress
        disp([num2str(s) ' setup out of ' num2str(numSetups)]);
        
        % Compute the setup characteristics
        [Runcorr,Rcorr,channelGaindB] = functionSetup(M,K,correlationFactor,stdLSF,alphaRange(alpha));
        
        % Compute the normalized channel gains, where the normalization is
        % by the noise power
        channelGainOverNoise = channelGaindB - noiseVariancedBm;
        
        % Select uncorrelated or correlated fading channels
        for chn = 1:2
            
            if chn == 1
                
                % Compute the channel estimates
                [H,Hhat_LS,Hhat_MMSE,~] = functionChannelEstimates(M,K,p,Runcorr,channelGainOverNoise,numRealizations);
                
            else
                
                % Compute the channel estimates
                [H,Hhat_LS,Hhat_MMSE,~] = functionChannelEstimates(M,K,p,Rcorr,channelGainOverNoise,numRealizations);
                
            end
            
            % Go through all channel realizations
            for n = 1:numRealizations
                
                % Extract channel estimate realizations from all UEs to the 
                % BS
                Hall(:,:,1) = reshape(H(:,n,:),[M K]);
                Hall(:,:,2) = reshape(Hhat_LS(:,n,:),[M K]);
                Hall(:,:,3) = reshape(Hhat_MMSE(:,n,:),[M K]);
    
                % Define the function to compute the sample probability
                sampleProbability = @(col) ((norm(Hall(:,col,1))^2) + (1/p))/((norm(Hall(:,:,1),'fro')^2) + (K/p));
                sampleProbability_LS = @(col) ((norm(Hall(:,col,2))^2) + (1/p))/((norm(Hall(:,:,2),'fro')^2) + (K/p));
                sampleProbability_MMSE = @(col) ((norm(Hall(:,col,3))^2) + (1/p))/((norm(Hall(:,:,3),'fro')^2) + (K/p));
                
                % Compute the PDFs
                pdfRowsRZF(:,:,1) = arrayfun(sampleProbability,(1:K)');
                pdfRowsRZF(:,:,2) = arrayfun(sampleProbability_LS,(1:K)');
                pdfRowsRZF(:,:,3) = arrayfun(sampleProbability_MMSE,(1:K)');
                
                % Save the pdf average value
                logSampleProb(:,s,1,chn,alpha) = logSampleProb(:,s,1,chn,alpha) + pdfRowsRZF(:,:,1)/numRealizations;
                logSampleProb(:,s,2,chn,alpha) = logSampleProb(:,s,2,chn,alpha) + pdfRowsRZF(:,:,2)/numRealizations;
                logSampleProb(:,s,3,chn,alpha) = logSampleProb(:,s,3,chn,alpha) + pdfRowsRZF(:,:,3)/numRealizations;
                
            end
            
        end
        
    end
    
end

% Reshape data & sort the data
reshap = reshape(logSampleProb,[K*numSetups,3,2,length(alphaRange)]);
sorted = sort(reshap,1);

%% Plot simulation results
figure;
hold on; box on; ax = gca;

% Preparing data
y = linspace(0,1,K*numSetups);

% Legend plots
plot(sorted(1,1,1,1),y(1),'k--','LineWidth',1) % Uncorrelated
plot(sorted(1,1,2,1),y(1),'k-.','LineWidth',1) % Correlated

plot(sorted(1,1,1,1),y(1),'ks','LineWidth',1) % True channel
plot(sorted(1,2,1,1),y(1),'ko','LineWidth',1) % LS
plot(sorted(1,3,1,1),y(1),'kd','LineWidth',1) % MMSE

% Defining a vector with the positions where the markers will be placed
markers = round(linspace(1,K*numSetups,10));

% Define a string containing the markers for true/LS/MMSE channels
markersStr = ['s','o','d'];

% Go through all simulations points
for alpha = 1:length(alphaRange)
    
    % Go through all channel estimate types
    for est = 1:3
        
        % Restarting color indexing of the figure
        ax.ColorOrderIndex = 1;
        
        % Select uncorrelated or correlated fading channels
        for chn = 1:2
            
            if chn == 1
                
                str = '--';
                
            elseif chn == 2
                
                str = '-.';
                
            end
            
            % Plotting process
            plot(sorted(:,est,chn,alpha),y,str,'LineWidth',1) 
            ax.ColorOrderIndex = ax.ColorOrderIndex-1;
            plot(sorted(markers,est,chn,alpha),y(markers),markersStr(est),'LineWidth',1)
            
        end
        
    end
    
end

% Plot settings
ylabel('CDF');
xlabel('Average sample probability ($\bar{P}_{r(t)}$)');

legend('Uncorrelated','Correlated','True channel','LS estimator','MMSE estimator','Location','NorthWest');
set(gca,'XScale','log');

xlim([4*10^(-4) 10^(0)])