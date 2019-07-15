% This Matlab script can be used to generate the Figure 4 in the article:
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

%% Define simulation scenario
% Select which figure to generate:
%   Set simulation = 1 to generate Figure 4a
%   Set simulation = 2 to generate Figure 4b
%
simulation = 1;

%% Simulation parameters

if simulation == 1
    
    % Correlation factor of the exponential correlation model
    correlationFactor = [0 0.25 0.5 0.75 1];
    
    % Standard deviation of the large-scale fading (shadowing) [dB]
    stdLSF = 0;
    
    % Number of points along the abcissa axis in final figure
    numPoints = length(correlationFactor);
    
elseif simulation == 2
    
    % Correlation factor of the exponential correlation model
    correlationFactor = 0;
    
    % Standard deviation of the large-scale fading (shadowing) [dB]
    stdLSF = [0 2 4 6 8];
    
    % Number of points along the abcissa axis in final figure
    numPoints = length(stdLSF);
    
end

% Number of BS antennas
M = 100;

% Number of UEs
K = 10;

% Number of rKA iterations
numIterationsRange = [1 100:100:1000];

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
meanSE_RZF  = zeros(3,numPoints,numSetups);
meanSE_hyb  = zeros(length(numIterationsRange),3,numPoints,numSetups);

% Go through all setups
for s = 1:numSetups
    
    % Output setups progress
    disp([num2str(s) ' setup out of ' num2str(numSetups)]);
    
    % Compute the setup characteristics
    [~,Rcorr,channelGaindB] = functionSetup(M,K,correlationFactor,stdLSF);
    
    % Compute the normalized channel gains
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    % Go through all points
    for pts = 1:numPoints
        
        % Extract the current correlation scenario
        RcorrMatrix = Rcorr(:,:,:,pts);
        
        % Compute the channel estimates for LS/MMSE channel estimators
        [H,Hhat_LS,Hhat_MMSE,tau_p] = functionChannelEstimates(M,K,p,RcorrMatrix,channelGainOverNoise,numRealizations);
        
        % Compute the UL SE of the RZF canonical scheme
        [SE_RZF,SE_RZF_LS,SE_RZF_MMSE] = functionComputeSE_UL_Canonical(M,K,p,tau_c,tau_p,numRealizations,H,Hhat_LS,Hhat_MMSE);

        % Save simulation results
        meanSE_RZF(1,pts,s) = mean(SE_RZF);
        meanSE_RZF(2,pts,s) = mean(SE_RZF_LS);
        meanSE_RZF(3,pts,s) = mean(SE_RZF_MMSE);
        
        % Run PARL rKA-based schemes with the hybrid initialization approach
        [V_hyb] = functionPARLrKA(M,K,p,numRealizations,numIterationsRange,H,Hhat_LS,Hhat_MMSE);

        % Compute the UL SEs for the hybrid approach
        [SE_hyb_true,SE_hyb_LS,SE_hyb_MMSE] = functionComputeSE_UL_PARLrKA(M,K,p,tau_c,tau_p,numRealizations,numIterationsRange,H,V_hyb);

        % Save simulation results
        meanSE_hyb(:,1,pts,s) = mean(SE_hyb_true,1);
        meanSE_hyb(:,2,pts,s) = mean(SE_hyb_LS,1);
        meanSE_hyb(:,3,pts,s) = mean(SE_hyb_MMSE,1);
        
    end
    
end

%% Plot simulation results

if simulation == 1
    
    figure;
    hold on; box on; grid on;
    
    surf(correlationFactor,numIterationsRange,100*(abs(mean(meanSE_hyb(1,:,:,:),4)-mean(meanSE_RZF,3)))./(mean(meanSE_RZF,3)),'FaceAlpha',0.2)
    surf(correlationFactor,numIterationsRange,100*(abs(mean(meanSE_hyb(2,:,:,:),4)-mean(meanSE_RZF,3)))./(mean(meanSE_RZF,3)),'FaceAlpha',0.2)
    surf(correlationFactor,numIterationsRange,100*(abs(mean(meanSE_hyb(3,:,:,:),4)-mean(meanSE_RZF,3)))./(mean(meanSE_RZF,3)),'FaceAlpha',0.2)
    
    colormap(autumn);
    shading interp
    
    set(gca,'ZScale','log');
    
    xlabel('$r$')
    ylabel('$\bar{T}_{\scriptscriptstyle\mathrm{rKA}}$');
    zlabel('Gap to the average UL RZF SE per UE bound [\%]');
    
    ylim([0 1000])
    
    view(135,10)
    
elseif simulation == 2
    
    figure;
    hold on; box on; grid on;
    
    surf(stdLSF,numIterationsRange,100*(abs(mean(meanSE_hyb(1,:,:,:),4)-mean(meanSE_RZF,3)))./(mean(meanSE_RZF,3)),'FaceAlpha',0.2)
    surf(stdLSF,numIterationsRange,100*(abs(mean(meanSE_hyb(2,:,:,:),4)-mean(meanSE_RZF,3)))./(mean(meanSE_RZF,3)),'FaceAlpha',0.2)
    surf(stdLSF,numIterationsRange,100*(abs(mean(meanSE_hyb(3,:,:,:),4)-mean(meanSE_RZF,3)))./(mean(meanSE_RZF,3)),'FaceAlpha',0.2)
    
    colormap(autumn);
    shading interp
    
    set(gca,'ZScale','log');
    
    xlabel('$\sigma$')
    ylabel('$\bar{T}_{\scriptscriptstyle\mathrm{rKA}}$');
    zlabel('Gap to the average UL RZF SE per UE bound [\%]');
    
    ylim([0 1000])
    zlim([10^(-3) 10^(2)])
    
    view(135,10)
    
end