function [H,Hhat_LS,Hhat_MMSE,tau_p] = functionChannelEstimates(M,K,p,R,channelGaindB,numRealizations)
% Compute the channel estimates by adopting both the LS and MMSE channel
% estimation approaches. Here, the pilot training phase of a Massive MIMO
% network comprised of a BS equipped with M antennas and serving K UEs is
% considered. Each UE has a uplink transmit power of p. To compute their
% channel estimates, it is necessary the channels' covariance matrices, Rs,
% and channel power characterististics (channelGaindB). numRealizations
% represents the realizations of small-scale fading assumed in a
% Monte-Carlo simulation context.
%
% This Matlab function is used in the articles:
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
%Input:
%   - M: number of BS antennas.
%   - K: number of UEs.
%   - p: total uplink transmit power per UE [mW].
%   - R: an M x M x K matrix with K x (M x M) covariance matrices of the K
%     UEs.
%   - channelGaindB: a K x 1 vector with the power of the K UEs connected
%     to the BS.
%   - numRealizations: number of small-scale fading realizations.
%
%Output:
%   - H: M x nbrOfRealizations x K matrix with the true channel values.
%   - Hhat_LS: M x nbrOfRealizations x K matrix with the channel estimates
%   when adopting the LS channel estimator.
%   - Hhat_MMSE: M x nbrOfRealizations x K matrix with the channel
%   estimates when adopting the MMSE channel estimator.
%   - tau_p: length of the pilot training phase.
%
% This is version 5.0 (Last edited: 2019-15-07)
%
% License: This code is licensed under the MIT license. If you in any way
% use this code for research that results in publications, please reference
% our original article as shown above.
%
% References:
% [1] Emil Bjornson, Jakob Hoydis, and Luca Sanguinetti (2017), "Massive
% MIMO Networks: Spectral, Energy, and Hardware Efficiency", Foundations
% and Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI:
% 10.1561/2000000093
%
%% Channel realizations

% Generates small-scale fading realizations
H = sqrt(0.5)*(randn(M,numRealizations,K)+1i*randn(M,numRealizations,K));

% Average large-scale fading (channel gain in [W])
betas = zeros(K,1);

% Go through all UEs
for k = 1:K
    
    % Extract the channel gain
    betas(k) = 10^(channelGaindB(k)/10);
    
    % Apply the channel gain to its respective and generated covariance
    % matrix
    R(:,:,k) = betas(k)*R(:,:,k);
    
    % Apply the channel structure containing all the large-scale effects
    % (pathloss, shadowing and spatial correlation - or not) to the
    % small-scale realizations
    Rsqrt = sqrtm(R(:,:,k)); %the square root means magnitude
    H(:,:,k) = Rsqrt*H(:,:,k);
    
end

%% Estimation Parameters

% Hold the M x M identity matrix
eyeM = eye(M);

% Length of the pilot sequences [samples]
tau_p = K;

% Prepare to store the different channel estimates
Hhat_LS = zeros(M,numRealizations,K);
Hhat_MMSE = zeros(M,numRealizations,K);

% Generate realizations of normalized noise (AWGN)
Np = sqrt(0.5)*(randn(M,numRealizations,K)+1i*randn(M,numRealizations,K));

%% Compute processed pilot signal for all UEs
yp = sqrt(p)*tau_p*H + sqrt(tau_p)*Np;

% Go through all UEs
for k = 1:K
    
    % Compute LS estimate of the channel between the BS and UE k
    A_LS = 1/(sqrt(p)*tau_p);
    Hhat_LS(:,:,k) = A_LS*yp(:,:,k);
    
    % Compute the matrix that is inverted in the MMSE estimator
    PsiInv = (p*tau_p*R(:,:,k) + eyeM);
    
    % Compute the MMSE estimate of the channel between the BS and UE k
    RPsi = R(:,:,k)/PsiInv;
    Hhat_MMSE(:,:,k) = sqrt(p)*RPsi*yp(:,:,k);
    
end

end
