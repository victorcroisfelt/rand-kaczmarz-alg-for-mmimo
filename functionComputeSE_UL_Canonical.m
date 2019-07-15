function [SE_true,SE_LS,SE_MMSE] = functionComputeSE_UL_Canonical(M,K,p,tau_c,tau_p,numRealizations,H,Hhat_LS,Hhat_MMSE)
% Compute the spectral efficiency (SE) for the uplink, relying on the
% use-and-then-forget (UatF) lower bound, defined in (4.14), p. 302 of [1].
%
% This function computes only the spectral efficiency of the canonical RZF
% scheme.
%
% This Matlab function is used in the article:
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
% Input:
%   - M: number of BS antennas.
%   - K: number of UEs.
%   - p: total uplink transmit power per UE [mW].
%   - tau_c: length of the coherence block [samples].
%   - tau_p: length of the pilot sequences [samples].
%   - numRealizations: number of small-scale fading realizations.
%   - H: M x nbrOfRealizations x K matrix with the true channel values.
%   - Hhat_LS: M x nbrOfRealizations x K matrix with LS channel estimates.
%   - Hhat_MMSE: M x nbrOfRealizations x K matrix with MMSE channel
%   estimates.
%
% Output:
%   - SE_true: K x 1 vector with the SE values of each kth UE for the
%   canonical RZF scheme using true channel values.
%   - SE_LS: K x 1 vector with the SE values of each kth UE for the
%   canonical RZF scheme using LS channel estimates.
%   - SE_MMSE: K x 1 vector with the SE values of each kth UE for the
%   canonical RZF scheme using MMSE channel estimates.
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
%% Preamble

% Hold the K x K identity matrix
eyeK = eye(K);

% Compute the pre-log factor (considering uplink transmission only)
prelogFactor = (tau_c-tau_p)/(tau_c);

% Signal gains
signal = zeros(K,3);

% Norms of combining vectors
combiningNorm = zeros(K,3);

% Sum of interference powers
interf = zeros(K,3);

%% Go through all channel realizations
for n = 1:numRealizations
    
    % Extract the channel realizations from all users to the BS
    Hall = reshape(H(:,n,:),[M K]);
    Hhatall_LS = reshape(Hhat_LS(:,n,:),[M K]);
    Hhatall_MMSE = reshape(Hhat_MMSE(:,n,:),[M K]);
    
    % Compute maximum-ratio combining (MRC)
    V_MR(:,:,1) = Hall;
    V_MR(:,:,2) = Hhatall_LS;
    V_MR(:,:,3) = Hhatall_MMSE;
    
    % Compute RZF combining
    V_RZF(:,:,1) = p*V_MR(:,:,1)/(p*(V_MR(:,:,1)'*V_MR(:,:,1))+eyeK);
    V_RZF(:,:,2) = p*V_MR(:,:,2)/(p*(V_MR(:,:,2)'*V_MR(:,:,2))+eyeK);
    V_RZF(:,:,3) = p*V_MR(:,:,3)/(p*(V_MR(:,:,3)'*V_MR(:,:,3))+eyeK);
    
    % Go through all users
    for k = 1:K
        
        % Go through all different channel estimates
        for est = 1:3
            
            % Extract the combining vector
            v_RZF = V_RZF(:,k,est);
            
            % Compute signal and interference + noise terms
            signal(k,est) = signal(k,est) + (v_RZF'*Hall(:,k))/numRealizations;
            combiningNorm(k,est) = combiningNorm(k,est) + (norm(v_RZF).^2)/numRealizations;
            interf(k,est) = interf(k,est) + p*sum(abs(v_RZF'*Hall).^2)/numRealizations;
            
            clear v_RZF
            
        end
        
    end
    
end

%% Compute the SEs

% Prepare to save the SE values
SE = zeros(K,3);

% Go through all different channel estimates
for est = 1:3
    
    SE(:,est) = prelogFactor*real(log2(1 + ((p*abs(signal(:,est)).^2) ./ (interf(:,est) - (p*abs(signal(:,est)).^2) + combiningNorm(:,est)))));
    
end

% Assigning outputs
SE_true = SE(:,1);
SE_LS   = SE(:,2);
SE_MMSE = SE(:,3);

