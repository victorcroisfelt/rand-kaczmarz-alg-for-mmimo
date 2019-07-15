function [SE_true,SE_LS,SE_MMSE] = functionComputeSE_UL_PARLrKA(M,K,p,tau_c,tau_p,numRealizations,numIterationsRange,H,V)
% Compute the spectral efficiency (SE) for the uplink, relying on the
% use-and-then-forget (UatF) lower bound, defined in (4.14), p. 302 of [1].
% This function do the computions regarding three different conditions of
% channel state information knownledge on the BS side: perfect knowledge,
% LS channel estimates, MMSE channel estimates when considering the
% PARL rKA-based scheme.
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
% Input:
%   - M: number of BS antennas.
%   - K: number of UEs.
%   - p: total uplink transmit power per UE [mW].
%   - tau_c: length of the coherence block [samples].
%   - tau_p: length of the pilot sequences [samples].
%   - numRealizations: number of small-scale fading realizations.
%   - H: M x nbrOfRealizations x K matrix with the true channel values.
%   - V_hyb: M x K x numRealizations x length(numIterationsRange) x 3
%   matrix with the receive combining vectors for the PARL rKA-based scheme
%   using true channel values. The last index is with respect to the 
%   channel estimation method; if 1 == true channel, if 2 == LS channel 
%   estimates, and if 3 == MMSE channel estimates.
%
% Output:
%   - SE_true: K x length(numIterationsRange) vector with the SE values of
%   each kth UE for the PARL rKA-based scheme using true channel values.
%   - SE_LS: K x length(numIterationsRange) vector with the SE values of
%   each kth UE for the PARL rKA-based scheme using LS channel estimates.
%   - SE_MMSE: K x length(numIterationsRange) vector with the SE values of
%   each kth UE for the PARL rKA-based scheme using MMSE channel estimates.
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

% Compute the pre-log factor (considering uplink transmission only)
prelogFactor = (tau_c-tau_p)/(tau_c);

% Signal gains
signal = zeros(K,length(numIterationsRange),3);

% Norms of combining vectors
combiningNorm = zeros(K,length(numIterationsRange),3);

% Sum of interference powers
interf = zeros(K,length(numIterationsRange),3);

% Prepare to save the SEs
SE = zeros(K,length(numIterationsRange),3);

%% Go through all channel realizations
for n = 1:numRealizations
    
    % Extracting channel realizations from all users to the BS
    Hall = reshape(H(:,n,:),[M K]);
    
    % Go through all users
    for k = 1:K
        
        % Go through all iterates
        for int = 1:length(numIterationsRange)
            
            % Go through all different channel estimates
            for est = 1:3
                
                v = V(:,k,n,int,est); % extract the combining vector
                
                % Compute signal and interference + noise terms
                signal(k,int,est) = signal(k,int,est) + (v'*Hall(:,k))/numRealizations;
                combiningNorm(k,int,est) = combiningNorm(k,int,est) + (norm(v).^2)/numRealizations;
                interf(k,int,est) = interf(k,int,est) + p*sum(abs(v'*Hall).^2)/numRealizations;
                
                clear v
                
            end
            
        end
        
    end
    
end

%% Go through all iterates
for int = 1:length(numIterationsRange)
    
    % Go through all different channel estimates
    for est = 1:3
        
        % Compute the SEs
        SE(:,int,est) = prelogFactor*real(log2(1 + ((p*abs(signal(:,int,est)).^2) ./ (interf(:,int,est) - (p*abs(signal(:,int,est)).^2) + combiningNorm(:,int,est)))));
        
    end
    
end

% Prepare the outputs
SE_true = SE(:,:,1);
SE_LS   = SE(:,:,2);
SE_MMSE = SE(:,:,3);

%% Treatitng NaN numbers
SE_true(isnan(SE_true)) = 0;
SE_LS(isnan(SE_LS))     = 0;
SE_MMSE(isnan(SE_MMSE)) = 0;
