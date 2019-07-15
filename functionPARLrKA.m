function [V_hyb,V_neat] = functionPARLrKA(M,K,p,numRealizations,numIterationsRange,H,Hhat_LS,Hhat_MMSE,neat)
% Perform the proposed rKA-based algorithm. This code parallelly computes
% PARL rKA-based RZF scheme for true, LS, and MMSE channel estimates, in
% order to reduce the simulation time. To run the neat randomized Kaczmarz
% algorithm only needs to set neat = 1
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
%   - numRealizations: number of small-scale fading realizations.
%   - numIterationsRange: range of the number of iterations to run the rKA.
%   Recall that this approach is iterative.
%   - H: M x nbrOfRealizations x K matrix with the true channel values.
%   - Hhat_LS: M x nbrOfRealizations x K matrix with LS channel estimates.
%   - Hhat_MMSE: M x nbrOfRealizations x K matrix with MMSE channel
%   estimates.
%   - neat: if neat is ommitted, only our proposed scheme will be computed;
%   otherwise, for near == 1, we also run the PARL rKA-based scheme w/o any
%   initialization procedure.
%
% Output:
%   - V_hyb: M x K x numRealizations x length(numIterationsRange) x 3
%   matrix with the receive combining vectors for the proposed PARL
%   rKA-based scheme using true channel values. The last index is with
%   respect to the channel estimation method; if 1 == true channel, if 2 ==
%   LS channel estimates, and if 3 == MMSE channel estimates.
%   - V_neat: M x K x numRealizations x length(numIterationsRange)
%   matrix with the receive combining vectors for the PARL rKA-based scheme
%   in [2] using true channel values. The last index is with respect to the
%   channel estimation method; if 1 == true channel, if 2 == LS channel 
%   estimates, and if 3 == MMSE channel estimates. This output depends on
%   the value of the neat input.
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
% Statistically Robust Precoder/Detector Computation for Massive MIMO
% Systems", IEEE Transactions on Wireless Communications, vol. 17, no. 10,
% pp. 6516ï¿½6530, 2018.
%
%% Preamble

% Extract the max number of iterations
maxNumIterations = max(numIterationsRange);

% Prepare to store the simulation results
V_hyb = zeros(M,K,numRealizations,length(numIterationsRange),3);

if nargin > 8 && neat == 1
    
    V_neat = zeros(M,K,numRealizations,length(numIterationsRange),3);
    
end

%% Go through all channel realizations
for n = 1:numRealizations
    
    % Extract channel estimate realizations from all UEs to the BS
    Hall(:,:,1) = reshape(H(:,n,:),[M K]);
    Hall(:,:,2) = reshape(Hhat_LS(:,n,:),[M K]);
    Hall(:,:,3) = reshape(Hhat_MMSE(:,n,:),[M K]);
    
    % Define the function to compute the sample probability
    sampleProbability = @(col) ((norm(Hall(:,col,1))^2) + (1/p))/((norm(Hall(:,:,1),'fro')^2) + (K/p));
    sampleProbability_LS = @(col) ((norm(Hall(:,col,2))^2) + (1/p))/((norm(Hall(:,:,2),'fro')^2) + (K/p));
    sampleProbability_MMSE = @(col) ((norm(Hall(:,col,3))^2) + (1/p))/((norm(Hall(:,:,3),'fro')^2) + (K/p));
    
    % Generate the random rows which will be acessed by the PARL rKA-based
    % RZF scheme using the sample probability given above
    randRows(:,1) = functionRandp(arrayfun(sampleProbability,(1:K)'),maxNumIterations,1);
    randRows(:,2) = functionRandp(arrayfun(sampleProbability_LS,(1:K)'),maxNumIterations,1);
    randRows(:,3) = functionRandp(arrayfun(sampleProbability_MMSE,(1:K)'),maxNumIterations,1);
    
    % Go through all users
    for k = 1:K
        
        % Define the canonical basis as the signal input
        e = zeros(K,1); e(k) = 1;
        
        % Prepare to store the rKA states
        u = zeros(M,maxNumIterations,3);
        z = zeros(K,maxNumIterations,3);
        
        if nargin > 8 && neat == 1
            
            u_neat = zeros(M,maxNumIterations,3);
            z_neat = zeros(K,maxNumIterations,3);
            
        end
        
        % Go through all iterations values
        for int = 1:maxNumIterations
            
            % Go through all different channel estimates
            for est = 1:3
                
                % Proposed PARL rKA-based RZF scheme
                
                % Hybrid initialization
                if int == 1, rowPick = k; else, rowPick = randRows(int,est); end
                
                % Pick the vector related to the chosen row and store it
                q_H = Hall(:,rowPick,est)'; q = q_H';
                
                % Compute the residual
                eta = (e(rowPick) - (q'*u(:,int,est)) - (z(rowPick,int,est)/p))/((norm(q)^2) + (1/p));
                
                % Update schedule
                u(:,int+1,est) = u(:,int,est) + eta*q;
                z(:,int+1,est) = z(:,int,est);
                z(rowPick,int+1,est) = z(rowPick,int+1,est) + eta;
                
                clear rowPick q_H q eta
                
                if nargin > 8 && neat == 1
                    
                    % PARL rKA-based RZF scheme w/o initialization [2]
                    rowPick = randRows(int,est);
                    
                    % Pick the vector related to the chosen row and store it
                    q_H = Hall(:,rowPick,est)'; q = q_H';
                    
                    % Compute the residual
                    eta = (e(rowPick) - (q'*u_neat(:,int,est)) - (z_neat(rowPick,int,est)/p))/((norm(q)^2) + (1/p));
                    
                    % Update schedule
                    u_neat(:,int+1,est) = u_neat(:,int,est) + eta*q;
                    z_neat(:,int+1,est) = z_neat(:,int,est);
                    z_neat(rowPick,int+1,est) = z_neat(rowPick,int+1,est) + eta;
                    
                    clear rowPick q_H q eta
                    
                end
                
            end
            
            % Check if the number of iterations reaches a desired number
            if sum(int == numIterationsRange)
                
                % Hold the position value
                index = find(int == (numIterationsRange));
                
                % Store the results given the index value
                V_hyb(:,k,n,index,1) = Hall(:,:,1)*z(:,(int+1),1);
                V_hyb(:,k,n,index,2) = Hall(:,:,2)*z(:,(int+1),2);
                V_hyb(:,k,n,index,3) = Hall(:,:,3)*z(:,(int+1),3);
                
                if nargin > 8 && neat == 1
                    
                    % Store the results given the index value
                    V_neat(:,k,n,index,1) = Hall(:,:,1)*z_neat(:,(int+1),1);
                    V_neat(:,k,n,index,2) = Hall(:,:,2)*z_neat(:,(int+1),2);
                    V_neat(:,k,n,index,3) = Hall(:,:,3)*z_neat(:,(int+1),3);
                    
                end
                
            end
            
        end
        
    end
    
end
