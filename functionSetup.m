function [Runcorr,Rcorr,channelGaindB] = functionSetup(M,K,correlationFactor,stdLSF,alpha)
% Generate a Massive MIMO system with a square-network layout comprised of
% a BS equipped an array of M antennas and serving K users equipments UEs.
% The code uniformly distributes the UEs within the cell area. The
% covariance matrices for uncorrelated and correlated fading channels are
% also generated as Runcorr and Rcorr, respectively. Correlated fading
% channels are being modeled using an exponential correlation model with
% large-scale fading variations over the array. The degree of antenna
% correlation is given by the correlationFactor and the stardard deviation
% of the shadowing by stdLSF. It is important to emphasize that, for the
% uncorrelated case, the large-scale fading is constant to all the elements
% of the antenna array. Moreover, the pathloss exponent can be set manually
% by changing the value of alpha. It is also possible to generate several
% realizations of the correlated fading matrices, under the same physical
% setup (distribution of the UEs in the coverage area), by given a vector
% of values for the correlationFactor or/and stdLSF.
%
% This Matlab function is used in the articles:
%
% Victor Croisfelt Rodrigues, José Carlos Marinello Filho, and Taufik Abrão,
% "Randomized Kaczmarz Algorithm for Massive MIMO Systems with Channel
% Estimation and Spatial Correlation", to appear, 2019.
%
% Download article:
%                   https://arxiv.org/abs/1904.04376
%
% Input:
%   - M: number of BS antennas.
%   - K: number of UEs.
%   - correlationFactor: physical correlation degree between the antenna
%     elements.
%   - stdLSF: large-scale fading (shadowing) standard deviation.
%   - alpha: pathloss exponent, default value is 3.76.
%
% Output:
%   - Runcorr: an M x M x K matrix with K x (M x M) covariance matrices of
%     the K UEs. Covariance matrices for uncorrelated channels.
%   - Rcorr: an M x M x K x matrix with K x (M x M) covariance matrices of
%     the K UEs. Covariance matrices for correlated channels. It is
%     possible to generate several statisticalRealizations by vectors of
%     correlationFactor and stdLSF.
%   - channelGaindB: a K x 1 vector with the power of the K UEs connected
%     to the BS.
%
% This is version 4.0 (Last edited: 2019-04-07)
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
%% Model parameters

% Set the length in meters of the total square area
squareLength = 250;

% Pathloss exponent
if nargin < 5
    
    alpha = 3.76; % default
    
end

% Average channel gain in dB at a reference distance of 1 meter. Note that
% -35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;

% Minimum distance between BSs and UEs
minDistance = 35;

% Deploy BSs on the grid
BSpositions = squareLength/2 + 1i*squareLength/2;

% Prepare to put out UEs in the cell
UEpositions = zeros(K,1);

% Check the channel condition
if (length(correlationFactor) == 1 && length(stdLSF) == 1)
    
    % #1 Scenario: correlationFactor and stdLSF are scalars
    correlationScenarios = 1;
    aux = 1;
    
else
    
    if (length(correlationFactor) > 1 && length(stdLSF) == 1)
        
        % #2 Scenario: correlationFactor is a vector and stdLSF is a scalar
        correlationScenarios = length(correlationFactor);
        aux = 2;
        
    elseif (length(correlationFactor) == 1 && length(stdLSF) > 1)
        
        % #3 Scenario: correlationFactor is a scalar and stdLSF is a vector
        correlationScenarios = length(stdLSF);
        aux = 3;
        
    end
    
end

% Prepare to store spatial correlation matrices
Runcorr = zeros(M,M,K);
Rcorr = zeros(M,M,K,correlationScenarios);

%% Put out K UEs in the cell, uniformly at random.
% The procedure is iterative since UEs that do not satisfy the minimum
% distance are replaced with new positions until contraint agreement.

% Starting loop variable
perBS = 0;

while perBS < K
    
    % Put out new UEs
    UEremaining = K - perBS;
    posX = rand(UEremaining,1)*squareLength - squareLength/2;
    posY = rand(UEremaining,1)*squareLength - squareLength/2;
    posXY = posX + 1i*posY;
    
    % Keep those that satisfy the minimum distance
    posXY = posXY(abs(posXY) >= minDistance);
    
    % Store new UEs
    UEpositions(perBS+1:perBS+length(posXY)) = posXY + BSpositions;
    perBS = perBS + length(posXY);
    
end

% Compute distances between UEs and the BS
distancesBSj = abs(UEpositions - BSpositions);

% Compute angles between UEs and the BS
angleBSj = angle(UEpositions - BSpositions);

% Compute distant-dependent path gains (in dB)
channelGaindB = constantTerm - alpha*10*log10(distancesBSj);

% Generate the spatial covariance matrices of all UEs
for k = 1:K
    
    % Check aux variable:
    if aux == 1
        
        % Generate the covariance matrix using the uncorrelated fading
        % model
        Runcorr(:,:,k) = 10.^(stdLSF*randn/10)*eye(M);
        
        % Generate realizations of large-scale fading variations over the
        % array for correlated Rayleigh fading channels
        fadingOverArray = randn(M,1);
        
        % Store in a diag matrix the square root of the power
        largeScaleFadingD = diag(10.^(stdLSF*fadingOverArray/20));
        
        % Generate the covariance matrix using the exponential correlation
        % model, including large-scale fading variations over the array
        Rcorr(:,:,k) = largeScaleFadingD*toeplitz((correlationFactor*exp(1i*angleBSj(k))).^(0:M-1))*largeScaleFadingD;
        
    elseif aux == 2
        
        % Go through correlation scenarios
        for r = 1:correlationScenarios
            
            % Generate realizations of large-scale fading variations over
            % the array for correlated Rayleigh fading channels
            fadingOverArray = randn(M,1);
            
            % Store in a diag matrix with the square root of the power
            largeScaleFadingD = diag(10.^(stdLSF*fadingOverArray/20));
            
            % Generate the covariance matrix using the exponential
            % correlation model, including large-scale fading variations
            % over the array
            Rcorr(:,:,k) = largeScaleFadingD*toeplitz((correlationFactor(r)*exp(1i*angleBSj(k))).^(0:M-1))*largeScaleFadingD;
            
        end
        
    elseif aux == 3
        
        % Go through correlation scenarios
        for r = 1:correlationScenarios
            
            % Generate realizations of large-scale fading variations over
            % the array for correlated Rayleigh fading channels
            fadingOverArray = randn(M,1);
            
            % Store in a diag matrix with the square root of the power
            largeScaleFadingD = diag(10.^(stdLSF(r)*fadingOverArray/20));
            
            % Generate the covariance matrix using the exponential
            % correlation model, including large-scale fading variations
            % over the array
            Rcorr(:,:,k,l) = largeScaleFadingD*toeplitz((correlationFactor*exp(1i*angleBSj(k))).^(0:M-1))*largeScaleFadingD;
            
        end
        
    end
    
end

end