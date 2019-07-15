% This Matlab script can be used to generate the Figure 5 in the article:
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

%% Simulation Parameters

% Range of the number of BS antennas
M = 50:300;

% Compute the number of users to obtain K/M = 0.1, do the same for K/M =
% 0.3 and K/M = 0.5
K01 = .1*M;
K03 = .3*M;
K05 = .5*M;

% Computing the upper bounds, as described by (23b)
TrKA01 = (K01.^3)./(3*M) + (K01.^2)/2 + (4*K01-(3*M.*K01))./(6*M);
TrKA03 = (K03.^3)./(3*M) + (K03.^2)/2 + (4*K03-(3*M.*K03))./(6*M);
TrKA05 = (K05.^3)./(3*M) + (K05.^2)/2 + (4*K05-(3*M.*K05))./(6*M);

%% Plot Figure 5

figure;
hold on; box on; ax = gca;

% Markes: Uncorr. - Error Bound of 10%
plot([136 105 104],[TrKA01(M==136) TrKA03(M==105) TrKA05(M==104)],'ro','LineWidth',1,'MarkerFaceColor','r')

% Markers: Corr. - Error Bound of 10%
plot([138 112 109],[TrKA01(M==138) TrKA03(M==112) TrKA05(M==109)],'b*','LineWidth',1,'MarkerFaceColor','r')

% Plot the upper bound curves
plot(M,TrKA01,'k-', 'LineWidth',1.0)
plot(M,TrKA03,'k--','LineWidth',1.0)
plot(M,TrKA05,'k-.','LineWidth',1.0)

% Trade-off line: Uncorr. - Error Bound of 10%
plot([136 105 104],[TrKA01(M==136) TrKA03(M==105) TrKA05(M==104)],'ro-','LineWidth',1,'MarkerFaceColor','r')

% Trade-off line: Corr. - Error Bound o 10%
plot([138 112 109],[TrKA01(M==138) TrKA03(M==112) TrKA05(M==109)],'b*-','LineWidth',1,'MarkerFaceColor','r')

% Trade-off line: Uncorr. - Error Bound of 1%
plot([239 184 173],[TrKA01(M==239) TrKA03(M==184) TrKA05(M==173)],'ro-','LineWidth',1,'MarkerFaceColor','r')

% Trade-off line: Corr. - Error Bound of 1%
plot([254 189 175],[TrKA01(M==254) TrKA03(M==189) TrKA05(M==175)],'b*-','LineWidth',1,'MarkerFaceColor','r')

% Plot settings
ylabel('Number of rKA iterations ($T_{\scriptscriptstyle\mathrm{rKA}}$)');
xlabel('Number of BS antennas ($M$)');

legend('Uncorrelated','Correlated','$\overline{T}^{\scriptscriptstyle\mathrm{RZF}}_{\scriptscriptstyle\mathrm{rKA}}$ for ${K}/{M}=0.1$','$\overline{T}^{\scriptscriptstyle\mathrm{RZF}}_{\scriptscriptstyle\mathrm{rKA}}$ for ${K}/{M}=0.3$','$\overline{T}^{\scriptscriptstyle\mathrm{RZF}}_{\scriptscriptstyle\mathrm{rKA}}$ for ${K}/{M}=0.5$','Location','SouthEast');

set(gca,'YScale','log');

xlim([min(M) max(M)])



