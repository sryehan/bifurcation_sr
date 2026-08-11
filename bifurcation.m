%% Bifurcation Analysis of the Companion Hamiltonian System
% Physical ODE (fully-integrated, traveling-wave form):
%   c*a0*W' = (a1-a4)*W - (c*a2/2)*W^2 - (c*a3/4)*W^4   =: g(W)
%
% Companion planar Hamiltonian system (Li Jibin / dynamical-system method):
%   W' = y
%   y' = g(W)/(c*a0) = k1*W - k2*W^2 - k3*W^4
%
%   k1 = (a1-a4)/(c*a0)
%   k2 = a2/(2*a0)
%   k3 = a3/(4*a0)
%
% Conserved Hamiltonian:
%   H(W,y) = 0.5*y^2 - Phi(W) = h
%   Phi(W) = (k1/2)*W^2 - (k2/3)*W^3 - (k3/5)*W^5
%
% Equilibria: y=0, W*(k1 - k2*W - k3*W^3) = 0
%   -> W = 0  and roots of  k3*W^3 + k2*W - k1 = 0   (depressed cubic)
%   -> up to 4 real equilibria total (richer than the C1,C2 case since
%      Phi(W) is NOT purely odd -- the k2 term breaks the symmetry,
%      which is exactly what separates "bell" from "antibell").
%
% Orbit classification (read off the contour plot):
%   - Homoclinic loop around a SADDLE  -> BELL / ANTIBELL soliton
%   - Heteroclinic line connecting two SADDLES -> KINK / ANTIKINK
%   - Closed orbits around a CENTER    -> periodic (cnoidal-type) wave family
%   - Orbit tangent to / passing through W where the leading coefficient
%     of the *original un-integrated* ODE vanishes -> singular (peakon/cuspon)
%     [only appears once a0 is un-frozen -- see Route 2 discussion]

clear; clc; close all;

%% ---- USER-EDITABLE PHYSICAL PARAMETERS ----
% Provide (a0, c, a1, a4, a2, a3) directly; k1,k2,k3 are derived.
% Six representative regimes are pre-loaded below to sweep "richness".
% Edit / add rows as needed: [a0, c, a1, a4, a2, a3]
paramSets = [ ...
    1,  1,  1,  0,   6,   -4;   % Regime 1: richest -> 4 real equilibria
    1,  1,  1,  0,   0,    4;   % Regime 2: pure quartic restoring (a2=0), symmetric
    1,  1,  1,  0,  -6,   -4;   % Regime 3: mirror of Regime 1 (sign of k2 flipped)
    1,  1, -1,  0,   6,   -4;   % Regime 4: k1<0
    1,  1,  1,  0,   6,    4;   % Regime 5: k3>0 (different well orientation)
    1,  1,  0,  0,   6,   -4];  % Regime 6: k1=0 (no linear restoring term)

nSets = size(paramSets,1);

% ---- Loop over regimes, each in a separate figure ----
for idx = 1:nSets
    a0 = paramSets(idx,1); c  = paramSets(idx,2);
    a1 = paramSets(idx,3); a4 = paramSets(idx,4);
    a2 = paramSets(idx,5); a3 = paramSets(idx,6);

    k1 = (a1-a4)/(c*a0);
    k2 = a2/(2*a0);
    k3 = a3/(4*a0);

    % ---- Equilibria: W=0 and roots of k3*W^3 + k2*W - k1 = 0 ----
    if abs(k3) > 1e-12
        cubicCoeffs = [k3, 0, k2, -k1];   % k3 W^3 + 0*W^2 + k2 W - k1
        Wroots = roots(cubicCoeffs);
    else
        Wroots = []; % degenerate: linear/quadratic case, handle separately
        if abs(k2) > 1e-12
            Wroots = roots([-k2, k1]); % -k2 W + k1 = 0  -> from k2 W = k1
        end
    end
    Wroots = Wroots(abs(imag(Wroots)) < 1e-8);
    Wroots = real(Wroots);
    Weq = unique([0; Wroots(:)]);

    % ---- Classify each equilibrium via Jacobian eigenvalues ----
    % J = [0 1; g'(W) 0],  g'(W) = k1 - 2*k2*W - 4*k3*W^3
    fprintf('\n--- Regime %d: a0=%.2g c=%.2g a1=%.2g a4=%.2g a2=%.2g a3=%.2g ---\n',...
        idx,a0,c,a1,a4,a2,a3);
    fprintf('   k1=%.3f  k2=%.3f  k3=%.3f\n',k1,k2,k3);
    types = strings(length(Weq),1);
    for j = 1:length(Weq)
        W0 = Weq(j);
        gprime = k1 - 2*k2*W0 - 4*k3*W0^3;
        if gprime > 1e-10
            types(j) = "SADDLE";
        elseif gprime < -1e-10
            types(j) = "CENTER";
        else
            types(j) = "DEGENERATE";
        end
        fprintf('   Eq at W=%7.4f (y=0)  ->  %s  (g''=%.4f)\n', W0, types(j), gprime);
    end

    % ---- Build Hamiltonian grid ----
    Wspan = max(2, 1.4*max(abs(Weq))+0.5);
    Wg = linspace(-Wspan, Wspan, 500);
    yg = linspace(-Wspan, Wspan, 500);
    [WW,YY] = meshgrid(Wg,yg);
    Phi = (k1/2)*WW.^2 - (k2/3)*WW.^3 - (k3/5)*WW.^5;
    H = 0.5*YY.^2 - Phi;

    % ---- Plot in a dedicated figure ----
    figure('Position', [100 + 50*(idx-1), 100 + 50*(idx-1), 700, 650]);
    contour(WW,YY,H,60,'LineWidth',0.6); hold on;
    colormap(gca,turbo);
    plot(Weq, zeros(size(Weq)), 'k.', 'MarkerSize', 18);
    for j = 1:length(Weq)
        text(Weq(j), 0.06*Wspan, sprintf('%s',types(j)), ...
            'FontSize',8,'HorizontalAlignment','center');
    end
    xlabel('W'); ylabel('y = W'''); axis square; grid on;
    title(sprintf('Regime %d:  k_1=%.2g, k_2=%.2g, k_3=%.2g',idx,k1,k2,k3),...
        'FontSize',10);
    drawnow;
end

%% ---- Save richest regime separately at high resolution, with annotation ----
figure('Position',[200 200 700 650]);
a0=1; c=1; a1=1; a4=0; a2=6; a3=-4;
k1 = (a1-a4)/(c*a0); k2 = a2/(2*a0); k3 = a3/(4*a0);
Wg = linspace(-3,3,700); yg = linspace(-3,3,700);
[WW,YY] = meshgrid(Wg,yg);
Phi = (k1/2)*WW.^2 - (k2/3)*WW.^3 - (k3/5)*WW.^5;
H = 0.5*YY.^2 - Phi;
contourf(WW,YY,H,80,'LineColor','none'); hold on;
contour(WW,YY,H,30,'k','LineWidth',0.4);
colormap(turbo); colorbar;
Wroots = roots([k3,0,k2,-k1]);
Wroots = real(Wroots(abs(imag(Wroots))<1e-8));
Weq = unique([0;Wroots]);
plot(Weq,zeros(size(Weq)),'wo','MarkerFaceColor','k','MarkerSize',8);
xlabel('W'); ylabel('y = W'''); axis square;
title('Richest regime: 4 equilibria \rightarrow kink + bell + antibell + periodic family coexist');

fprintf('\nDone. Read the separatrix (thick contour through each saddle) as:\n');
fprintf('  - closed loop around a CENTER enclosing one saddle  -> BELL/ANTIBELL soliton\n');
fprintf('  - open branch connecting two SADDLES along y=0       -> KINK/ANTIKINK\n');
fprintf('  - nested closed contours inside a loop                -> periodic wave family\n');
