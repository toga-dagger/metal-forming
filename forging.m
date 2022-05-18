clear;
clc;
%% Material Properties
mu = 0.25; % Friction Coefficient
%% Material Geometry 
R = 150; %mm
h = 50; %mm
%% Transition
r_lim = round(R + ((h/(2*mu))*log(1.732*mu)),0);
%% Dry Friction
p2k_dry = @(r)(0.866*exp((2*mu/h)*(R-r)));
%% Sticking Friction
p2k_sticking = @(r)(0.287 + 0.333*(R-r)/h);
p2k_sticking_trans = @(r)(((1/(2*mu))*(1+log(1.732*mu))) + ((R-r)/h));
size = (R) - r_lim + 1;
%% plotter - main 
fplot(p2k_sticking_trans, [0,r_lim], 'color', 'black', 'LineWidth', 3); % Sticking Friction from 0 -> r_lim
hold on
fplot(p2k_dry, [r_lim, R], 'color', 'black', 'LineWidth', 3); % Dry Friction from r_lim -> r
%% plotter - pseudo
fplot(p2k_dry, [0,R], ':', 'color', 'black'); %Dry Friction from 0 -> r_lim
hold on 
fplot(p2k_sticking, [0, R], ':', 'color', 'black'); %Sticking Friction from 0 -> r_lim
% plot([r_lim, r_lim], [p2k_sticking(r_lim), p2k_dry(r_lim)], ':', 'color', 'black'); %Jump
set(gca,'xtick',[],'ytick',[])