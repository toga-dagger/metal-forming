clear;
clc;
%% Material Properties
mu = 0.2; % Friction Coefficient
%% Material Geometry 
l = 100; %mm
h = 10; %mm
%% Transition
x_lim = round((h/(2*mu))*log(1/(2*mu)),0);
%% Dry Friction
p2k_dry = @(x)(exp((2*mu/h)*x));
%% Sticking Friction
p2k_sticking = @(x)(1 + x/h);
p2k_sticking_trans = @(x)(((1/(2*mu))*(1-log(1/(2*mu)))) + (x/h));
size = (l/2) - x_lim + 1;
%% plotter - main 
fplot(p2k_dry, [0,x_lim], 'color', 'black', 'LineWidth', 3); % Dry Friction from 0 -> x_lim
hold on
fplot(p2k_sticking_trans, [x_lim, (l/2)], 'color', 'black', 'LineWidth', 3); % Sticking Friction from x_lim -> l/2
hold on
plot(linspace(l/2, l-x_lim+1, size), flip(p2k_sticking_trans(x_lim:l/2)), 'color', 'black', 'LineWidth', 3); %Sticking Friction from l/2 -> l/2 + x_lim
hold on
plot(linspace(l - x_lim + 1, l, x_lim+1), flip(p2k_dry(0:x_lim)), 'color', 'black', 'LineWidth', 3); %Dry Friction from l/2 + x_lim -> l
hold on 
%% plotter - pseudo
fplot(p2k_sticking, [0,l/2], ':', 'color', 'black'); %Sticking Friction from 0 -> l/2
% hold on 
% plot([x_lim, x_lim], [p2k_sticking(x_lim), p2k_dry(x_lim)], ':', 'color', 'black'); %Left hand Jump
hold on 
fplot(p2k_dry, [x_lim, l/2], ':', 'color', 'black'); %Dry Friction from x_lim -> l/2
hold on
plot(linspace(l/2, l-x_lim+1, size),  flip(p2k_dry(x_lim:l/2)), ':', 'color', 'black'); %Dry Friction from l/2 -> l/2 + x_lim
hold on
% plot([l-x_lim+1, l-x_lim+1], [p2k_sticking(x_lim), p2k_dry(x_lim)], ':', 'color', 'black'); %Right hand Jump
% hold on 
plot(linspace(l/2, l, l/2 + 1), flip(p2k_sticking(0:l/2)), ':', 'color', 'black'); % Sticking Friction from l/2 -> l 
set(gca,'xtick',[],'ytick',[])

