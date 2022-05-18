clear
clc
%% Default Material Model - Strain Hardening (General Eqn)
global Y0
% Y0 = input('Enter Yield Strength of Material (MPa) ');
Y0 = 75;
global n
% n = input('Enter strain hardening exponent ');
n = 0.55;
global K
% K = input('Enter strength coeff. (MPa) ');
K = 300;
%% Die Geometry 
global alpha
alpha = 8; %degrees
global mu
mu = 0.10;
%% Material Geometry (mm)
% d_o = input('Enter original diameter of wire');
global d_o
d_o = 10;
global l_o
% l_o = input('Enter original length of the wire');
l_o = 25;
%% Draw Ratio
global d_r
% d_r = input('Enter desired draw ratio ');
d_r = 1.22;
% die_sz = d_r/(2*tand(alpha));
die_sz = d_o*(1 - (1/d_r))/(2*tand(alpha));
die_sz_dp = round(die_sz, 1);
%% Euler Analysis
ys_data = zeros((10*die_sz_dp) + 1, 1);
die_pres = zeros((10*die_sz_dp) + 1, 1);
count = 1;
count_vec = linspace(0, die_sz_dp, (10*die_sz_dp) + 1);
for x = 0:0.1:die_sz_dp
    if x == 0
     ys_data(1) = Y0;   
    else
     ys_data(count) = yield_strength(x); 
    end
    count = count + 1;
end
%% Euler - With Friction
[x, sig_x] = ode45(@sigma_x, (count_vec'), 0); 
for i = 1:size(count_vec')
   die_pres(i) = ys_data(i) - sig_x(i); 
end
%% Euler - Without Friction
nofric_die_pres = zeros((10*die_sz_dp) + 1, 1);
[x_nofric, fsig_x] = ode45(@nofric_sigma_x, (count_vec'), 0);

for i = 1:size(count_vec')
   nofric_die_pres(i) = ys_data(i) - fsig_x(i); 
end
%% Plotter 
subplot(1, 2, 1);
plot(count_vec', ys_data, '--', 'color', 'black');
hold on 
plot(count_vec', sig_x, '-.', 'color', 'black');
hold on
plot(count_vec', die_pres, 'color', 'black');
% title('With friction, \mu = 0.10')
hold off
xlim([0, die_sz_dp])
xlabel('x (mm)')
ylabel('Pressure, Stress, Yield Strength (MPa)')
legend("Yield Strength", "Axial Stress", "Die Pressure", 'location', 'northwest')
legend boxoff
a = gca;
a.Position(2) = 0.17;
l1 = ['With Friction, \mu = ' num2str(mu)];	
l2 = ['\alpha = ', num2str(alpha), '^{\circ}'];
l3 = ['y = ', num2str(Y0), '+', num2str(K), '\epsilon^{', num2str(n), '}', 'MPa'];
str = {l1,l2,l3};
t = annotation('textbox', [0.05,0.02, 0.1,0.1], 'String', str, 'LineStyle', 'none');
sz1 = t.FontSize;
t.FontSize = 8;

subplot(1,2,2);
plot(count_vec', ys_data, '--', 'color', 'black');
hold on 
plot(count_vec', fsig_x, '-.', 'color', 'black');
hold on
plot(count_vec', nofric_die_pres, 'color', 'black');
% title('Without friction')
hold off
xlim([0, die_sz_dp])
xlabel('x (mm)')
legend("Yield Strength", "Axial Stress", "Die Pressure", 'location', 'northwest')
legend boxoff
b = gca;
b.Position(2) = 0.17;
l1 = ['Frictionless, \mu = ' num2str(0)];	
l2 = ['\alpha = ', num2str(alpha), '^{\circ}'];
l3 = ['y = ', num2str(Y0), '+', num2str(K), '\epsilon^{', num2str(n), '}', 'MPa'];
str = {l1,l2,l3};
t = annotation('textbox', [0.5,0.02, 0.1,0.1], 'String', str, 'LineStyle', 'none');
sz2 = t.FontSize;
t.FontSize = 8;

%% Functions
function eps = instant_eps(x)
    global d_o
    global alpha
    eps = (2*log(d_o/(d_o - 2*x*tand(alpha))));
end

function s_ys = yield_strength(x)
    global n
    global K
    global Y0
    s_ys = Y0 + K*((instant_eps(x))^n);
end

function dsigma_x = sigma_x(x, sig_x)
    global alpha
    global mu
    global d_o
    dsigma_x = (4*tand(alpha)*((yield_strength(x)*(1 + (mu*cotd(alpha)))) - (sig_x*mu*cotd(alpha))))/((d_o - 2*x*tand(alpha)));
end

function fsigma_x = nofric_sigma_x(x, fsig_x)
    global alpha
    global d_o
    fsigma_x = (4*tand(alpha)*(yield_strength(x)))/((d_o - 2*x*tand(alpha)));
end



