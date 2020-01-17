%% AUTHOR    : Eric Graves
%% $DATE     : 07-Dec-2018 21:18:03 $
%% $Revision : 1.00 $
%% DEVELOPED : 9.4.0.813654 (R2018a)
%% FILENAME  : single_driver_bass_response.m

% Engineering Acoustics Homework 7 parts b) and c) : modelling volume velocities
% and SPL for a cabinet with and without a port.
% Includes modified code from HW6 solutions.
clear
clc
close all;

rho = 1.2;
co = 344;
Zo = rho*co;

diam = 0.107; %diameter of the driver
Ad = pi * diam^2 / 4; %cross sectional area of the driver
a = diam/2;

m = 0.009; %moving mass is 9 grams
Cm = 740e-6; %compliance in meters / Newton
omres_driver = sqrt(1/(m*Cm)); %mech resonance of driver
fres_driver = 1/(2*pi) * omres_driver
Qfact = 2.21; %quality factor of resonance for driver
Rm = omres_driver * m/Qfact; %obtaining Rm from Qfact (see course notes)

Vs = 2.83; % RMS voltage

BL = 5.9; %BL constant for the driver system
Re = 6; %resistance in ohms
L = 0.27e-3; %inductance in Henry
Vencl = 0.0638; % enclosed volume in m^2

Ca = Vencl / (Zo * co); % acoustic compliance
Zarp = 885;
Zard = 885;
Rap = 0;
ma = 35.21;

f = logspace(0,3,500);
Qd = zeros(length(f), 1);
Qp = zeros(length(f), 1);
Qtot = zeros(length(f), 1);
P = zeros(length(f),1);

for n = 1:length(f)

    om = 2*pi*f(n);
    F = BL*Vs / (Re + i*om*L);
    Pd = F / Ad;

    Zmd = Rm + i*om*m + 1/(i*om*Cm) + Zpar2(BL^2/(i*om*L),BL^2/Re);

    Z1 = Zard + Zmd / (Ad^2);
    Z2 = Zarp + Rap + i*om*ma; % 1e9 used for closed cabinet
    Z3 = 1 / (i*om*Ca);

    Z_MAT = [(Z1 + Z3) Z1; Z1 (Z1 + Z2)]
    P_MAT = [Pd; Pd];

    Q_MAT = Z_MAT \ P_MAT;

    Qd(n) = Q_MAT(1) + Q_MAT(2);
    Qp(n) = - Q_MAT(2);

    Qtot(n) = Qd(n) + Qp(n);

    P(n) = i*rho*om*Qtot(n) / (4*pi*1); % at 1 meter radius

end

figure (1)
loglog(f,abs(Qtot),'b');
grid on;
hold on;
loglog(f,abs(Qp),'r');
loglog(f,abs(Qd),'m');
xlim([1 1000]);
xlabel('Frequency (Hz)');
ylabel('Volume Velocity [m^3 / sec]');
legend('Qtot', 'Qp', 'Qd');
title('Volume Velocities vs Frequency');
pedit;

figure(2);
%SPL = 20*log10(abs(P)/20e-6); % Used for calculating bass reflex SPL
load('bass_response_spl.mat'); % load bass reflex SPL for plotting
SPL_closed = 20*log10(abs(P)/20e-6);
semilogx(f,SPL,'r'); hold on;
semilogx(f,SPL_closed,'k');
grid on;
xlim([20 1000]);
ylim([50 100]);
xlabel('Frequency (Hz)');
ylabel('SPL dB');
legend('Bass Reflex', 'Closed Cabinet');
title('SPL vs Frequency');
pedit;

% ===== EOF ====== [hw7.m] ======
