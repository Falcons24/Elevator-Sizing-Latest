clc
clear

% Conversion factors
meters_to_feet = 3.28084;
y_axis_stretch_factor = 2.5;

% Original data in SI units
Vmax = 22; % m/s
Vs = 9.33; % m/s
Sw = 1.41; % m^2
Sh = 0.313; % m^2
Cbar = 0.58; % m
Tmax = 49.5; % N
rho = 1.225; % kg/m^3
Cmo = 0.0012;
zT = 0.025; % m
CLa = 0.068 * 57.3; % 1/rad
CLah = 0.061 * 57.3; %  1/rad
CLa_wf = 0.069 * 57.3;
g = 9.81; % m/s^2
m = 16; % kg
CLo = 0.956;
taw = 0.6; % Effectiveness parameter
etha_h = 0.9; % Tail Efficiency
lh = 1.108; % m from main landing gear
AR_w = 4.48; % Aspect Ratio of wing
de_da = ((CLa * 2) / (3.14 * AR_w)); % Downwash
CLdE = -CLah * etha_h * Sh * taw / Sw;

% Most aft cg
xcg = 0; % m from main landing gear
h_to_ho = (0.4 - 0.3684) / Cbar; % m
l_h1 = lh + xcg; %m
VH1 = (l_h1 * Sh) / (Sw * Cbar);
CmdE1 = -CLah * etha_h * VH1 * taw;
Cma1 = CLa_wf * h_to_ho - CLah * etha_h * Sh * (l_h1 / Cbar) * (1 - de_da) / Sw;

% Most forward cg
xcg = 0.025; % m from main landing gear
h_to_ho = (0.4 - 0.3684 - 0.025) / Cbar; % m
l_h2 = lh + xcg; % m
VH2 = (l_h2 * Sh) / (Sw * Cbar);
CmdE2 = -CLah * etha_h * VH2 * taw;
Cma2 = CLa_wf * h_to_ho - CLah * etha_h * Sh * (l_h2 / Cbar) * (1 - de_da) / Sw;

i = 1;
for U1 = Vs - 2 : 0.5 : Vmax + 2
    qbar = 0.5 * rho * U1^2;
    CL1 = (m * g) / (qbar * Sw);
    f1 = ((Tmax * zT) / (qbar * Sw * Cbar)) + Cmo;
    dE1(i) = -((f1 * CLa) + (CL1 - CLo) * Cma1) / (CLa * CmdE1 - Cma1 * CLdE);
    dE2(i) = -((f1 * CLa) + (CL1 - CLo) * Cma2) / (CLa * CmdE2 - Cma2 * CLdE);
    V(i) = U1;
    i = i + 1;
end

% Stretch and 2.5x all y-axis graph line values
dE1 = dE1 * y_axis_stretch_factor;
dE2 = dE2 * y_axis_stretch_factor;

% Plotting
figure;
plot(V * meters_to_feet, dE1 * 57.3, 'o', 'Color', [0 0.2784 0.5373], 'MarkerFaceColor', [0 0.2784 0.5373]);
hold on;
plot(V * meters_to_feet, dE2 * 57.3, '*', 'Color', [0 0.5451 0.8157]);

grid off;
xlabel('Speed (ft/s)')
ylabel('\delta_E (deg)')

% Find the values of dE1 and dE2 at Vmax and Vs
dE1_at_Vmax = interp1(V, dE1, Vmax, 'linear', 'extrap');
dE2_at_Vmax = interp1(V, dE2, Vmax, 'linear', 'extrap');
dE1_at_Vs = interp1(V, dE1, Vs, 'linear', 'extrap');
dE2_at_Vs = interp1(V, dE2, Vs, 'linear', 'extrap');

% Add vertical lines at Vmax and Vs
line([Vmax * meters_to_feet, Vmax * meters_to_feet], [0, dE2_at_Vmax * 57.3], 'Color', [0 0.1922 0.6235], 'LineStyle', '--');
line([Vs * meters_to_feet, Vs * meters_to_feet], [0, dE2_at_Vs * 57.3], 'Color', [0 0.1922 0.6235], 'LineStyle', '--');

% Add labels for the vertical lines
text(Vmax * meters_to_feet, dE2_at_Vmax * 57.3, [' Vdive = ' num2str(round(Vmax*meters_to_feet,1)) ' ft/s'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', [0 0.1922 0.6235]);
text(Vs * meters_to_feet, dE2_at_Vs * 57.3, [' Vstall = ' num2str(round(Vs*meters_to_feet,1)) ' ft/s'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', [0 0.1922 0.6235]);

% Update the legend with the correct labels
legend('Most aft cg', 'Most forward cg', 'Location', 'Best');

hold off;
