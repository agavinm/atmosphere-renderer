close all;
%% Load data z
data = load('dataZ.txt');
x = data(:,1); %% z
y1 = data(:,2); %% rayleigh_scattering
y2 = data(:,3); %% aerosol_scattering
y3 = data(:,4); %% ozone_absorption
y4 = data(:,5); %% aerosol_absorption
y5 = data(:,6); %% extinction
y6 = data(:,7); %% rayleigh_density
y7 = data(:,8); %% ozone_density_1
y8 = data(:,9); %% ozone_density_2
y9 = data(:,10); %% ozone_density_3
y10 = data(:,11); %% ozone_density_4
y11 = data(:,12); %% ozone_density_5
y12 = data(:,13); %% ozone_density_6
y13 = data(:,14); %% ozone_density_7
y14 = data(:,15); %% ozone_density_8
y15 = data(:,16); %% ozone_density_9
y16 = data(:,17); %% ozone_density_10
y17 = data(:,18); %% ozone_density_11
y18 = data(:,19); %% ozone_density_12
y19 = data(:,20); %% aerosol_density

%% Plot
figure;
hold on;
ylabel('rayleigh scattering (km^{-1})'); xlabel('z (km)');
plot(x, y1, 'LineWidth', 2);

figure;
hold on;
ylabel('aerosol scattering (km^{-1})'); xlabel('z (km)');
plot(x, y2, 'LineWidth', 2);

figure;
hold on;
ylabel('ozone absorption (km^{-1})'); xlabel('z (km)');
plot(x, y3, 'LineWidth', 2);

figure;
hold on;
ylabel('aerosol absorption (km^{-1})'); xlabel('z (km)');
plot(x, y4, 'LineWidth', 2);

figure;
hold on;
ylabel('extinction (km^{-1})'); xlabel('z (km)');
plot(x, y5, 'LineWidth', 2);


figure;
hold on;
ylabel('km^{-1}'); xlabel('z (km)');
plot(x, y1, 'LineWidth', 2, x, y2, 'LineWidth', 2, x, y3, 'LineWidth', 2, x, y4, 'LineWidth', 2, x, y5, 'LineWidth', 2);
ylim([0 0.015]);
xlim([0 60]);
legend('rayleigh scattering (km^{-1})', 'aerosol scattering (km^{-1})', 'ozone absorption (km^{-1})', 'aerosol absorption (km^{-1})', 'extinction (km^{-1})');

figure;
hold on;
ylabel('km^{-1}'); xlabel('z (km)');
plot(x, y1, 'LineWidth', 2, x, y2, 'LineWidth', 2, x, y3, 'LineWidth', 2, x, y4, 'LineWidth', 2, x, y5, 'LineWidth', 2);
%ylim([0 0.0015]);
legend('rayleigh scattering (km^{-1})', 'aerosol scattering (km^{-1})', 'ozone absorption (km^{-1})', 'aerosol absorption (km^{-1})', 'extinction (km^{-1})');

figure;
hold on;
ylabel('rayleigh density (km^{-3})'); xlabel('z (km)');
plot(x, y6, 'LineWidth', 2);

figure;
hold on;
ylabel('ozone density (km^{-3})'); xlabel('z (km)');
plot(x, y7, 'LineWidth', 2, x, y8, 'LineWidth', 2, x, y9, 'LineWidth', 2, x, y10, 'LineWidth', 2, x, y11, 'LineWidth', 2, x, y12, 'LineWidth', 2, x, y13, 'LineWidth', 2, x, y14, 'LineWidth', 2, x, y15, 'LineWidth', 2, x, y16, 'LineWidth', 2, x, y17, 'LineWidth', 2, x, y18, 'LineWidth', 2);
%ylim([0 10e15]);
legend('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');

figure;
hold on;
ylabel('aerosol density (km^{-3})'); xlabel('z (km)');
plot(x, y19, 'LineWidth', 2);
%ylim([0 10e15]);


%% Load data wl
data = load('dataWl.txt');
x = data(:,1); %% wl
y1 = data(:,2); %% rayleigh_cross_section_scattering
y2 = data(:,3); %% ozone_cross_section_absorption
y3 = data(:,4); %% aerosol_cross_section_scattering
y4 = data(:,5); %% aerosol_cross_section_absorption

%% Plot
figure;
hold on;
ylabel('rayleigh cross section scattering (km^{2})'); xlabel('wl (nm)');
plot(x, y1, 'LineWidth', 2);

figure;
hold on;
ylabel('ozone cross section absorption (km^{2})'); xlabel('wl (nm)');
plot(x, y2, 'LineWidth', 2);

figure;
hold on;
ylabel('aerosol cross section scattering (km^{2})'); xlabel('wl (nm)');
plot(x, y3, 'LineWidth', 2);

figure;
hold on;
ylabel('aerosol cross section absorption (km^{2})'); xlabel('wl (nm)');
plot(x, y4, 'LineWidth', 2);

figure;
hold on;
ylabel('km^{2}'); xlabel('wl (nm)');
plot(x, y3, 'LineWidth', 2, x, y4, 'LineWidth', 2);
legend('aerosol cross section scattering (km^{2})', 'aerosol cross section absorption (km^{2})');
