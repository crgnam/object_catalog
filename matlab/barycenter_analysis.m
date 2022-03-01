matlabrc; clc; close all;
addpath(genpath('mice'));

% Furnsh all of the planet/satellite ephems:
cspice_furnsh('meta_kernel.tm')

% Define planetary radii:
Rmercury = 2440;
Rvenus = 6052;
Rearth = 6371;
Rmars = 3390;
Rjupiter = 69911;
Rsaturn = 58232;
Ruranus = 25362;
Rneptune = 24622;
Rpluto = 1188;

% Generate a time range:
num_dates = 1000;
start_date = datetime('01-Jan-2022');
end_date = datetime('31-Dec-2022');

et_start = cspice_str2et(datestr(start_date));
et_end = cspice_str2et(datestr(end_date));
et_range = linspace(et_start, et_end, num_dates);

% Calculate the motion of the planets relative to their barycenter:
r = cspice_spkpos('MERCURY',et_range,'J2000','NONE','MERCURY_BARYCENTER');
[~,r_mercury] = normc(r);

r = cspice_spkpos('VENUS',et_range,'J2000','NONE','VENUS_BARYCENTER');
[~,r_venus] = normc(r);

r = cspice_spkpos('EARTH',et_range,'J2000','NONE','EARTH_BARYCENTER');
[~,r_earth] = normc(r);

r = cspice_spkpos('MARS',et_range,'J2000','NONE','MARS_BARYCENTER');
[~,r_mars] = normc(r);

r = cspice_spkpos('JUPITER',et_range,'J2000','NONE','JUPITER_BARYCENTER');
[~,r_jupiter] = normc(r);

r = cspice_spkpos('SATURN',et_range,'J2000','NONE','SATURN_BARYCENTER');
[~,r_saturn] = normc(r);

r = cspice_spkpos('URANUS',et_range,'J2000','NONE','URANUS_BARYCENTER');
[~,r_uranus] = normc(r);

r = cspice_spkpos('NEPTUNE',et_range,'J2000','NONE','NEPTUNE_BARYCENTER');
[~,r_neptune] = normc(r);

% r = cspice_spkpos('PLUTO',et_range,'J2000','NONE','PLUTO_BARYCENTER');
% [~,r_pluto] = normc(r);


% Plot the results:
LW = 2;
t_plt = linspace(start_date, end_date, num_dates);
semilogy(t_plt,100*r_mercury/Rmercury,'LineWidth',LW); hold on
semilogy(t_plt,100*r_venus/Rvenus,'LineWidth',LW);
semilogy(t_plt,100*r_earth/Rearth,'LineWidth',LW);
semilogy(t_plt,100*r_mars/Rmars,'LineWidth',LW);
semilogy(t_plt,100*r_jupiter/Rjupiter,'LineWidth',LW);
semilogy(t_plt,100*r_saturn/Rsaturn,'LineWidth',LW);
semilogy(t_plt,100*r_uranus/Ruranus,'LineWidth',LW);
semilogy(t_plt,100*r_neptune/Rneptune,'LineWidth',LW);
% semilogy(t_plt,100*r_pluto/Rpluto,'LineWidth',LW);

legend('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune')

grid on
ylabel('Percentage of Radius')
xlabel('Date')
title('Distance away from Barycenter')