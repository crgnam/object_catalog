matlabrc; clc; close all;

%% Load
data = readmatrix('../test_data.csv');

%% Plot
x = data(:,1);
y = data(:,2);
z = data(:,3);
bools = boolean(data(:,4));
plot3(x(~bools), y(~bools), z(~bools),'.r','MarkerSize',0.1); hold on
plot3(x(bools), y(bools), z(bools)+1e9,'.g','MarkerSize',10);
plot3(3e8, 0, 1e9, '.k','MarkerSize',20);
axis equal
grid on
xlim([-1 1]*1e9)
ylim([-1 1]*1e9)

view([0 90])