matlabrc; clc; close all;

%% Load
data = readmatrix('../test_data.csv');
data2 = readmatrix('../planet_positions.csv');

%% Plot
x = data(:,1);
y = data(:,2);
z = data(:,3);
bools = boolean(data(:,4));
plot3(x(~bools), y(~bools), z(~bools),'.r','MarkerSize',0.1); hold on
% plot3(x(bools), y(bools), z(bools)+1e9,'.g','MarkerSize',10);
plot3(x(bools), y(bools), z(bools)+1e9,'og','MarkerSize',10,'LineWidth',2);
plot3(data2(:,1),data2(:,2),data2(:,3)+1e9,'.b','MarkerSize',10);
plot3(3e8, 0, 2e9, 'xk','MarkerSize',10,'LineWidth',3);
axis equal
grid on
xlim([-1 1]*5e9)
ylim([-1 1]*5e9)

view([0 90])