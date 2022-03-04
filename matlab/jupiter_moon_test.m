matlabrc; clc;
addpath(genpath('mice'));
addpath(genpath('tools'));

% Furnsh all of the planet/satellite ephems:
cspice_furnsh('meta_kernel.tm')

% Obtain the gravitational constants:
mu_jupiter    = cspice_bodvrd( 'JUPITER', 'GM', 1 );
radii_jupiter = cspice_bodvrd( 'JUPITER', 'RADII', 3 );
radius_jupiter = max(radii_jupiter);

% J2 coefficient for Jupiter:
J2_jupiter = 0.01475;

% Generate a time range:
num_dates = 1000;
start_date = datetime('01-Jan-2022');
end_date = datetime('01-Feb-2022');
et_start = cspice_str2et(datestr(start_date));
et_end = cspice_str2et(datestr(end_date));
et_range = linspace(et_start, et_end, num_dates);

tsince = et_range - et_start;

t_plt = tsince/86400;
XLABEL = 'Time since Epoch (days)';
moon_names = {'IO','EUROPA','GANYMEDE','CALLISTO'};

inc_rate  = [1.6834e-11, 5.3868e-11, -2.1560e-12, 2.6934e-12];
node_rate = [-3.3668e-12, 4.0401e-11, 6.7335e-12, -6.7335e-12];
peri_rate = [0, 0, 0, 0];
arrayfun(@cla,findall(0,'type','axes'))

%% Loop through all Moons:
for ii = 1:length(moon_names)
    % Convert to orbital elements:
    rv = cspice_spkezr(moon_names{ii},et_range,'J2000','NONE','JUPITER');
    r = rv(1:3,:);
    v = rv(4:6,:);
    
    % Generate the orbital elements:
    [a,e,i,peri,node,M0] =  rv2kep(mu_jupiter, r(:,1),v(:,1));
    
    % Apply general perturbation models:
%     dNode = -(3/2)*J2_jupiter*(radius_jupiter/(a*(1-e^2)))^2*sqrt(mu_jupiter/(a^3))*cos(i);
    node = node - node_rate(ii)*tsince;
    i = i - inc_rate(ii)*tsince;
    
    % Evaluate orbit:
    r_kep = kep2rv(mu_jupiter, a,e,i,peri,node,M0, tsince);

    % Calculate error metrics:
    rsw = zeros(size(r_kep));
    angular_error = zeros(1,length(r_kep));
    a2 = zeros(1,length(r_kep));
    e2 = zeros(1,length(r_kep));
    i2 = zeros(1,length(r_kep));
    peri2 = zeros(1,length(r_kep));
    node2 = zeros(1,length(r_kep));
    M02 = zeros(1,length(r_kep));
    for jj = 1:length(tsince)
        % Calcualte RSW error:
        rsw(:,jj) = eci2hill(r(:,jj), v(:,jj), r_kep(:,jj));

        % Calculate angular error in position:
        angular_error(jj) = acos(dot(r(:,jj), r_kep(:,jj))/(norm(r(:,jj))*norm(r_kep(:,jj))) );
        
        % Calculate the time history of orbital elements:
        [a2(jj),e2(jj),i2(jj),peri2(jj),node2(jj),M02(jj)] =  rv2kep(mu_jupiter, r(:,jj),v(:,jj));
    end
    r_error = r - r_kep;
    [~,mag_error] = normc(r_error);
    
    % Plot the results:
    figure(1)
        plot3(r(1,:),r(2,:),r(3,:),'k'); hold on
        plot3(r_kep(1,:),r_kep(2,:),r_kep(3,:),'--r');
        plot3(r(1,end),r(2,end),r(3,end),'.k','MarkerSize',20);
        plot3(r_kep(1,end),r_kep(2,end),r_kep(3,end),'.r','MarkerSize',20);
        title('Jovian Moons')
        axis equal
        grid on

    figure(2)
    subplot(3,1,1)
        plot(t_plt,rsw(1,:)); hold on
        title('Errors in Kepler vs SPK (RSW)')
        ylabel('Radial (km)')
        grid on
    subplot(3,1,2)
        plot(t_plt,rsw(2,:)); hold on
        ylabel('Cross-track (km)')
        grid on
    subplot(3,1,3)
        plot(t_plt,rsw(3,:)); hold on
        ylabel('In-track (km)')
        grid on
        xlabel(XLABEL)

    figure(3)
    subplot(3,1,1)
        plot(t_plt,r_error(1,:)); hold on
        title('Errors in Kepler vs SPK (Cartesian)')
        ylabel('X (km)')
        grid on
    subplot(3,1,2)
        plot(t_plt,r_error(2,:)); hold on
        ylabel('Y (km)')
        grid on
    subplot(3,1,3)
        plot(t_plt,r_error(3,:)); hold on
        ylabel('Z (km)')
        grid on
        xlabel(XLABEL)

    figure(4)
        plot(t_plt,rad2deg(angular_error)); hold on
        grid on
        title('Angular Error Relative to Jupiter')
        ylabel('Error (deg)')
        xlabel(XLABEL);
    
    figure(5)
        plot(t_plt,mag_error); hold on
        grid on
        title('Magnitude of Error Relative to Jupiter')
        ylabel('Error (km)')
        xlabel(XLABEL);
        
    figure(6)
        subplot(3,2,1)
            plot(t_plt,a2-a); hold on
            ylabel('a (km)')
            title('Changes in Elements')
            xlabel(XLABEL)
            grid on
        subplot(3,2,2)
                plot(t_plt,e2-e); hold on
                ylabel('e')
                xlabel(XLABEL)
                grid on
        subplot(3,2,3)
                plot(t_plt,rad2deg(i2-i)); hold on
                ylabel('i (deg)')
                xlabel(XLABEL)
                grid on
        subplot(3,2,4)
                plot(t_plt,rad2deg(peri2-peri)); hold on
                ylabel('peri (deg)')
                xlabel(XLABEL)
                grid on
        subplot(3,2,5)
                plot(t_plt,rad2deg(node2-node)); hold on
                ylabel('node (deg)')
                xlabel(XLABEL)
                grid on
end

legend(moon_names)