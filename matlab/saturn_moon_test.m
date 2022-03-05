matlabrc; clc;
addpath(genpath('mice'));
addpath(genpath('tools'));

% Furnsh all of the planet/satellite ephems:
cspice_furnsh('meta_kernel.tm')

% Planet/Moons to look at:
planet_name = 'SATURN';
moon_names = {'TITAN','RHEA','TETHYS','DIONE','IAPETUS','ENCELADUS'};
moon_i = deg2rad([0.349, 0.327, 0.168, 0.002, 15.470, 0.01, 1.566]); % TODO:  Figure out how to calculate these properly

% Obtain the gravitational constants:
mu = cspice_bodvrd( planet_name, 'GM', 1 );
J2 = cspice_bodvrd( planet_name, 'J2', 1 );
radii = cspice_bodvrd( planet_name, 'RADII', 3 );
radius = max(radii);

% Generate a time range:
num_dates = 1000;
start_date = datetime('01-Jan-2022');
end_date = datetime('01-May-2022');
et_start = cspice_str2et(datestr(start_date));
et_end = cspice_str2et(datestr(end_date));
et_range = linspace(et_start, et_end, num_dates);

tsince = et_range - et_start;

t_plt = tsince/86400;
XLABEL = 'Time since Epoch (days)';
arrayfun(@cla,findall(0,'type','axes'))

h = [];
%% Loop through all Moons:
for ii = 1:length(moon_names)
    % Convert to orbital elements:
    rv = cspice_spkezr(moon_names{ii},et_range,'J2000','NONE',planet_name);
    r = rv(1:3,:);
    v = rv(4:6,:);
    
    % Generate the orbital elements:
    [a,e,i,peri,node,M0] =  rv2kep(mu, r(:,1),v(:,1));
    
    % Apply a fix:
    n = sqrt(mu/(a^3));
    dPeri = 3*n*(radius^2)*J2*(4 - 5*sin(moon_i(ii))^2)/(4*a^2);
    peri = peri + dPeri*tsince;
    
    % Evaluate orbit:
    r_kep = kep2rv(mu, a,e,i,peri,node,M0, tsince);

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
        [a2(jj),e2(jj),i2(jj),peri2(jj),node2(jj),M02(jj)] =  rv2kep(mu, r(:,jj),v(:,jj));
    end
    r_error = r - r_kep;
    [~,mag_error] = normc(r_error);
    
    % Plot the results:
    figure(1)
        h(ii) = plot3(r_kep(1,:),r_kep(2,:),r_kep(3,:)); hold on
        plot3(r(1,:),r(2,:),r(3,:),'--k');
        plot3(r(1,end),r(2,end),r(3,end),'.k','MarkerSize',20);
        plot3(r_kep(1,end),r_kep(2,end),r_kep(3,end),'.r','MarkerSize',20);
        title([planet_name,' Moons'])
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
        ylabel('In-track (km)')
        grid on
    subplot(3,1,3)
        plot(t_plt,rsw(3,:)); hold on
        ylabel('Cross-track (km)')
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
        title(['Angular Error Relative to ',planet_name])
        ylabel('Error (deg)')
        xlabel(XLABEL);
    
    figure(5)
        plot(t_plt,mag_error); hold on
        grid on
        title(['Magnitude of Error Relative to ',planet_name])
        ylabel('Error (km)')
        xlabel(XLABEL);
        
%     figure(6)
%         subplot(3,2,1)
%             plot(t_plt,a2-a); hold on
%             ylabel('a (km)')
%             title('Changes in Elements')
%             xlabel(XLABEL)
%             grid on
%         subplot(3,2,2)
%                 plot(t_plt,e2-e); hold on
%                 ylabel('e')
%                 xlabel(XLABEL)
%                 grid on
%         subplot(3,2,3)
%                 plot(t_plt,rad2deg(i2-i)); hold on
%                 ylabel('i (deg)')
%                 xlabel(XLABEL)
%                 grid on
%         subplot(3,2,4)
%                 plot(t_plt,rad2deg(peri2-peri)); hold on
%                 ylabel('peri (deg)')
%                 xlabel(XLABEL)
%                 grid on
%         subplot(3,2,5)
%                 plot(t_plt,rad2deg(node2-node)); hold on
%                 ylabel('node (deg)')
%                 xlabel(XLABEL)
%                 grid on
end
figure(1)
    legend(h,moon_names)
    
figure(4)
    legend(moon_names)
    
figure(5)
    legend(moon_names)