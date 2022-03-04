function [a,e,i,peri,node,M0] = rv2kep(mu,r,v)
    % Calculate Angular Momentum:
    H = cross(r,v);

    % Calculate node vector:
    n = cross([0 0 1],H);

    % Calculate eccentricity:
    evec = ((norm(v)^2-mu/norm(r))*r - dot(r,v)*v)/mu;
    e    = norm(evec);

    % Calculate specific energy:
    energy = (norm(v)^2)/2 - mu/norm(r);

    % Calculate semi-major axis:
    a = -mu/(2*energy);

    % Calculate inclination:
    i = acos(H(3)/norm(H));

    % Calculate RAAN:
    node = acos(n(1)/norm(n));
    if n(2) < 0
       node = 2*pi-node;
    end

    % Calculate Argument of Periapsis:
    peri = acos(dot(n,evec)/(norm(n)*e));
    if evec(3) < 0
       peri = 2*pi-peri;
    end

    % Calculate true anomaly:
    nu = acos(dot(evec,r)/(e*norm(r)));
    if dot(r,v)<0
       nu = 2*pi - nu;
    end

    % Calculate mean anomaly (Approximation):
    M0 = nu - 2*e*sin(nu) + ((3/4)*(e^2) + (1/8)*e^4)*sin(2*nu) -...
         (1/3)*(e^3)*sin(3*nu) + (5/32)*(e^4)*sin(4*nu);
end