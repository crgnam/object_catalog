function [rHill] = eci2hill(rTgt, vTgt, rChase)
    rTgtMag   = sqrt(sum(rTgt.^2,1));
    rChaseMag = sqrt(sum(rChase.^2,1));

    RSW = eci2rsw(rTgt, vTgt);

    % Use RSW rotation matrix to convert rChase and vChase to RSW
    r_Chase_RSW = RSW*rChase;

    % Find Rotation angles to go from target to interceptor
    phi_chase    = asin(r_Chase_RSW(3)/rTgtMag);
    lambda_chase = atan2(r_Chase_RSW(2),r_Chase_RSW(1));

    % Find Position component rotations
    rHill = cat(1, rChaseMag - rTgtMag, ...
                   lambda_chase*rTgtMag, ...
                   phi_chase*rTgtMag);
end