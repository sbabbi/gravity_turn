%% Compute the gravity turn profile
% This formulation of the equations is singular when v^2 = mu/r.
% Not a big deal since you are far from optimal when this happens (ideally you want to have this condition at the very end of the turn)
% Drag is not taken into account, but should not be hard to modify accordingly

% Adimensionalization:
% distances are adimensionalized with the body radius
% accelerations are adimensionalized with the body surface gravity
% velocity is adimensionalized with sqrt(gR)
% time (and ISP) is adimensionalized with sqrt(R/g)

% @param beta0 Initial angle
% @param TWR_0 Initial *surface* TWR
% @param v0 Initial velocity (default is 0.0329 that is 80m/s on kerbin)
% @param iISP inverse of the ISP (default is zero, that means mass loss is ignored)
% @param r0 Initial radius (that is initial altitude + body radius). Default is whatever you get by burning with TWR_0 until you reach your initial velocity input
% @return beta The angle between your velocity and the azimuth. Will try to get all the profile up to pi/2 unless we hit the singularity (see above)
% @return r The radius profile (body radius + altitude)
% @return theta The longitude profile (starts at 0)
% @return v The velocity profile
% @return TWR the TWR profile

% all angles in radians

function [beta, r, theta, v, TWR] = gravity_turn(beta0, TWR_0, v0, iISP, r0)
    
    if nargin < 3
        v0 = 0.032974624;
    end
    
    if nargin < 4
        iISP = 0;
    end

    vopt = odeset('RelTol', 1e-4, 'AbsTol', 1e-3, 'InitialStep', 1e-4, 'MaxStep', 1e-3, 'Events', @(beta,r) term(beta, r(1), r(3) ) );
    
    dr = @(beta,r) grav(beta,r(1), r(3), r(4), iISP);
    
    warning('off', 'OdePkg:InvalidArgument' );
    [beta, res] = ode45(dr, [beta0 0.5*pi], [1 + .5 * v0^2, 0, v0, TWR_0], vopt );
    warning('on', 'OdePkg:InvalidArgument' );
    r = res(:,1);
    theta = res(:,2);
    v = res(:,3);
    TWR = res(:, 4);
end
    
function ans = grav( beta, r, v, TWR, iISP )
    k = 1/r - v^2;
    dv = (TWR - cos(beta) / r^2) * (r * v) / (sin(beta) * k);

    dr = v^2 * r / (tan(beta) * k);
    dtheta = v^2 / k;
    dTWR = TWR^2 * iISP * (r * v) / (k* sin(beta));
    
    ans = [dr; dtheta; dv; dTWR];
end

function [value, isterminal, direction] = term(beta, r, v)
    value = 1/r - v^2;
    isterminal = 1;
    direction = 0;
end
