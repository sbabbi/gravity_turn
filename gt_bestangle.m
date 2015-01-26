
%%
% Compute the optimal starting angle for the input initial TWR, ISP and initial velocity
function [beta, r, theta, v] = gt_bestangle( TWR_0, iISP, v0 )

    if nargin < 2
        iISP = 0;
    end

    if nargin < 3
        v0 = 0.032974624;
    end
    
    beta = fzero( @(beta) gt_overshoot(beta, TWR_0, v0, iISP), [0.1*pi/180, 20*pi/180] );
    
    [beta, r, theta, v, TWR] = gravity_turn( beta, TWR_0, v0, iISP );
    
    figure;
    plot( theta .* r, r );
    xlabel('x','Interpreter','tex');
    ylabel('r','Interpreter','tex');

    figure;
    plot( beta * 180 / pi, v );
    xlabel('\beta','Interpreter','tex');
    ylabel('v','Interpreter','tex');
    
    figure;
    plot( theta .* r, TWR, '-r' );
    xlabel('x','Interpreter','tex');
    ylabel('TWR','Interpreter','tex');
    
    disp( sprintf( 'Optimal angle is %fdeg', beta(1) * 180 / pi) )

end

%% Get the extra v with respect to a circular orbit a the end of the turn.
% penalize if we hit the singularity
function ans = gt_overshoot(beta0, TWR0, v0, iISP, penalize)
    
    if nargin < 4
        iISP = 0
    end
    
    [beta, r, theta, v] = gravity_turn(beta0,TWR0,v0,iISP);
    
    if nargin < 5
        penalize = 16.0;
    end
    
    ans = v(end)^2 - 1 / r(end) + penalize * (.5 * pi - beta(end) );
end
