% script name: "positivity_verification_v2"
% ----------------------------------------
% For: "Geometric Hermite Interpolation in $\mathbb{R}^n$"
%
% We exhaustively check the domain of theta, theta_0, and theta_1
% to verify that
% f = alpha*(x^2+y^2) - F(theta, theta_0, theta_1)>0 
% where 0<alpha<1
%
% ----------------------------------------
tic

% main parameters
gradient = 1000; %6.3*10^4;       % THE CURRENT proven BOUND IS 6.3*10^4;
alpha = .99;
fact  = 1/(gradient*3^0.5);
% we must make sure: step <= (f_ref/(gradient*3^0.5))
safe_ball = 0.1;    % The neighbourhood (ball) around zero where we do not need to check

% main loop
counter = 0;
current_x    = safe_ball;
while current_x<(.75*pi)
    current_y = safe_ball;
    while current_y<(.75*pi)
        % the ball constraint
        r = current_x^2 + current_y^2;
        if (r^0.5>3*pi/4)
            break   % current_y is too large
        end
        current_z = safe_ball;
        while current_z<(pi)
            %[current_x, current_y, current_z]
            counter = counter + 1;
            if and(current_z+current_x>=current_y, current_z+current_y>=current_x)
                if (current_x+current_y>=current_z)
                    pos_val = eval_diff(alpha, r, current_x, current_y, current_z);
                    
                    if (pos_val<0)
                        [current_x,current_y,current_z, pos_val]
                        ME = 'We found a negative point!';
                        error(ME)
                    end
                end  
            else
                break % z is too large
            end
            step_size = pos_val*fact;
            current_z = current_z + step_size;
        end
        current_y = current_y + step_size;
    end   
    current_x = current_x + step_size;
end
% search is over:
fprintf('We checked overall %d values \n', counter)

fprintf('Success!\n')
toc()




%=============================================
function [d] = eval_diff(alpha, r, x, y, z)
c = cos((x+y)/4);
cx = cos(x);
cy = cos(y);
cz = cos(z);
theta00 = (4*c^2*cx+1-cz)/(16*c^4+2-2*cz+8*c^2*(cx-cy))^0.5;
theta01 = (12*c^2+cx-5*cos(y))*c^2;
theta01 = theta01/(8*c^4+1-cz+4*c^2*(cx-cy))^0.5;
theta01 = theta01/(18*c^4+1+cz-6*c^2*(cx+cy))^0.5;

d = alpha*r-(acos(theta00)^2+acos(theta01)^2);
end
%=============================================
