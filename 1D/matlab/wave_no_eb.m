clear

% Physical parameters
L = 3;
Tend = 5.7;

CFL = 0.9;
nx = 10;

nplt = 10;

h = L/(nx-1);
k = CFL*h;
nsteps = ceil(Tend/k);
k = Tend / nsteps

% The grid
x = h*(1:nx-1)';

% Initial data
u = init_cond(x);
um = u - k*init_velocity(x) + k^2/2*(compute_uxx(u,h)+forcing(x,0));

um = update_bc(um,x,-k);

% Start the time loop
for it = 1:nt
    t = (it-1)*k;
    
    % make sure the boundary conditions are 
    % set
    u = update_bc(u,x,t);
    % compute terms in the right hand side
    uxx = compute_uxx(u,h);
    f = forcing(x,t);
    % Leap-frog to the next time level
    up = 2*u-um + k^2*(uxx+f);
    % Swap the time levels
    um = u;
    u = up;
    
    %
    if(mod(it,nplt) == 0)
        u = update_bc(u,x,t+k);
        plot(x,u,'k','linewidth',2)
        axis([0 L -1.2 1.2])
        drawnow
    end
        
    
end

function u=init_cond(x);
   u = sin(pi*x/3.0);
end
function u=init_velocity(x);
    u = sin(2*pi*x);
end

function u = update_bc(u,x,t);
% This function returns the solution array with correctly imposed
% boundary conditions 
    u(1) = h_bc(x(1),t);
    u(end) = h_bc(x(end),t);
end

function u = h_bc(x,t)
% This function returns u at x and t.
% x is supposed to be on the boundary
% h_0, h_1 are our boundary conditions
    
    L = x(end);
    
    if x < L/2
        u = h_0(t);
    else
        u = h_1(t);
    end
    
end


 function u = h_0(t)
 % This function returns u at t on the boundary x = 0

 end

 function u = h_1(t)
 % This function returns u at t on the boundary x = L

 end


function u = h_0(t)
% This function returns u at t on the boundary x = 0
    u = 5*sin(t)-pi*cos(t);
end

function u = h_1(t)
% This function returns u at t on the boundary x = L
    u = 4*sin(t-pi/3)-2*pi*cos(t-pi/3);

end

function f = forcing(x,t)
% This function returns the right hand side forcing 
% to the wave equation 
f = sin(10*t)*exp(-(x-1/2).^2);
end

function uxx = compute_uxx(u,h)
% This function returns the second derivative 
% at all interior points 

% Get the dimension of u
    [dim_x, dim_y]  = size(u);
    % initialize the uxx matrix
    uxx = zeros(dim_x, dim_x);
    % Here I was more thinking something along the lines 
    uxx = zeros(dim_x,1);  
    ih2 = 1/h^2;
    for ix = 2:dim_x-1
        uxx(ix) = ih2*(u(ix+1)-2*u(ix)+u(ix-1));
    end

end
