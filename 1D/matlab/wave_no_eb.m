
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
% AG


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



function u = update_bc(u,x,t);
% This function returns the solution array with correctly imposed
% boundary conditions 
   % DEAA 
end

function u = h_bc(x,t)
% This function returns u at x and t.
% x is supposed to be on the boundary
% h_0, h_1 are our boundary conditions
if x < L/2
    u = h_0(t);
else
    u = h_1(t);
end

end

function f = forcing(x,t)
% This function returns the right hand side forcing 
% to the wave equation 
f=sin(10*t)*exp(-(x-1/2)^2);
end

function uxx = compute_uxx(u,h)
% This function returns the second derivative 
% at all interior points 
% (but uxx has the same dimension as u)   <-- I think we need the x, not
% the u
  %YH  
  %Get the dimension of u
  [dim_x, dim_y]  = size(u);
  %initialize the uxx matrix
  uxx = zeros(dim_x, dim_x);

  % Here I was more thinking something along the lines 
  uxx = zeros(dim_x,1);  
  ih2 = 1/h^2;
  for ix = 2:dim_x-1
        uxx(ix) = ih2*(u(ix+1)-2*u(ix)+u(ix-1));
  end
  uxx(1) = ih2 *(-2*u(1) + u(2));
  uxx(end) = ih2 *(-2u(end) + u(end-1) );
  
  
  
        
      
end
