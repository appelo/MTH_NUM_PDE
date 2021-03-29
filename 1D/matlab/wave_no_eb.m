
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
  
  %change this condition for more complicated shape
  boundary_left = 1;
  boundary_right = dim_x;
  
  for ix = 1:dim_x
      if ix ~= boundary_left && ix ~= boundary_right
          %if the point is regular
          uxx(ix, ix-1:ix+1) = [1,-2, 1];
      elseif ix == boundary_left
          uxx(ix, ix:ix+1) = [-2, 1];
          %impose some boundary condition for the left point
          uxx(ix, end) = 1;   %<-- For circular boundary condition
          %uxx(ix, end) = 0;   %<-- For dirchlet condition
          %uxx(ix, 1:s) == []  %<-- if using some interpolation
      elseif ix == boundary_right
          uxx(ix, ix-1:ix) = [1, -2];
          uxx(ix, 1) = 1; %<-- For circular boundary condition
          %uxx(ix, 1) = 0 %<-- For dirchlet condition
          %uxx(ix, ix-s:ix) = [] %<--- if using some interpolation
      end
  end
  uxx = uxx / h^2;
        
      
end
