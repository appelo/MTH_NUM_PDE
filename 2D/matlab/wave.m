clear

Xl = -1/3;
Xr = 1/3;
Yb = -1/2;
Yt = 1/2;

Tend = 1.23;

CFL = 0.5;
% Use Twilight?
mms = 1;
% Plotting frequency
nplot = 100000;

nx = 31;
ny = 21;

hx = (Xr-Xl)/(nx-1);
hy = (Yt-Yb)/(ny-1);

dt = CFL*min([hx,hy]);

nt = ceil(Tend/dt);
dt = Tend/nt;


% Grid
x = zeros(nx,ny);
y = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        x(i,j) = Xl + (i-1)*hx;
        y(i,j) = Yb + (j-1)*hy;
    end
end

t = 0;

% Initial data
if(mms == 1)
    u = mmsfun(x,y,t,0,0,0);
    ut = mmsfun(x,y,t,1,0,0);
    um = mmsfun(x,y,t-dt,0,0,0);
else
    u = zeros(nx,ny);
    ut = zeros(nx,ny);
    lapu = compute_lap(u,x,y,t,mms);
    % FIX ME
end

if (1==1)
    err = um-mmsfun(x,y,t-dt,0,0,0);
    mesh(x,y,err)
    title('Error in initial data for um')
    disp(['Max error in initial data for um is ' num2str(max(max(abs(err))))])
    pause
end

if (1==1)
    [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);
    u = update_bc(u,ulef,urig,ubot,utop);
    err = u-mmsfun(x,y,t,0,0,0);
    mesh(x,y,err)
    title('Error in BC')
    disp(['Max error in BC is ' num2str(max(max(abs(err))))])
    pause
end

if (1==1)
    lapu = compute_lap(u,x,y,t,mms);
    err = lapu-mmsfun(x,y,t,0,2,0)-mmsfun(x,y,t,0,0,2);
    mesh(x,y,err)
    title('Error in Lap')
    disp(['Max error in Laplace ' num2str(max(max(abs(err))))])
    pause
end

if (1==1)
    [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);
    u = update_bc(u,ulef,urig,ubot,utop);
    lapu = compute_lap(u,x,y,t,mms);
    f = compute_forcing(x,y,t,mms);
    up = 2*u-um+dt^2*(lapu+f);
    err = up-mmsfun(x,y,t+dt,0,0,0);
    mesh(x,y,err)
    title('Error in first timestep')
    disp(['Max error in first timestep ' num2str(max(max(abs(err))))])
    pause
end    
    

% Main timestepping loop

for it = 1:nt
    t = (it-1)*dt;
    % Compute forcing on Boundaries 
    [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);
    % Update Boundary conditions;
    u = update_bc(u,ulef,urig,ubot,utop);
    % Compute Laplace(u)
    lapu = compute_lap(u,x,y,t,mms);
    % Compute Forcing
    f = compute_forcing(x,y,t,mms);
    % Leap frog solution
    up = 2*u-um+dt^2*(lapu+f);
    % Swap variables
    um = u;
    u = up;
    
    t = t+dt;
    if(mod(it,nplot)==0)
        % Compute forcing on Boundaries 
        [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);
        % Update Boundary conditions
        u = update_bc(u,ulef,urig,ubot,utop);
        mesh(x,y,u)
        drawnow

    end

    % Add option of error calculation every now and then.
       
    
end

% Add option of error calculation at the final time.






function u = update_bc(u,ulef,urig,ubot,utop);
    
% Update Dirichlet boundary conditions for u 
    
    u(1,:) = ulef;
    u(end,:) = urig;
    u(:,1) = ubot;
    u(:,end) = utop;
end


function [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);

% Compute Dirichletboundary condition values for u

    if(mms == 1)
        ulef = mmsfun(x(1,:),y(1,:),t,0,0,0);
        urig = mmsfun(x(end,:),y(end,:),t,0,0,0);
        ubot = mmsfun(x(:,1),y(:,1),t,0,0,0);
        utop = mmsfun(x(:,end),y(:,end),t,0,0,0);
    else
        disp(['option not implemented yet'])
        return
    end
end    


function [lapu] = compute_lap(U,x,y,t,mms);

% Compute approximate Laplacian of u.
    
    lapu = zeros(size(U));
    if(mms ==1)
        lapu = mmsfun(x,y,t,0,2,0)+mmsfun(x,y,t,0,0,2);
    else
        disp(['option not implemented yet'])
        return
    end
end

function f = compute_forcing(x,y,t,mms);

% Compute forcing
    
    if(mms ==1)
        % f = u_tt - u_xx - u_yy
        f = mmsfun(x,y,t,1,0,0)-mmsfun(x,y,t,0,2,0)-mmsfun(x,y,t,0,0,2);
    else
        disp(['option not implemented yet'])
        return
    end
    

end

function f = mmsfun(x,y,tin,td,xd,yd);
% Manufactured solution
    t = tin*ones(size(x));
    if (td == 0)
        T=1+t/2+(t.^2)/3;
    elseif(td == 1)
        T=1/2+2*t/3;
    elseif(td == 2)
        T=0*t+2/3; % 0*t so that function vectorizes
    else
        T=0*t;
    end
    if (xd == 0  && yd == 0)
        f=(x.^2+2*x.*y+y.^2).*T;
    end
    if (xd == 1  && yd == 0)
        f=(2*x+2*y).*T;
    end
    if (xd == 0  && yd == 1)
        f=(2*x+2*y).*T;
    end
    if (xd == 1  && yd == 1)
        f=(2).*T;
    end
    if (xd == 2  && yd == 0)
        f=(2).*T;
    end
    if (xd == 0  && yd == 2)
        f=(2).*T;
    end
    if (xd+yd > 2)
        f=(0*(x+y)).*T;
    end
end
