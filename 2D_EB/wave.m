clear
clf
Xl = -1;
Xr = 1;
Yb = -1;
Yt = 1;

Tend = 1.23;

CFL = 0.25;
% Use Twilight?
mms = 1;
% Plotting frequency
nplot = 1;

nx = 21;
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

mask = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        if(level(x(i,j),y(i,j)) < 0) % if inside the computational domain
            mask(i,j) = 1;
        end
    end 
end

maskin = mask;
signal = 0;
nbp = 0; %number of boundary points
for i = 1:nx
    for j = 1:ny
        if(mask(i,j)==1)
            signal=mask(i-1,j)+mask(i+1,j)+mask(i,j-1)+mask(i,j+1);
            if (signal < 4)
                maskin(i,j)=-1;
                nbp=nbp+1;
            end
        end
    end
end
xm = x./(mask==1);
ym = y./(mask==1);


xb = x./(maskin==-1);
yb = y./(maskin==-1);
%clf
plot(xm,ym,'k',xm',ym','k',xm,ym,'bo',xb,yb,'r*','linewidth',2)
axis equal
axis([Xl Xr Yb Yt])


[xx,yy] = meshgrid(linspace(-1,1,501));
zz = level(xx,yy);
hold on
contour(xx,yy,zz,[0,0],'r','linewidth',2)

%data structure that contains the ghost and interior point for
%interpolation
ghost_points=zeros(nbp,4);
count=1;
for i = 1:nx
    for j = 1:ny
        if(maskin(i,j)==-1)
            % Index coordinates for ghost points
           ghost_points(count,1)=i;
           ghost_points(count,2)=j;
           count = count+1;
        end
    end
end

bcfs = zeros(nbp,4); %data structure containing x_gam, y_gam, c_int ,c_gam

for i = 1:nbp
    igp = ghost_points(i,1); 
    jgp = ghost_points(i,2);
    if (level(x(igp-1,jgp),y(igp-1,jgp))>=0)
        % Look to the left for a root
        if (level(x(igp,jgp-1),y(igp,jgp-1))>0)
            % look down for root
            fun = @(x)level(x,y(igp,jgp));
            xgam1 = fzero(fun,[x(igp-1,jgp) x(igp,jgp)]);
            ygam1 = y(igp,jgp);
            fun = @(y)level(x(igp,jgp),y);
            xgam2 = x(igp,jgp);
            ygam2 = fzero(fun,[y(igp,jgp-1) y(igp,jgp)]);
            if(abs(xgam1-x(igp,jgp)) < abs(ygam2-y(igp,jgp)))
                plot(xgam1,ygam1,'rp','markersize',14,'linewidth',2)
                xgam = xgam1; ygam = ygam1;
                ghost_points(i,3)=igp+1;
                ghost_points(i,4)=jgp;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(xgam-x(igp,jgp))/(xgam-xin);
                bcfs(i,4)=(x(igp,jgp)-xin)/(xgam-xin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            else
                plot(xgam2,ygam2,'gp','markersize',14,'linewidth',2)
                xgam = xgam2; ygam = ygam2;
                ghost_points(i,3)=igp;
                ghost_points(i,4)=jgp+1;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(ygam-y(igp,jgp))/(ygam-yin);
                bcfs(i,4)=(y(igp,jgp)-yin)/(ygam-yin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            end
        elseif (level(x(igp,jgp+1),y(igp,jgp+1))>0)            
            % look up for root
            fun = @(x)level(x,y(igp,jgp));
            xgam1 = fzero(fun,[x(igp-1,jgp) x(igp,jgp)]);
            ygam1 = y(igp,jgp);
            fun = @(y)level(x(igp,jgp),y);
            %here is looking up for the root
            xgam2 = x(igp,jgp);
            ygam2 = fzero(fun,[y(igp,jgp) y(igp,jgp+1)]);
            if(abs(xgam1-x(igp,jgp)) < abs(ygam2-y(igp,jgp)))
                xgam = xgam1; ygam = ygam1;
                plot(xgam1,ygam1,'rs','markersize',14,'linewidth',2)
                ghost_points(i,3)=igp+1;
                ghost_points(i,4)=jgp;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(xgam-x(igp,jgp))/(xgam-xin);
                bcfs(i,4)=(x(igp,jgp)-xin)/(xgam-xin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            else
                plot(xgam2,ygam2,'gs','markersize',14,'linewidth',2)
                xgam = xgam2; ygam = ygam2;
                ghost_points(i,3)=igp;
                ghost_points(i,4)=jgp-1;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(ygam-x(igp,jgp))/(ygam-yin);
                bcfs(i,4)=(y(igp,jgp)-xin)/(ygam-yin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            end
        else
            % If all the above fails just find root to the left! 
            fun = @(x)level(x,y(igp,jgp));
            xgam = fzero(fun,[x(igp-1,jgp) x(igp,jgp)]);
            ygam = y(igp,jgp);
            plot(xgam,ygam,'bp','markersize',14)
            ghost_points(i,3)=igp+1;
            ghost_points(i,4)=jgp;
            bcfs(i,1)=xgam;
            bcfs(i,2)=ygam;
            xin=x(ghost_points(i,3),ghost_points(i,4));
            yin=y(ghost_points(i,3),ghost_points(i,4));
            bcfs(i,3)=(xgam-x(igp,jgp))/(xgam-xin);
            bcfs(i,4)=(x(igp,jgp)-xin)/(xgam-xin);
            plot(xin,yin,'g*','markersize',14,'linewidth',2)
       end
    elseif (level(x(igp+1,jgp),y(igp+1,jgp))>=0)
        % Look to the right for a root
        if (level(x(igp,jgp-1),y(igp,jgp-1))>0)
            % look down for root
            fun = @(x)level(x,y(igp,jgp));
            xgam1 = fzero(fun,[x(igp,jgp) x(igp+1,jgp)]);
            ygam1 = y(igp,jgp);
            fun = @(y)level(x(igp,jgp),y);
            xgam2 = x(igp,jgp);
            ygam2 = fzero(fun,[y(igp,jgp-1) y(igp,jgp)]);
            if(abs(xgam1-x(igp,jgp)) < abs(ygam2-y(igp,jgp)))
                plot(xgam1,ygam1,'rp','markersize',14,'linewidth',2)
                xgam = xgam1; ygam = ygam1;
                ghost_points(i,3)=igp-1;
                ghost_points(i,4)=jgp;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(xgam-x(igp,jgp))/(xgam-xin);
                bcfs(i,4)=(x(igp,jgp)-xin)/(xgam-xin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            else
                plot(xgam2,ygam2,'gp','markersize',14,'linewidth',2)
                xgam = xgam2; ygam = ygam2;
                ghost_points(i,3)=igp;
                ghost_points(i,4)=jgp+1;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(ygam-y(igp,jgp))/(ygam-xin);
                bcfs(i,4)=(y(igp,jgp)-yin)/(ygam-yin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            end
           
        elseif (level(x(igp,jgp+1),y(igp,jgp+1))>0)            
            % look up for root
            fun = @(x)level(x,y(igp,jgp));
            xgam1 = fzero(fun,[x(igp,jgp) x(igp+1,jgp)]);
            ygam1 = y(igp,jgp);
            fun = @(y)level(x(igp,jgp),y);
            xgam2 = x(igp,jgp);
            ygam2 = fzero(fun,[y(igp,jgp) y(igp,jgp+1)]);
            if(abs(xgam1-x(igp,jgp)) < abs(ygam2-y(igp,jgp)))
                plot(xgam1,ygam1,'rs','markersize',14,'linewidth',2)
                xgam = xgam1; ygam = ygam1;
                ghost_points(i,3)=igp-1;
                ghost_points(i,4)=jgp;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(xgam-x(igp,jgp))/(xgam-xin);
                bcfs(i,4)=(x(igp,jgp)-xin)/(xgam-xin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            else
                plot(xgam2,ygam2,'gs','markersize',14,'linewidth',2)
                xgam = xgam2; ygam = ygam2;
                ghost_points(i,3)=igp;
                ghost_points(i,4)=jgp-1;
                bcfs(i,1)=xgam;
                bcfs(i,2)=ygam;
                xin=x(ghost_points(i,3),ghost_points(i,4));
                yin=y(ghost_points(i,3),ghost_points(i,4));
                bcfs(i,3)=(ygam-y(igp,jgp))/(ygam-yin);
                bcfs(i,4)=(y(igp,jgp)-yin)/(ygam-yin);
                plot(xin,yin,'g*','markersize',14,'linewidth',2)
            end
        else
            % Find root to the right! 
            fun = @(x)level(x,y(igp,jgp));
            xgam = fzero(fun,[x(igp,jgp) x(igp+1,jgp)]);
            ygam = y(igp,jgp);
            plot(xgam,ygam,'bp','markersize',14)
            ghost_points(i,3)=igp-1;
            ghost_points(i,4)=jgp;
            bcfs(i,1)=xgam;
            bcfs(i,2)=ygam;
            xin=x(ghost_points(i,3),ghost_points(i,4));
            yin=y(ghost_points(i,3),ghost_points(i,4));
            bcfs(i,3)=(xgam-x(igp,jgp))/(xgam-xin);
            bcfs(i,4)=(x(igp,jgp)-xin)/(xgam-xin);
            plot(xin,yin,'g*','markersize',14,'linewidth',2)
        end
    else
        if (level(x(igp,jgp-1),y(igp,jgp-1))>=0)
            % look down for root
            fun = @(y)level(x(igp,jgp),y);
            xgam = x(igp,jgp);
            ygam = fzero(fun,[y(igp,jgp-1) y(igp,jgp)]);
            plot(xgam,ygam,'gp','markersize',14,'linewidth',2)
            ghost_points(i,3)=igp;
            ghost_points(i,4)=jgp+1;
            bcfs(i,1)=xgam;
            bcfs(i,2)=ygam;
            xin=x(ghost_points(i,3),ghost_points(i,4));
            yin=y(ghost_points(i,3),ghost_points(i,4));
            bcfs(i,3)=(ygam-x(igp,jgp))/(ygam-yin);
            bcfs(i,4)=(y(igp,jgp)-yin)/(ygam-yin);
            plot(xin,yin,'g*','markersize',14,'linewidth',2)
        elseif (level(x(igp,jgp+1),y(igp,jgp+1))>=0)            
            % look up for root
            fun = @(y)level(x(igp,jgp),y);
            xgam = x(igp,jgp);
            ygam = fzero(fun,[y(igp,jgp) y(igp,jgp+1)]);
            plot(xgam,ygam,'gs','markersize',14,'linewidth',2)
            ghost_points(i,3)=igp;
            ghost_points(i,4)=jgp-1;
            bcfs(i,1)=xgam;
            bcfs(i,2)=ygam;
            xin=x(ghost_points(i,3),ghost_points(i,4));
            yin=y(ghost_points(i,3),ghost_points(i,4));
            bcfs(i,3)=(ygam-x(igp,jgp))/(ygam-yin);
            bcfs(i,4)=(y(igp,jgp)-yin)/(ygam-yin);
            plot(xin,yin,'g*','markersize',14,'linewidth',2)
        %else
        %    fun = @(y)level(x(igp,jgp),y);
        %    disp(level(x(igp,jgp-1),y(igp,jgp-1)))
         %   xgam = x(igp,jgp);
          %  ygam = fzero(fun,y(igp,jgp));
           % plot(xgam,ygam,'m*','markersize',14,'linewidth',2)
        end
    end        

   % if()
    %    bcfs(i,1)= xgam;
     %   bcfs(i,2) = ygp;
      %  bcfs(i,3)=(xgam-xgp)/(xgam-xsp);
       % bcfs(i,4)=(xgp-xsp)/(xgam-xsp);
    %else
    %end

end
return



return


%ghost_points
return

% TASK 1 mark the outermost point by changing maskin(i,j) = -1 if
% (i,j) is a boundary point. Check that you got this correct by
% addding some red stars to the above plot. DONE. 

% TASK 3 Set up a data structure with interpolation information for
% all the maskin = -1 points. 
% This can for example be one array of size nbp x 4 where the first
% two entries on each row correspond to the (i_igp,j_igp) coordinate of the
% inner ghostpoint and the second pair is the (i_in,j_in) coordinate of
% the interior point that is included in the interpolation (You
% need some logic to check that this point exists)   
% You also need a double array of size nbp x 4 with the values c_in and c_gam such
% that u_IGP = c_in*u(i_in,j_in) + c_gam*eb_boundary_fun(x_gam,y_gam,t); 
% you also need 2 values for x_gam y_gam.

% TASK 2 to do TASK 3 you first need a function that: given
% (i_igp,j_igp) and (i_in,j_in) finds x_gam and y_gam
% NOTE that one of these is trivial as it is the same as, say
% y_gam = y(i_igp,J) = y(i_in,J) when J = j_in = j_gam
% the other coordinate must then be found by root finding
% If level(z,y_gam) = 0, then x_gam = z;
% Note that you know that if, say, i_in < i_igp then 
% x(i_gp,J) <= z <= x(i_gp+1,J)


% TASK 4 Write a function that given the two arrays above updates
% the values at the interior ghost points at a given time t and
% assuming the function eb_boundary_fun(x,y,t) is known

% TASK 5 use twilight to impose "fake" exact boundary conditions in all
% interior ghostpoints (i.e. where maskin == -1) and make sure 
% the timesteping works. Once the boundary conditions have been
% updated you can still use the same lapu = compute_lap(u,x,y,t,0);
% followed by lapu = lapu.*mask; (this sets lapu = 0 on all outside
% points). The rest of the timestepping loop should not need to
% change. When you do the error computation make sure you only use
% mask == 1 points!




return


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
    lapu = lapu.*mask;
end

if (1==12)
    err = um-mmsfun(x,y,t-dt,0,0,0);
    mesh(x,y,err)
    title('Error in initial data for um')
    disp(['Max error in initial data for um is ' num2str(max(max(abs(err))))])
    pause
end

if (1==12)
    [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);
    u = update_bc(u,ulef,urig,ubot,utop);
    err = u-mmsfun(x,y,t,0,0,0);
    mesh(x,y,err)
    title('Error in BC')
    disp(['Max error in BC is ' num2str(max(max(abs(err))))])
    pause
end

if (1==12)

    lapu = compute_lap(u,x,y,t,0);
    err = lapu-mmsfun(x,y,t,0,2,0)-mmsfun(x,y,t,0,0,2);
    err(1,:) = 0;
    err(end,:) = 0;
    err(:,1) = 0;
    err(:,end) = 0;

    mesh(x,y,err)
    title('Error in Lap')
    disp(['Max error in Laplace ' num2str(max(max(abs(err))))])
    pause
end

if (1==12)
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
ERR = 0; 
TIME = 0;
for it = 1:nt
    t = (it-1)*dt;
    % Compute forcing on Boundaries 
    [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);
    % Update Boundary conditions;
    u = update_bc(u,ulef,urig,ubot,utop);
    % Impose B.C.s at ghost points here?
    % Compute Laplace(u)
    lapu = compute_lap(u,x,y,t,mms);
    lapu = lapu.*mask;
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
    err = u-mmsfun(x,y,t,0,0,0);
    title(['Error in ', num2str(it)  ,' timestep'])
    disp(['Max error in ' num2str(it) ' timestep ' num2str(max(max(abs(err))))])
    ERR = [ERR ; max(max(abs(err)))];
    TIME = [TIME ; t];
end


% Add option of error calculation at the final time.
err = u - mmsfun(x,y,Tend,0,0,0);
title('Error in first timestep')
disp(['Max error in final timestep ' num2str(max(max(abs(err))))])



function f = level(x,y);
    
    a = 0.82;
    b = 0.9;
    f = (x/a).^2+(y/b).^2 - 1;
    
end


function u = update_bc(u,ulef,urig,ubot,utop);
    
% Update Dirichlet boundary conditions for u 
    
    u(1,:) = ulef;
    u(end,:) = urig;
    u(:,1) = ubot;
    u(:,end) = utop;
end


function [ulef,urig,ubot,utop] = get_bc(x,y,t,mms);

% Compute Dirichlet boundary condition values for u

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

    [x_size, y_size] = size(U);
    hx = x(2,2) - x(1,1);
    hy = y(2,2) - y(1,1);

    lapu = zeros(size(U));
    if(mms ==1)
        lapu = mmsfun(x,y,t,0,2,0)+mmsfun(x,y,t,0,0,2);
    else

        for ix = 2:x_size-1
            for jy = 2:y_size-1
                lapu(ix,jy) = 1/hx^2 *(U(ix+1,jy)-2*U(ix,jy)+U(ix-1,jy)) + 1/hy^2*(U(ix,jy+1)-2*U(ix,jy)+U(ix,jy-1));
            end
        end

    end
end

function f = compute_forcing(x,y,t,mms);

% Compute forcing
    
    if(mms ==1)
        % f = u_tt - u_xx - u_yy
        f = mmsfun(x,y,t,2,0,0)-mmsfun(x,y,t,0,2,0)-mmsfun(x,y,t,0,0,2);
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

function u = update_IGP(u,array1,array2,t);
    % Inputs:
        % u
        % array1: array of size (nbp,4) where each row is
        % (i_igp,j_igp,i_in,j_in)
        % array2: array of size (nbp,4) where each row is
        % (c_in,c_gam,x_gam,y_gam)
        % t: time
    
    % Output: 
        % u: u with updated "boundary conditions" at interior ghost points
        % at time t
        
    nbp = size(array1,1); % number of interior ghost points
    u_IGP = u;
    
    for i = 1:nbp
        % check all these indices based on actual arrays from Andres/Yutao!
        i_igp = array1(i,1);
        j_igp = array1(i,2);
        i_in  = array1(i,3);
        j_in  = array1(i,4);
        
        c_in  = array2(i,1);
        c_gam = array2(i,2);
        x_gam = array2(i,3);
        y_gam = array2(i,4);
        
        u_IGP(i_igp,j_igp) = c_in*u(i_in,j_in) + c_gam*eb_boundary_fun(x_gam,y_gam,t);
        % not sure if the i_igp and j_igp are coordinates or indices -
        % check once Andres/Yutao part updated
        
    end
    
    u = u_IGP;
    
end
