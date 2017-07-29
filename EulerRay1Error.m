%%%%%% Euler Method and Analytical Solution Ray Tracing %%%%%%

clear
close all

%%%%%%%% Create Velocity and Slowness Models %%%%%%%%

v0 = 0.5;
vf = 5;
dv = .01; 
nv = floor((vf - v0)/dv);

s0_2 = 1/v0^2;
sf_2 = 1/vf^2;

gradSlo2 = (sf_2 - s0_2) / nv;

nx = nv; 
nz = nv;

            
%%%%%%%% Euler's Method Implementation %%%%%%%%%

t0 = 0; 
tf = 200;
%nt = nv; 

xPos0 = 0; 
zPos0 = 0; 
alpha=30; 
gVx = 0;
gVz = gradSlo2; 
xDir0 = sind(alpha); 
zDir0 = cosd(alpha); 


%%%%%%% Analytical Solution Parameters
dtVec=1:10; 



for itf=1:10
    dt=(0.5)^itf;
    nt = (tf - t0)/dt;
    vals=1:nt;
    xPos = vals;
    zPos = vals;
    xPlt = vals;
    zPlt = vals;
    xDir = vals;
    zDir = vals;
    tVec = vals;
    slope = vals;
    d_t = vals;
    error = vals; 
    
    d_t(itf) = dt;
    
    xDir(1) = xDir0; 
    zDir(1) = zDir0; 
    xPos(1) = xPos0;
    zPos(1) = zPos0;
    
    for ix=1:nt-1
        xPos(ix+1) = xPos(ix) + xDir(ix)*dt;
        xDir(ix+1) = xDir(ix) + gVx*dt; 

        zPos(ix+1) = zPos(ix) + zDir(ix)*dt;
        zDir(ix+1) = zDir(ix) + gVz*dt; 

        if zPos(ix+1) < 0
            zPos(ix+1) = 0; 
        end
    end
    %%%%%%%%%%% Analytical Solution %%%%%%%%%%%%%%
    for it=1:nt               % fill tVec with times
        tVec(it) = t0 + (it-1) * dt;
    end
    for it=1:nt
        t = tVec(it);
        xPos_ = xPos0 + xDir0*t + gVx * t*t * 0.5;
        zPos_ = zPos0 + zDir0*t + gVz * t*t * 0.5;

        if zPos_ < 0
            zPos_ = 0;
        end

        xPlt(it) = xPos_; 
        zPlt(it) = zPos_; 
    end
    diffX = xPos(4)-xPlt(4); 
    diffZ = zPos-zPlt; 
    magDiff = sqrt(diffX.^2 + diffZ.^2);
    maxmD = max(magDiff);
    slope(itf) = maxmD/find(magDiff==maxmD);
    
    %%%%% for comparison %%%%%
    error(itf) = -gradSlo2 * dt*dt * 0.5; 
    
end

subplot(3,1,1)
plot(d_t,slope,'LineWidth',1.25)
title('Error Values for Different Time Steps','fontsize',16)
xlabel('Time Step (dt)')
ylabel('Error')
subplot(3,1,2)
plot(d_t,error,'LineWidth',1.25)
title('Expected Error for First Order Scheme','fontsize',16)
xlabel('Time Step (dt)')
h1=ylabel('Error (g_{0}dt^{2}/2)');
subplot(3,1,3)
plot(log(d_t),log(slope),'LineWidth',1.25)
title('Natural Log of Error','fontsize',16)
xlabel('Natural Log of dt')
ylabel('Natural Log of Error')
hold off
set([h1], 'interpreter', 'tex')




