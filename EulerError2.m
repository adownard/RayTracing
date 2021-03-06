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
vals = [1:nv];
            


t0 = 0; 
tf = 100:10:250;
nt = nv; 
slope = 1:numel(tf);
d_t = 1:numel(tf);
error = 1:numel(tf); 
xPos0 = 0; 
zPos0 = 0; 

xPos = vals;
zPos = vals;
xDir = vals;
zDir = vals;
tVec = vals;

xPos(1) = xPos0;
zPos(1) = zPos0;

alpha=30; 

gVx = 0;
gVz = gradSlo2; 
xDir0 = sind(alpha); 
zDir0 = cosd(alpha); 
xDir(1) = xDir0; 
zDir(1) = zDir0; 

%%%%%%% Analytical Solution Parameters

xPlt = [1:nx];
zPlt = [1:nz];

for itf=1:numel(tf)
    dt = (tf(itf) - t0)/nt;
    d_t(itf) = dt; 
    %%%%%%%% Euler's Method Implementation %%%%%%%%%
    for ix=1:nx-1
        xPos(ix+1) = xPos(ix) + xDir(ix)*dt + gVx * dt*dt * 0.5;
        xDir(ix+1) = xDir(ix) + gVx * dt;                          %% since gVx is constant only need gVx not 2*gVx/2 
   
        zPos(ix+1) = zPos(ix) + zDir(ix)*dt + gVz * dt*dt * 0.5;
        zDir(ix+1) = zDir(ix) + gVz * dt; 
       
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
    diffX = xPos-xPlt; 
    diffZ = zPos-zPlt; 
    magDiff = sqrt(diffX.^2 + diffZ.^2);
    maxmD = max(magDiff);
    slope(itf) = maxmD/find(magDiff==maxmD,1);
    
end
scatter(d_t,slope)
title('Error Values for Different Time Steps')
xlabel('Time Step (dt)')
ylabel('Error')
hold off