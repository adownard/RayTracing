% Ray tracing written by Luke Decker
clear

dt = .01;
nt = 1000;




t0 = 0;

tf = t0+nt*dt;
tVec = [1:nt+1];
x0 = [0 0];

%%%%% Time Vector for analytic solution %%%%%
for it=0:nt
    tVec(it+1) = t0 + dt*(it);
end

% velocity field

v0 = 0.5;
vf = 1;

vz0 = 0;
vzf = 1;
nvz = 101;


s0 = 1/v0;
s02 = s0*s0;
sf = 1/vf;
sf2 = sf*sf;

dels = sf2-s02;
ds = dels/(nvz-1);
dvz = (vzf-vz0)/(nvz-1);

slow2=zeros(nvz,1);



for i=1:nvz
   slow2(i) =  s02 + ds*(i-1);  
end

vel = 1./sqrt(slow2);

gX = 0;
gZ = dels;

initialp = 0;
np = 81 ; 
dp = 2;



ymin=0;
ymax=0;
hold on

for ip=1:np
    
    pindex = (ip-1)*(nt+1);
    
    theta0 = initialp+dp*(ip-1);
    p0 = [cosd(theta0) sind(theta0)];
    X = zeros(nt+1,2);
    P = X;
    X(1,:) = x0;

    P(1,:) = p0;
    
    for it=2:nt+1
        
        
        
        z = X(it-1,1);
        pz = P(it-1,1);

        x = X(it-1,2);
        px = P(it-1,2);

        iz = floor((z-vz0)/dvz)+1;
        if iz<1.5
            iz = 2;
        end
        
        if iz>nvz-1
            iz = nvz-1;
        end
        dXt = x * px * dt + px * gX * dt*dt + gX*gX * dt*dt*dt * (1./3);
        dZt = z * pz * dt + pz * gZ * dt*dt + gZ*gZ * dt*dt*dt * (1./3);

        z1 = z + pz*vel(iz)*dt;
        x1 = x + px*vel(iz)*dt;

        gradvZ = (vel(iz+1)-vel(iz-1))/(dvz*2);
        deriv = -gradvZ/vel(iz);

        pz10 = pz + deriv*dt; % derivitave of p in packet

        % normalize
        amp = sqrt(pz10*pz10+px*px); % normalization constant for ray parameter


        X(it,1) = z1; % write out new position
        X(it,2) = x1;

        P(it,1) = pz10/amp; % normalized ray parameter
        P(it,2) = px/amp;
        

        if z1<vz0 % if breach surface, write out rest at same position
            for it2 = it+1:nt+1
                X(it2,1) = z1;
                X(it2,2) = x1;

                P(it2,1) = pz10/amp;
                P(it2,2) = px/amp;

            end
            break
        end
        if z1>vzf % if beneath velocity file, write out position
            for it2 = it+1:nt+1
                X(it2,1) = z1;
                X(it2,2) = x1;

                P(it2,1) = pz10/amp;
                P(it2,2) = px/amp;

            end
            break
        end
    end
    if max(X(:,2)) > ymax
       ymax =  max(X(:,2));
    end
    if min(X(:,2)) < ymin
       ymin =  min(X(:,2));
    end
    
    plot(X(:,2),X(:,1))
end


title('Explicit Euler Ray Tracing')
ylabel('depth (km)')
xlabel('horizontal (km)')
axis([ymin,ymax,vz0,vzf])
axis ij
pbaspect([(ymax-ymin) 1 1])
hold off

figure

%%%%%%%%%%%% Analytic Solution %%%%%%%%%%%%%%
xPlt = [1:nt+1];
zPlt = [1:nt+1];



hold on;
for ip= 1:np
    p = initialp + dp*ip;
    xDir0 = sind(p);
    zDir0 = cosd(p);
  
    for i=1:nt+1
        t = tVec(i); 
        xPos = x0(2) + xDir0*t + gX * t*t * 0.5;
        zPos = x0(1) + zDir0*t + gZ * t*t * 0.5;
        if zPos < 0
            zPos = 0;
        end
        
        xPlt(i) = xPos;
        zPlt(i) = zPos;
    end

    plot(xPlt,zPlt)
        
end
title('Analytic Solution Ray Tracing')
ylabel('depth (km)')
xlabel('horizontal (km)')
axis([ymin,ymax,vz0,vzf])
axis ij
pbaspect([(ymax-ymin) 1 1])
hold off
        
    






















