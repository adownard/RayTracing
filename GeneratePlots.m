%%%%%% Code to Generate Nice Looking Plots %%%%%%

clear
close all

%%%%%%%% Create Velocity and Slowness Models %%%%%%%%

v0 = 0.5;
vf = 5;
dv = .005; 
nv = floor((vf - v0)/dv);

s0_2 = 1/v0^2;
sf_2 = 1/vf^2;

gradSlo2 = (sf_2 - s0_2) / nv;

nx = nv; 
nz = nv;
vals = [1:nv];

%%%%%%%% Structs for Evaluating the Difference Between Solutions %%%%%%%%
allPosVec = struct('x_Pos',vals,'z_Pos',vals);
allPltVec = struct('x_Plt',vals,'z_Plt',vals);

%%%%%%%% Euler's Method Implementation %%%%%%%%%

t0 = 0; 
tf = 430;
nt = nv; 
dt = (tf-t0)/nt;

xPos0 = 0; 
zPos0 = 0; 

xPos = vals;
zPos = vals;
xDir = vals;
zDir = vals;
tVec = vals;

xPos(1) = xPos0;
zPos(1) = zPos0;

a0 = 0;
af = 90;
da = 4; 
na = floor((af - a0)/da);

gVx = 0;
gVz = gradSlo2; 

hold on 
for a=1:na+1
    alpha = a0 + (a-1)*da;
    xDir0 = sind(alpha); 
    zDir0 = cosd(alpha); 
    xDir(1) = xDir0; 
    zDir(1) = zDir0; 

    for ix=1:nx-1
        xPos(ix+1) = xPos(ix) + xDir(ix)*dt;
        xDir(ix+1) = xDir(ix) + gVx*dt; 
    end
    for ix=1:nx-1
        zPos(ix+1) = zPos(ix) + zDir(ix)*dt;
        zDir(ix+1) = zDir(ix) + gVz*dt; 
        if zPos(ix+1) < 0
            zPos(ix+1) = -1; 
        end
        if zPos(ix+1) > 100
            zPos(ix+1:nx) = 101;
            break
        end
    end
    allPosVec(a).x_Pos = xPos; 
    allPosVec(a).z_Pos = zPos; 
end

%%%%%%%%%%% Analytical Solution %%%%%%%%%%%%%%

for it=1:nt               % fill tVec with times
    tVec(it) = t0 + (it-1) * dt;
end

xPlt = [1:nx];
zPlt = [1:nz];

hold on 
for a=1:na+1
    alpha = a0 + (a-1)*da;
    xDir0 = sind(alpha); 
    zDir0 = cosd(alpha); 

    for it=1:nt
        t = tVec(it);
        xPos_ = xPos0 + xDir0*t + gVx * t*t * 0.5;
        xPlt(it) = xPos_; 
    end
    for it=1:nt
        t = tVec(it); 
        zPos_ = zPos0 + zDir0*t + gVz * t*t * 0.5; 
        if zPos_ < 0
            zPos_ = -1;
        end
        if zPos_ > 100
            zPlt(it:nt) = 101;
            break
        else
            zPlt(it) = zPos_;
        end 
    end
    allPltVec(a).x_Plt = xPlt; 
    allPltVec(a).z_Plt = zPlt; 
end


%%%%%% Plot Velocity Colormap in the Background %%%%%%%%
maxx=0; 
maxz=0;
for ia=1:na
    tempx=max(allPltVec(ia).x_Plt); 
    tempz=max(allPltVec(ia).z_Plt); 
    if tempx > maxx
        maxx=tempx; 
    end
    if tempz > maxz
        maxz=tempz; 
    end
end

scalex=floor(maxx)+1;
scalez=floor(maxz)+1;
map=colormap(gray((maxx+1)*(maxz+1))); 
slo2 = [1:maxz];

for is=1:maxz
    slo2(is) = s0_2 + is*gradSlo2; 
end

slo = sqrt(slo2); 
vel=1./slo; 

slo_=repmat(slo,nv,1);
t_slo=transpose(slo_);

for ia=1:na
    plot(allPltVec(ia).x_Plt,allPltVec(ia).z_Plt,'LineWidth',1.35)
end
title('Analytical Solution','fontsize',18)
ylabel('Depth of Ray','fontsize',14)
xlabel('Horizontal Distance Traveled','fontsize',14)
axis([.25,230,0,maxz-1])
b=surf(t_slo,'EdgeColor','none');
hold on 
z=get(b,'Zdata');
set(b,'Zdata',z-10);
set(gca,'YDir','reverse');
c = colorbar;
c.Label.String= 'Slowness'
c.FontSize=12
hold off

figure;
hold on 
map=colormap(gray((maxx+1)*(maxz+1))); 
for ia=1:na
    plot(allPosVec(ia).x_Pos,allPosVec(ia).z_Pos,'LineWidth',1.35)
end
title('Euler''s Method Numerical Solution','fontsize',18)
ylabel('Depth of Ray','fontsize',14)
xlabel('Horizontal Distance Traveled','fontsize',14)
axis([.25,230,0,maxz-1])
b=surf(t_slo,'EdgeColor','none');
z=get(b,'Zdata');
set(b,'Zdata',z-10);
set(gca,'YDir','reverse');
c = colorbar;
c.Label.String= 'Slowness'
c.FontSize=12
hold off


