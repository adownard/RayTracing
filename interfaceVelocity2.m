%%%%%% Code to Generate Nice Looking Plots %%%%%%
clear
close all
%%%%%%%% Create Velocity and Slowness Models %%%%%%%%
v0 = 0.5;
vf = 5;
dv = .01; 
nv = floor((vf - v0)/dv);

nx = nv; 
nz = nv;
vals = 1:nv;
vel=1:nv; 
vel(1:100)=1; 
vel(101:240)=8; 
vel(241:379)=.5; 
vel(379:450)=4; 
slo=1./vel; 

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
xPlt=vals;
zPlt=vals; 

xPos(1) = xPos0;
zPos(1) = zPos0;

a0 = 0;
af = 90;
da = 4; 
na = floor((af - a0)/da);

gsz=vals; 
%dsz=(zsloVec(end)-zsloVec(1))/nv; 
 
for a=1:na+1
    alpha = a0 + (a-1)*da;
    xDir0 = sind(alpha); 
    zDir0 = cosd(alpha); 
    xDir(1) = xDir0/vel(1); 
    zDir(1) = zDir0/vel(1); 

    for ix=1:nv-1
        xPos(ix+1) = xPos(ix) + xDir(ix)*dt;
        xDir(ix+1) = xDir(ix) + 0*dt; 
    end
    for ix=2:nv-1
        if ix == 0
            gsz(ix) = (slo(ix+1) - slo(ix))/dt;
        elseif ix == nv
            gsz(ix) = (slo(ix) - slo(ix-1))/dt;
        else
            gsz(ix) = (slo(ix+1) - slo(ix-1))/(2*dt);
        end
        zPos(ix+1) = zPos(ix) + zDir(ix)*dt;
        zDir(ix+1) = zDir(ix) + gsz(ix)*dt; 
        if zPos(ix+1) < 0
            zPos(ix+1) = 0; 
        end
    end
    allPosVec(a).x_Pos = xPos; 
    allPosVec(a).z_Pos = zPos; 
end


%%%%%%%%%%%%% Analytic Solution %%%%%%%%%%%%%%%%%
gVz=1:nv;
g1z=vel(100)/101; 
g2z=vel(240)/140; 
g3z=vel(379)/138; 
g4z=vel(450)/71; 
gVz(1:100)=g1z; 
gVz(101:240)=g2z; 
gVz(241:379)=g3z; 
gVz(380:450)=g4z;
gVx=0; 
tVec=1:nt; 

for it=1:nt               % fill tVec with times
    tVec(it) = t0 + (it-1) * dt;
end
for a=1:na+1
    alpha = a0 + (a-1)*da;
    for it=1:100
        xDir0 = sind(alpha)/vel(it); 
        zDir0 = cosd(alpha)/vel(it);
        t = tVec(it);
        xPos_ = xPos0 + xDir0*t + gVx * t*t * 0.5;
        xPlt(it) = xPos_; 
        zPos_ = zPos0 + zDir0*t + gVz(it) * t*t * 0.5; 
        zPlt(it) = zPos_;
    end
    for it=101:240
        xDir0 = sind(alpha)/vel(it); 
        zDir0 = cosd(alpha)/vel(it); 
        xPos0=xPlt(100); 
        zPos0=zPlt(100); 
        t = tVec(it);
        xPos_ = xPos0 + xDir0*t + gVx * t*t * 0.5;
        xPlt(it) = xPos_; 
        
        zPos_ = zPos0 + zDir0*t + gVz(it) * t*t * 0.5; 
        zPlt(it) = zPos_;
    end
    for it=241:379
        xDir0 = sind(alpha)/vel(it); 
        zDir0 = cosd(alpha)/vel(it); 
        xPos0=xPlt(240); 
        zPos0=zPlt(240); 
        t = tVec(it);
        xPos_ = xPos0 + xDir0*t + gVx * t*t * 0.5;
        xPlt(it) = xPos_; 
        zPos_ = zPos0 + zDir0*t + gVz(it) * t*t * 0.5; 
        zPlt(it) = zPos_;
    end
    for it=380:450
        xDir0 = sind(alpha)/vel(it); 
        zDir0 = cosd(alpha)/vel(it); 
        xPos0=xPlt(379); 
        zPos0=zPlt(379); 
        t = tVec(it);
        xPos_ = xPos0 + xDir0*t + gVx * t*t * 0.5;
        xPlt(it) = xPos_; 
        zPos_ = zPos0 + zDir0*t + gVz(it) * t*t * 0.5; 
        zPlt(it) = zPos_;
    end
    allPltVec(a).x_Plt = xPlt; 
    allPltVec(a).z_Plt = zPlt; 
    %plot(xPlt,zPlt);        
end
    
set(gca,'YDir','reverse'); 
%%%%% Plot Velocity Colormap in the Background %%%%%%%%
maxx=0; 
maxz=0;

for ia=1:na
    tempx=max(allPosVec(ia).x_Pos); 
    tempz=max(allPosVec(ia).z_Pos); 
    if tempx > maxx
        maxx=tempx; 
    end
    if tempz > maxz
        maxz=tempz; 
    end
end

scalex=nv/maxx;
scalez=nv/maxz;
colormap gray(4); 
slo2D=repmat(slo,nv,1);
vel2D=transpose(repmat(vel,nv,1)); 
t_slo=transpose(slo2D); 
hold on 
for ia=1:na
    aPosx=allPosVec(ia).x_Pos * scalex; 
    aPosz=allPosVec(ia).z_Pos * scalez;
    plot(aPosx,aPosz,'LineWidth',1.35)
end
title('Euler''s Method in an Anisotropic Veclocity Medium','fontsize',18)
ylabel('Depth of Ray','fontsize',14)
xlabel('Horizontal Distance Traveled','fontsize',14)
axis([1,nv,0,nv])
b=surf(vel2D,'EdgeColor','none');
view(0,90);
z=get(b,'Zdata');
set(b,'Zdata',z-10);
set(gca,'YDir','reverse');
c = colorbar;
c.Label.String= 'Slowness';
c.FontSize=12;
hold off