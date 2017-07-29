%%%%%% Code to Generate Nice Looking Plots %%%%%%

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

slo2Vec=vals;
for islo=1:nv
    slo2Vec(islo)=s0_2+gradSlo2*islo; 
end
slo=sqrt(slo2Vec); 
zv0=.7; 
zvf=7; 
zslo0=1/zv0; 
zslof=1/zvf;
gsz=(zslof-zslo0)/nv; 
xsloVec=vals; 
zsloVec=vals; 

for ixz=1:nv
    zsloVec(ixz)=zslo0+ixz*gsz; 
    xsloVec(ixz)=sqrt(slo(ixz)^2-zsloVec(ixz)^2);
end

gsx=vals; 
dsx=(xsloVec(end)-xsloVec(1))/nv; 
for a=1:na+1
    alpha = a0 + (a-1)*da;
    xDir0 = sind(alpha); 
    zDir0 = cosd(alpha); 
    xDir(1) = xDir0; 
    zDir(1) = zDir0; 
    
    for ix=1:nv-1
        if ix == 1
            gsx(ix) = (xsloVec(ix+1) - xsloVec(ix))/dsx; 
        else
            gsx(ix) = (xsloVec(ix+1) - xsloVec(ix-1))/(2*dsx);
        end
        gsx(nv)=gsx(nv-1); 
        xPos(ix+1) = xPos(ix) + xDir(ix)*dt;
        xDir(ix+1) = xDir(ix) + gsx(ix)*dt; 
    end
    for ix=1:nv-1
        zPos(ix+1) = zPos(ix) + zDir(ix)*dt;
        zDir(ix+1) = zDir(ix) + gsz*dt; 
        if zPos(ix+1) < 0
            zPos(ix+1) = -1; 
        end
    end
    allPosVec(a).x_Pos = xPos; 
    allPosVec(a).z_Pos = zPos; 
end


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
map=colormap(gray((maxx+1)*(maxz+1))); 
slo2 = [1:maxz];

slo2D=zeros(nv,nv);
for ix=1:nv
    for iz=1:nv
        slo2D(ix,iz)=(xsloVec(ix)^2+zsloVec(iz)^2);
    end
end
vel=1./slo2D;

hold on
for ia=1:na
    aPosx=allPosVec(ia).x_Pos * scalex; 
    aPosz=allPosVec(ia).z_Pos * scalez;
    plot(aPosx,aPosz,'LineWidth',1.35)
end
title('Euler''s Method in an Anisotropic Veclocity Medium','fontsize',18)
ylabel('Depth of Ray','fontsize',14)
xlabel('Horizontal Distance Traveled','fontsize',14)
axis([1,nv,0,nv]);
ax=gca;
ax.Box='on';
b=surf(slo2D,'EdgeColor','none');
view(0,90);
z=get(b,'Zdata');
set(b,'Zdata',z-10);
set(gca,'YDir','reverse');
c = colorbar;
c.Label.String= 'Slowness';
c.FontSize=12;
hold off

%%%%%%%%% Error %%%%%%%%%
error=dt*dt*0.5*sqrt(gsz^2 + gsx.^2); 
figure; 
plot(error)
xlabel('Number of Steps (dt)')
ylabel('Magnitude of Error')
title('Error for Inhomogenous Velocity Model')






