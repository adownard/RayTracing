% %%%%%% Euler Method and Analytical Solution Ray Tracing %%%%%%
% 
% clear
% close all
% 
% %%%%%%%% Create Velocity and Slowness Models %%%%%%%%
% 
% v0 = 0.5;
% vf = 5;
% dv = .01; 
% nv = floor((vf - v0)/dv);
% 
% s0_2 = 1/v0^2;
% sf_2 = 1/vf^2;
% 
% gradSlo2 = (sf_2 - s0_2) / nv;
% 
% % slo2 = [1:nv];
% % 
% % for is=1:nv
% %     slo2(is) = s0_2 + is*gradSlo2; 
% % end
% % 
% % slo = sqrt(slo2); 
% 
% nx = nv; 
% nz = nv;
% vals = [1:nv];
% 
% %%%%%%%% Structs for Evaluating the Difference Between Solutions %%%%%%%%
% allPosVec = struct('x_Pos',vals,'z_Pos',vals);
% allPltVec = struct('x_Plt',vals,'z_Plt',vals);
%             
% %%%%%%%% Euler's Method Implementation %%%%%%%%%
% 
% t0 = 0; 
% tf = 200;
% nt = nv; 
% dt = (tf-t0)/nt;
% 
% xPos0 = 0; 
% zPos0 = 0; 
% 
% xPos = vals;
% zPos = vals;
% xDir = vals;
% zDir = vals;
% tVec = vals;
% 
% xPos(1) = xPos0;
% zPos(1) = zPos0;
% 
% a0 = 0;
% af = 90;
% da = 4; 
% na = floor((af - a0)/da);
% 
% gVx = 0;
% gVz = gradSlo2; 
% 
% 
% hold on 
% for a=1:na+1
%     alpha = a0 + (a-1)*da;
%     xDir0 = sind(alpha); 
%     zDir0 = cosd(alpha); 
%     xDir(1) = xDir0; 
%     zDir(1) = zDir0; 
% 
%     for ix=1:nx-1
%         xPos(ix+1) = xPos(ix) + xDir(ix)*dt;
%         xDir(ix+1) = xDir(ix) + gVx*dt; 
% 
%         zPos(ix+1) = zPos(ix) + zDir(ix)*dt;
%         zDir(ix+1) = zDir(ix) + gVz*dt; 
% 
%         if zPos(ix+1) < 0
%             zPos(ix+1) = 0; 
%         end
% 
%     end
%     allPosVec(a).x_Pos = xPos; 
%     allPosVec(a).z_Pos = zPos; 
% 
%     plot(xPos,zPos)
% 
% end
% title('Euler''s Method Ray Tracing','fontsize',16)
% ylabel('Depth of Ray')
% xlabel('Horizontal Distance Traveled')
% %axis([0,120,0,50])
% axis ij
% hold off
% 
% %%%%%%%%%%% Analytical Solution %%%%%%%%%%%%%%
%  figure 
% 
% for it=1:nt               % fill tVec with times
%     tVec(it) = t0 + (it-1) * dt;
% end
% 
% xPlt = [1:nx];
% zPlt = [1:nz];
% 
% hold on 
% for a=1:na+1
%     alpha = a0 + (a-1)*da;
%     xDir0 = sind(alpha); 
%     zDir0 = cosd(alpha); 
% 
%     for it=1:nt
%         t = tVec(it);
%         xPos = xPos0 + xDir0*t + gVx * t*t * 0.5;
%         zPos = zPos0 + zDir0*t + gVz * t*t * 0.5;
% 
%         if zPos < 0
%             zPos = 0;
%         end
% 
%         xPlt(it) = xPos; 
%         zPlt(it) = zPos; 
%     end
% 
%     allPltVec(a).x_Plt = xPlt; 
%     allPltVec(a).z_Plt = zPlt; 
% 
%     plot(xPlt,zPlt)
% end
% 
%     title('Analytical Solution','fontsize',16)
%     ylabel('Depth of Ray')
%     xlabel('Horizontal Distance Traveled')
%     %axis([0,120,0,50])
%     axis ij
%     %colorbar(1/slo)
%     hold off
% 
% diffX = struct('difference',[1:nv]); 
% diffZ = struct('difference',[1:nv]); 
% magDiff = struct('magDif',[1:nv]);
% for ia=1:na
%     diffX(ia).difference = allPosVec(ia).x_Pos - allPltVec(ia).x_Plt; 
%     diffZ(ia).difference = allPosVec(ia).z_Pos - allPltVec(ia).z_Plt; 
%     magDiff(ia).magDif = sqrt((diffX(ia).difference).^2 + (diffZ(ia).difference).^2);
%     end
% figure 
% plot(magDiff(1).magDif)
% title('Magnitude of the Difference Between Rays')



%%%%%% First Order Error Analysis %%%%%%
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
            
%%%%%%%% Euler's Method Implementation %%%%%%%%%

t0 = 0; 
tf = 100:10:200;
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
    for ix=1:nx-1
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
a=log(d_t); 
b=log(slope); 
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
plot(a,b,'LineWidth',1.25)
title('Natural Log of Error','fontsize',16)
xlabel('Natural Log of dt')
ylabel('Natural Log of Error')
axis([a(1),a(end),b(1),b(end)])
hold off
set([h1], 'interpreter', 'tex')