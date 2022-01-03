%Compact 2D magnetic vertical componet modelling (case : fluxgate magnetometer)
%Mohammad Rheza Zamani
%Reference : Stocco, Stefano & Godio, A. & Sambuelli, Luigi. (2009). Modelling and compact inversion of magnetic data: A Matlab code. Computers & Geosciences. 35. 2111-2118. 10.1016/j.cageo.2009.04.002. 
clear all;
clc;
Fe = 40000; %Eartch magnetic field intensity (nT) (25,000 nT - 65,000 nT)
betha = 10; %strike angle of the prism relative to themagnetic north in degree
I =  30; % the earth's magnetic field inclination in degree
lengthx = 1000;
lengthz = 1000;
%Block dimenssion
dx = 20;
dz = 20;
%Block middle point 
dmx = dx/2;
dmz = dz/2;
%Number of vertical and horizontal block
nx = lengthx/dx;
nz = lengthz/dz;
%Total number of block
nb = nx*nz;
%Suscepbility contrast
V = zeros(nz,nx);
%V(24:25,1:12) = 0.9*10^-3;
%V(26:27,12:22) = 1*10^-3;
%V(28:29,22:32) = 0.7*10^-3;
%V(30:31,32:42) = 0.6*10^-3;
%V(32:33,42:50) = 0.5*10^-3;
V(1:12,24:25) = 1*10^-3;
V(12:22,26:27) = 0.8*10^-3;
%Make block model
for i = 1 : nx
    x(i) = dx*i - dmx;
end
xx = repmat(x,nz,1);
for j = 1 : nz
    z(j) = dz*j - dmz;
end
z1 = z';
zz = repmat(z1,1,nx);

%Kernel matrix
for i=1:nx
    for j = 1:nb
        r1 = sqrt((zz(j)-dz/2).^2 + (x(i)-xx(j)+dx/2).^2);
        r2 = sqrt((zz(j)+dz/2).^2 + (x(i)-xx(j)+dx/2).^2);
        r3 = sqrt((zz(j)-dz/2).^2 + (x(i)-xx(j)-dx/2).^2);
        r4 = sqrt((zz(j)+dz/2).^2 + (x(i)-xx(j)-dx/2).^2);
        theta1 = atan((x(i)-xx(j)+dx/2)/(zz(j)-dz/2));
        theta2 = atan((x(i)-xx(j)+dx/2)/(zz(j)+dz/2));
        theta3 = atan((x(i)-xx(j)-dx/2)/(zz(j)-dz/2));
        theta4 = atan((x(i)-xx(j)-dx/2)/(zz(j)+dz/2));    
        Kernell(i,j) = 2.*(((cosd(I))*(sind(betha))*log((r2*r3)/(r4*r1))-sind(I.*(theta1-theta2-theta3+theta4))));
    end 
end
%Calculated magnetic response
V_rs = reshape(V,nb,1);
dBdz = Fe.*Kernell*V_rs;

%Plot kernel matrix
figure(1)
imagesc(Kernell)
set(gcf, 'Position', get(0, 'Screensize'));
ylabel('Observation Points','FontWeight','bold','FontSize',10)
xlabel('The Blocks','FontWeight','bold','FontSize',10)
cb = colorbar;
cb.Label.String = 'Value';
cb.Location = 'southoutside';
colormap(jet)

figure(2)
subplot(2,1,1)
plot(x,dBdz,'*-b')
xlabel('Distance (m)','FontWeight','bold','FontSize',10)
ylabel('Magnetic Anomaly(nT)','FontWeight','bold','FontSize',10)
title('Geomagnetic Response','FontWeight','bold','FontSize',10)
grid on
subplot(2,1,2)
s = pcolor(x,z,V);
s.FaceColor = 'interp';
xlabel('Distance (m)','FontWeight','bold','FontSize',10)
ylabel('Depth(m)','FontWeight','bold','FontSize',10)
title('Subsurface Model','FontWeight','bold','FontSize',10)
cb = colorbar;
cb.Label.String = 'Contrast Susceptibility (SI) ';
cb.Location = 'southoutside';
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'Ydir','reverse')
colormap(jet)

data = [x' dBdz];
saveas(figure(1),'Kernell compact magnetic.png')
saveas(figure(2),'Model compact Magnetic.png')
writematrix(data,'Data forward compact magnetic.dat','Delimiter','tab')
