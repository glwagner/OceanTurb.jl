dz=2;
depth=80;
tmax=45;
lat=31;

q0=-300; 
tau0=0.1;

z=[-dz/2:-dz:-depth+dz/2]';
lz=length(z);
zmid=[z(1:lz-1)-dz/2]';
%temp=4+11*exp(z/500)+10*exp(z/10);
temp=15+10*exp(z/10);
u=0*z;
v=0*z;

as=[14 26 -depth 0];

rho0=1000;
cp=3900;
alpha=2e-4;
g=9.81;

dt=1/32;

f=2*2*pi*sin(lat/180*pi);
sa=sin(f*dt/2);
ca=cos(f*dt/2);

%us=u(1);
%vs=v(1);
fac1=0.25/dz/dz;
fac2=alpha*g/dz;

frm=-1;
plot(temp,z,u*10+20,z,v*5+10,z)
legend("T","u","v")
axis(as);
title(sprintf('t = %f',0));
drawnow;
frm=mkpng(frm);

ts=temp;

for t=dt:dt:tmax
 q=q0;
 if t>15 & t<30
  q=-q;
  end
 if rem(t,1)==0 && t<=18
  dc=1;
%  purge_tmp_files;
 else
  dc=0;
 end
 temp(1)=temp(1)+q*dt*86400/rho0/cp/dz;
 tav=temp(1);
 uav=u(1);
 vav=v(1);
 ind=2;
 while tav<temp(ind)
  tav=(tav*(ind-1)+temp(ind))/ind;
  uav=(uav*(ind-1)+u(ind))/ind;
  vav=(vav*(ind-1)+v(ind))/ind;
  ind=ind+1;
  if ind>lz
    break;
    end
  end
 if dc & ind>2
  plot(temp,z,tav,(1-ind)*dz,'o');
  title(sprintf('Convective mixing step -- t = %f',t));
  axis(as);
  drawnow;
  frm=mkpng(frm);
%  pause;
 end
 for j=1:ind-1
  temp(j)=tav;
  u(j)=uav;
  v(j)=vav;
  end
 
  mld=(ind-1)*dz;

  tau=tau0*24*3600/rho0/mld*dt;
  ut=u*ca+v*sa;
  v=v*ca-u*sa;
  u=ut;
  du=tau*(z > -mld);
  u=u+du;
  ut=u*ca+v*sa;
  v=v*ca-u*sa;
  u=ut;

  i1=1:lz-1;
  i2=2:lz;
  rind=fac2*(temp(i1)-temp(i2))-fac1*(u(i1)-u(i2)).^2-fac1*(v(i1)-v(i2)).^2;
  [rind0,indx]=min(rind);
  nmix=0;
  while rind0<0
   if dc
    plot(temp,z,u*10+20,z,v*10+20,z,[15,25],[zmid(indx),zmid(indx)]);
    title(sprintf('Ri mixing -- t = %f',t));
    axis(as);
    drawnow;
  frm=mkpng(frm);
    if nmix<4
%     pause;
    end
   end
%    ['mix ',num2str(indx)]
%    [temp(indx),temp(indx+1)];
    delta=1/8-3/2*rind0*dz*dz/((u(indx)-u(indx+1)).^2+(v(indx)-v(indx+1)).^2);
    mixmat=[1-delta,delta;delta,1-delta];
    temp(indx:indx+1)=mixmat*temp(indx:indx+1);
    u(indx:indx+1)=mixmat*u(indx:indx+1);
    v(indx:indx+1)=mixmat*v(indx:indx+1);
%    i1=indx-1+(indx==1);
%    i2=indx+1-(indx==lz);
    if indx==1
     i1=1:2;i2=2:3;
    elseif indx==lz-1
     i1=lz-2:lz-1;i2=lz-1:lz;
    else
     i1=indx-1:indx+1;i2=indx:indx+2;
    end
    rind(i1)=fac2*(temp(i1)-temp(i2))-fac1*(u(i1)-u(i2)).^2-fac1*(v(i1)-v(i2)).^2;
    [rind0,indx]=min(rind);
    nmix=nmix+1;
    end
% if rem(t,1)==0
   plot(temp,z,u*10+20,z,v*10+20,z);
   legend("T","u","v");
   axis(as);
   title(sprintf('Adjusted Profile -- t = %f',t));
   drawnow;
  frm=mkpng(frm);
%   ts=[ts,temp];
%   return
%   end
 end
