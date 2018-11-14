tmax=3*365;
t0=100;
dz=2;
depth=80;%200;
lat=31;
dt=1/32;
dts=dt*86400;

z=-[0.5*dz:dz:depth]';
tau0=0.1;

heatinginit;

lz=length(z);
zmid=[z(1:lz-1)-dz/2]';

%T=4+11*exp(z/500)+10*exp(z/10);
%T=18+5*exp(z/30);
T=0*z+18.5;
u=0*z;
v=0*z;

as=[14 26 -depth 0];
rho0=1000;
cp=3900;
alphag=2e-4*9.81;

dt=1/32;

f=2*2*pi*sin(lat/180*pi);
sa=sin(f*dt/2);
ca=cos(f*dt/2);

%us=u(1);
%vs=v(1);
fac1=0.25/dz/dz;
fac2=alphag/dz;

function [mixmat,jf]=cadj(b)
  l=length(b);
  jf=1;
  while jf<l
    if mean(b(1:jf))>b(jf+1)
      break;
    end
    jf=jf+1;
  end
  if jf==l
    mixmat=ones(l,l)/l;
  else
    mixmat=[ones(jf,jf)/jf,zeros(jf,l-jf);zeros(l-jf,jf),eye(l-jf)];
  end
endfunction

function [mixmat,jf]=bulkri(b,u,v,j0,dz)
  Ri_b=0.65;
  l=length(b);
  jf=j0;
% bulk ri
  while jf<l
    mld=jf*dz;
    rind=Ri_b*(mean(u(1:jf))-u(jf+1))^2+Ri_b*(mean(v(1:jf))-v(jf+1))^2-(b(jf+1)-mean(b(1:jf)))*mld;
    if rind>0
      break;
    end
    jf=jf+1;
  end
  if jf==l
    mixmat=ones(l,l)/l;
  else
    mixmat=[ones(jf,jf)/jf,zeros(jf,l-jf);zeros(l-jf,jf),eye(l-jf)];
  end
endfunction


plot(T,z,u*10+20,z,v*5+10,z)
title(['t = ',num2str(t0)])
axis(as);
drawnow;

day_sv=t0;
T_sv=T;
u_sv=u;
v_sv=v;
mld_sv=0;

ca_cnt=[0];
bri_cnt=[0];
ri_cnt=[0];

for t=t0+dt:dt:t0+tmax

  [qi,ql,lgt]=heating(t,qdm);
  [qip,qlp,lgtp]=heating(t+dt,qdm);

  % Forward Euler timestep for temperature
  T=T+qi*ab*dts; % flux divergence of shortwave solar insolation
  T(1)=T(1)+ql/rhocpw/dz*dts; % flux divergence of longwave solar insolation 

  [mixmat,ind]=cadj(T); % convective adjustment. mixmat: homogenizer; ind: mixed layer depth index
  T=mixmat*T;
  u=mixmat*u;
  v=mixmat*v;

  mld=ind*dz; % mixed layer depth

  ca_cnt=[ca_cnt,ind];


  tau=tau0*24*3600/rho0/mld*dt; % constant wind stress

  ut=u*ca+v*sa; % u at next timestep (analytical integration of u_t + im f u = 0)

  % Update v
  v=v*ca-u*sa;
  u=ut;

  % distribute windstress over mixed layer (may need to divide by mld)
  du=tau*(z > -mld);
  u=u+du;

  ut=u*ca+v*sa;
  v=v*ca-u*sa;

  u=ut;

  [mixmat,ind]=bulkri(alphag*T,u,v,ind,dz);
  T=mixmat*T;
  u=mixmat*u;
  v=mixmat*v;
%  b=mixmat*b;

  bri_cnt=[bri_cnt,ind];

  i1=1:lz-1;
  i2=2:lz;
  rind=fac2*(T(i1)-T(i2))-fac1*(u(i1)-u(i2)).^2-fac1*(v(i1)-v(i2)).^2;
  [rind0,indx]=min(rind);
  nmix=0;
  while rind0<0
%    ['mix ',num2str(indx)]
%    [T(indx),T(indx+1)];
    delta=1/8-3/2*rind0*dz*dz/((u(indx)-u(indx+1)).^2+(v(indx)-v(indx+1)).^2);
    mixmat=[1-delta,delta;delta,1-delta];
    T(indx:indx+1)=mixmat*T(indx:indx+1);
    u(indx:indx+1)=mixmat*u(indx:indx+1);
    v(indx:indx+1)=mixmat*v(indx:indx+1);
    if indx==1
      i1=1:2;i2=2:3;
    elseif indx==lz-1
      i1=lz-2:lz-1;i2=lz-1:lz;
    else
      i1=indx-1:indx+1;i2=indx:indx+2;
    end
    rind(i1)=fac2*(T(i1)-T(i2))-fac1*(u(i1)-u(i2)).^2-fac1*(v(i1)-v(i2)).^2;
    [rind0,indx]=min(rind);
    nmix=nmix+1;
  end
  pt_cnt=[ri_cnt,nmix];

  if rem(t,1)==0
%    plot(T,z,u*10+20,z,v*10+20,z);
%    axis(as);
%    title(['t = ',num2str(t)])
%    drawnow;;
    t
    fflush(1);
    day_sv=[day_sv,t];
    T_sv=[T_sv,T];
    u_sv=[u_sv,u];
    v_sv=[v_sv,v];
    mld_zv=[mld_sv,mld];
				%   return
  end
%  ts=[ts,T];
%  us=[us,u];
%  vs=[vs,v];
end

imagesc(T_sv);
