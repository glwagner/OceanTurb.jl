global qisum qiwin qlsum qlwin qoffsum qoffwin

 rhocpw = 1.024e3*3990;
    rs1 = 0.62;         % red
    rs2 = 1-rs1;        % blue
  qisum = 1400;         % noon amplitude of insolation (W/m**2)
  qiwin = 700;          % noon amplitude winter
  qlsum = -144;         % long wave heat loss (W/m**2)
  qlwin = -117;         % long wave loss winter
qoffsum = 0.0;          % below horizon
qoffwin = -0.4;         % below horizon (winter)

  beta1 = 0.6;          % longwave extinction coefficient (m)
  beta2 = 20;           % shortwave extinction coefficient (m)

    qdm = 0.28885;      % mean of daily cycle
    qdm = 0;            % use daily cycle

zt = z+0.5*dz; % top of cell
zb = z-0.5*dz; % bottom of cell

qb = rs1*exp(-depth/beta1) + rs2*exp(-depth/beta2);

%     r₁ ( e^{zᵀ/β₁) - e^{zᴮ/β₁} )        + r₂ ( e^{zᵀ/β₂} - e^{zᴮ/β₂} ) /          ( (1-qᵇ)*ρ₀*Cᵖ*dz )  
ab = (rs1*(exp(zt/beta1) - exp(zb/beta1)) + rs2*(exp(zt/beta2) - exp(zb/beta2))) / ((1-qb)*rhocpw*dz);
lg = exp(z/beta2);

function [qi, ql, lgt, qd]=heating(t, qd)
  global qisum qiwin qlsum qlwin qoffsum qoffwin
  ct = cos(pi/180 * rem(t, 365.0));
  ql = (qlsum+qlwin)/2 - (qlsum-qlwin)/2*ct;
  qa = (qisum+qiwin)/2 - (qisum-qiwin)/2*ct;
  qo = (qoffsum+qoffwin)/2 - (qoffsum-qoffwin)/2*ct;
  if qd == 0
    qd = (1-qo).*cos(2*pi*(rem(t,1)-0.5)) + qo;
    qd = qd.*(qd>0);
  end
  qi = qa.*qd;
%  qi=qa*0.3;
  lgt = qi./qisum;
endfunction

t = 0:dt:365-dt/2;
[qi, ql, lgt, qd] = heating(t, qdm);
qlsum = -qlsum*sum(qi) / sum(ql);
qlwin = -qlwin*sum(qi) / sum(ql);
mean(qd)
