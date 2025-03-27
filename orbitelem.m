function [e,a,inclination,RAAN,AOP,T_anomaly,Tp]=orbitelem(r,v)
mu=398600*10^9;
i=[1 0 0];
j=[0 1 0];
k=[0 0 1];
h=cross(r,v);
e=(1/mu)*(((dot(v,v)-(mu/sqrt(dot(r,r))))*r)-dot(r,v)*v);
n=cross(k,h);
RAAN=acosd(dot(i,n)/sqrt(dot(n,n)));
inclination=acosd(dot(k,h)/sqrt(dot(h,h)));
AOP=acosd(dot(e,n)/(sqrt(dot(n,n))*sqrt(dot(e,e))));
T_anomaly=acosd(dot(e,r)/(sqrt(dot(r,r))*sqrt(dot(e,e))));
if n(1,2)>0
    RAAN=360-RAAN;
end
if e(1,3)<0
    AOP=360-AOP;
end
if dot(r,v)<0
    T_anomaly=360-T_anomaly;
end
evector=e;
e=sqrt(dot(e,e));
hvector=h;
h=sqrt(dot(h,h));
nvector=n;
n=sqrt(dot(n,n));
a=(h^2)/(mu*(1-e^2));
Tp=2*pi*sqrt((a^3)/mu);
end