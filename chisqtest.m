%%Chi-Square Test by Rizve Chowdhury
%Description: The Chi-Square Test culls disruptive measurements that might interfere with the UKF's operation.
%Input Variables:
%1) r is the measured attitude matrix at time step K+1
%2) M is the propagated mean of the attitude matrix at time step K+1
%3) R is the covariance of the measured attitude matrix at time step K+1
%Output Variable
%t is a Boolean that outputs a 1 if alterations to the UKF need to be made and a 0 if
%alteration to the UKF are not required.
function t= chisqtest(r,M,R)
q=(r-M)*inv(R)*(r-M)';
qref=12.592;
if q>qref
    t=1;
else
    t=0;
end
end