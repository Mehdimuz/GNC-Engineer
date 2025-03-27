function drv=orbit(t,rv)
    mu=398600;
    r=sqrt((rv(1)^2)+(rv(2)^2)+(rv(3)^2));
    drv=[rv(4); rv(5); rv(6); -(mu/r^3)*rv(1); -(mu/r^3)*rv(2); -(mu/r^3)*rv(3)];
end