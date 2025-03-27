function psi = Psi(q)
psi = [q(4)*eye(3)-cpm(q(1:3)); -q(1) -q(2) -q(3)]; 