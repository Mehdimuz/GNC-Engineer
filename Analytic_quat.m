function [q_hat] = Analytic_quat(Vecs, bVecs, sig)
%ANALYTIC_QUAT lalalalalala
%   Fun, simple function that returns the optimal quaternion for if you
%   have two references frame vectors and two body frame vectors. If there
%   are more than two vectors, Davenport-q/OLAE/FOAM etc. should be used
%   for an optimal solution.
%
%   ARGS:
%       Vecs:  3xn of vectors as observed in the reference frame
%             (i.e., [v_1, v_2, ..., v_n]
%       bVecs: 3xn of vectors as observed in the body frame
%             (i.e., [b_1, b_2, ..., b_n]
%       sig:   1xn uncertainty of each measurement
%   OUTPUT:
%       q_hat: 4x1 quaternion relating reference vectors to body vectors
%
%   NOTE: if more than 2 vectors are provided, the algorithm picks the two
%   vectors with the lowest uncertainty. However, it can often perform
%   better if two vectors with the greatest orthogonality is provided. I
%   haven't done a covariance analysis so I don't know which is better. Not
%   a terribly big deal to figure it out though.
cpm = @(r) [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
Xi  = @(q) [q(4)*eye(3)+cpm(q(1:3)); -q(1) -q(2) -q(3)];
Psi = @(q) [q(4)*eye(3)-cpm(q(1:3)); -q(1) -q(2) -q(3)];

[~, n] = size(Vecs);
if nargin == 3
    m1 = 1; m2 = 2;
else
    [s, ind] = sort(sig);
    m1 = ind(1); m2 = ind(2);
end

for i = 1:n
    Vecs(:,i) = Vecs(:,i) ./ norm(Vecs(:,i));
    bVecs(:,i) = bVecs(:,i) ./ norm(bVecs(:,i));    
end

v1 =  Vecs(:,m1); v2 =  Vecs(:,m2);
b1 = bVecs(:,m1); b2 = bVecs(:,m2);

q_hat = [cross(b1-v1,b2-v2);  b2'*v1 - b1'*v2];
q_hat = q_hat./norm(q_hat).*sign(q_hat(4));

end

