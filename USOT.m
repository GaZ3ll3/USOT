[p,t,e, P,T,E] = MeshGen('prq32.0a0.0001zeB',[0. 0.; 1. 0.; 1. 1.; 0. 1.]');
[I,J,U,V] = MatrixAssem(p,t,P,T);
S = sparse(I,J,V);
M = sparse(I,J,U);

k = 50.0;
j = sqrt(-1);
sigma = 0.2;
Q =  S - (k^2 + j*k*sigma) * M;

% figure(1);
% triplot(T' + 1, P(1,:), P(2,:));

% figure(2);
% spy(M);

% figure(3);
% spy(S);

b = unique(E + 1);


freenodes=setdiff(1:p,b); 

X = Q(freenodes, freenodes);

Y = Q(freenodes, b);

G = X\Y;

local_radius = 0.02;

loc_x = 0.5;
loc_y = 0.5;

location = [];

for i = 1:size(freenodes, 2)
    if (P(1,freenodes(i)) - loc_x)^2 + (P(2,freenodes(i)) - loc_y)^2 < local_radius^2
        location(size(location) + 1) = i;
    end
end


Sol=zeros(p,1);

L = G(location,:);

[v, d] = eig((eye(size(b, 1)) +  G'*G)\(L'*L));

Sol(b)=v(:,1);

F= -Q*Sol;

Sol(freenodes) = Q(freenodes, freenodes)\F(freenodes);


Estimate = zeros(p, 1);
for i = 1:p
    Estimate(i) = besselj(0.,k*sqrt((P(1,i) - loc_x)^2 + (P(2,i) - loc_y)^2));
end


ratio = (Sol'*Sol)/(Estimate'*Estimate);

figure(1)
trisurf(T' + 1,P(1,:),P(2,:),real(Sol.*conj(Sol)),'EdgeColor','none');
shading interp
colorbar;
figure(2);
trisurf(T' + 1,P(1,:),P(2,:),ratio*real(Estimate.*conj(Estimate)),'EdgeColor','none');
shading interp
colorbar;
figure(3);
trisurf(T' + 1,P(1,:),P(2,:),ratio*real(Estimate.*conj(Estimate)) -real(Sol.*conj(Sol)) ,'EdgeColor','none');
shading interp
colorbar;
Energy = sigma*Sol.*conj(Sol);

