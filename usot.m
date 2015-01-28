function usot(prec, min_area)


fem = FEM([0 0 1 0 1 1 0 1]', prec, min_area);

N = size(fem.Promoted.nodes, 2);
numofnodes = fem.Num_nodes;


boundary = Boundary();


boundary.set_boundary('x - 1');
boundary.set_boundary('y - 1');
boundary.set_boundary('x');
boundary.set_boundary('y');


[bc1.bc, bc2.bc, bc3.bc, bc4.bc] = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 4);

boundary.setDirichlet(bc1.bc);
boundary.setDirichlet(bc2.bc);
boundary.setDirichlet(bc3.bc);
boundary.setDirichlet(bc4.bc);


[dofs, ndofs] = boundary.dofs(N);


M = fem.assema(1);
S = fem.assems(1);

k = 50.0;
% sigma = 0.2;


Q =  S - (k^2) * M;

X = Q(dofs, dofs);
Y = Q(dofs, ndofs);

G = X\Y;


local_radius = 0.04;

loc_x = 0.5;
loc_y = 0.5;

location = [];

for i = 1:size(dofs, 1)
    if (fem.Promoted.nodes(1,dofs(i)) - loc_x)^2 + (fem.Promoted.nodes(2,dofs(i)) - loc_y)^2 < local_radius^2
        location(size(location, 2) + 1) = i;
    end
end

Sol=zeros(N, 1);

L = G(location,:);

[v, d] = eig((eye(size(ndofs, 1)) +  G'*G)\(L'*L));

disp(d(1, 1))

Sol(ndofs)=v(:,1);

F= -Q*Sol;

Sol(dofs) = Q(dofs, dofs)\F(dofs);



Estimate = zeros(N, 1);
for i = 1:N
    Estimate(i) = besselj(0.,k*sqrt((fem.Promoted.nodes(1,i) - loc_x)^2 + (fem.Promoted.nodes(2,i) - loc_y)^2));
end

tar = Sol(dofs(location));

disp((tar'*tar)/(Sol'*Sol))



min_scale = min(Sol);
max_scale = max(Sol);

if max_scale + min_scale < 0
    scale = min_scale;
else
    scale = max_scale;
end


figure(1);
trisurf(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), Sol(1:numofnodes)/scale);


figure(2);
trisurf(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), Estimate(1:numofnodes));

end

