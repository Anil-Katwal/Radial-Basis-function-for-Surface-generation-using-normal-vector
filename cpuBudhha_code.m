% Data processing column by column
load('buddha_2.mat');
norms_data = normals(1:20:end,:);
dsites = dsites(1:20:end,:);
boundary_dsites = dsites(1:100:end,:);  
boundary_norms = norms_data(1:100:end,:);
outside_data = boundary_dsites + boundary_norms * 0.0005;  

x = outside_data(:,1);
y = outside_data(:,2);
z = outside_data(:,3);

xc = dsites(:,1);
yc = dsites(:,2);
zc = dsites(:,3);

ctre = [xc, yc, zc];

N = length(ctre);
cntrs = [x y z]; M = length(cntrs);
ctrs = [ctre; cntrs];
rhs = [zeros(N, 1); -ones(M, 1)];
c = 700;
%Distance calculation........
DM_data = pdist2(ctrs, ctrs);
IM = sqrt((c * DM_data).^2 + 1);
coe = IM \ rhs;

bmin = min(ctrs, [], 1);
bmax = max(ctrs, [], 1);

% Define the mesh grid
neval = 70;
xgrid = linspace(bmin(1), bmax(1), neval);
ygrid = linspace(bmin(2), bmax(2), neval);
zgrid = linspace(bmin(3), bmax(3), neval);
[xe1, ye1, ze1] = meshgrid(xgrid, ygrid, zgrid);
epoints = [xe1(:) ye1(:) ze1(:)];
%partition for jobs..................
S = 500;
part = ceil(size(epoints, 1) / S);
Pf = zeros(size(epoints, 1), 1);
for i = 1:part
    Start = (i - 1) * S + 1;
    End = min(i * S, size(epoints, 1));
    epointsp = epoints(Start:End, :);
    DM_eval = pdist2(epointsp, ctrs);
   rbf = @(c, r) sqrt((c * r).^2 + 1);
    EM = rbf(c, DM_eval);
    Pf(Start:End) = EM * coe;
end

% Transform Pf for 3D visualization.
Pf = reshape(Pf, neval, neval, neval);

% Plot 3D interpolation
figure(2);
p = patch(isosurface(xe1, ye1, ze1, Pf, 0));
isonormals(xe1, ye1, ze1, Pf, p);
set(p, 'FaceColor', 'cyan', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);
axis off;
camlight;
lighting gouraud;
title('3D Buddha surface reconstruction');

% Plot the data points
figure;
plot3(x, y, z, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
hold on;
plot3(xc, yc, zc, 'b.');
