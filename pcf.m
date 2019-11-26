
clear
% Particle setup
domain_size = 22.0;
num_particles = 500;

% Calculation setup
dr = 0.1;

% Random arrangement of particles
rMax = domain_size / 4;

%n=readmatrix('problem4_data_set1.dat');
%x=n(:,1);
%y=n(:,2);
%z=n(:,3);

n=readmatrix('problem4_data_set2.dat');
x=n(:,2);
y=n(:,3);
z=n(:,4);

%ee=100;
%rmin=2;
%n = randi([rmin domain_size*ee],num_particles,3)/ee;
%x = n(:,1);
%y = n(:,2);
%z = n(:,3);

f=pairCorrelationFunction_3D(x, y, z, domain_size, rMax, dr);

function f = pairCorrelationFunction_3D(x, y, z, domain_size, rMax, dr)
%{
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        domain_size     length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell
        %}

        % Find particles which are close enough to the box center that a circle of radius
        % rMax will not cross any edge of the box
        bools1 = gt(x , rMax);
        bools2 = lt(x , (domain_size - rMax));
        bools3 = gt(y , rMax);
        bools4 = lt(y , (domain_size - rMax));
        bools5 = gt(z , rMax);
        bools6 = lt(z , (domain_size - rMax));

        interior_indices = find(bools1 .* bools2 .* bools3 .* bools4 .* bools5 .* bools6);
        num_interior_particles = length(interior_indices);

        edges = 0:dr:rMax + 1.1 * dr;
        num_increments = length(edges) - 1;
        g = zeros(num_interior_particles, num_increments);
        radii = zeros(num_increments,1);
        numberDensity = length(x) / (domain_size)^3;


        %
        for p = 1:num_interior_particles
            index = interior_indices(p);
            d = sqrt((x(index) - x).^2 + (y(index) - y).^2 + (z(index) - z).^2);
            d(index) = 2 * rMax;
            [result,]= histcounts(d, (edges),'Normalization','probability'); %'count', 'countdensity', 'cumcount', 'probability', 'pdf', 'cdf'
            %result(:,1)=result(1,:);
            g(p,:) = result / numberDensity;
        end
        %
        g_average = zeros(num_increments,1);

        for i = 1: num_increments
            radii(i) = (edges(i) + edges(i+1)) / 2.;
            rOuter = edges(i + 1);
            rInner = edges(i);
            g_average(i) = mean(g(:, i)) / (4.0 / 3.0 * pi * (rOuter^3 - rInner^3));
        end
        f=g_average;
        plot(f);
        xlabel('distance');
        ylabel('P(r)');
end