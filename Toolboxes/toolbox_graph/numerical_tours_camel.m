getd = @(p)path(p,path); % scilab users must *not* execute this

getd('toolbox_signal/');
getd('toolbox_general/');
getd('toolbox_graph/');

name = 'camel';
options.name = name;
[V,F] = read_mesh(name);
N = size(V,2);

clf;
plot_mesh(V,F, options);

i = 1;
U = perform_fast_marching_mesh(V, F, i);

options.method = 'continuous';
J = randperm(N); J = J(1:50);
paths = compute_geodesic_mesh(U, V, F, J, options);

clf;
plot_fast_marching_mesh(V, F, U, paths, options);