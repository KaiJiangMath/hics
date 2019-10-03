function save_simplice( DIM, m )

format long

tic;
simplex = regular_simplex(DIM);

rotation_simplice = generate_rotation_simplex( DIM, m, simplex );

fname = sprintf('./simplice/%dD_%d_simplice.mat', DIM, m);
save(fname, 'rotation_simplice');

gen_simplice = toc


end

