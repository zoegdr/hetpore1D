function params = create_parameters

clear params

params.N = input("Space discretisation = ");
params.p = input("Number of periods to run = ");
params.Phi0 = input("Initial porosity = ");
params.nu = input("Poisson ratio = ");
params.Astar = input("Magnitude of applied load = ");
params.omega = input("Frequency of applied load = ");
params.perm_law = input("Permeability law = ");
params.stress_law = input("Stress law = ");
params.damage = input("Type of damage (local, dec_neg, dec_pos, inc_neg, inc_pos) = ");
params.lp = input("Location of permeability damage = ");
params.dp = input("Magnitude of permeability damage = ");
params.ls = input("Location of stiffness damage = ");
params.ds = input("Magnitude of stiffness damage = ");
params.v = input("Variance of Gaussian = ");

end