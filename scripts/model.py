import sys
import numpy
from scipy import integrate
from scipy import special


n_aa_per_phage = 1225786
n_proteins_per_cell = 3.5e6
mean_mass_aa = 118.9 # daltons (g/mol)
#mean_n_aa_per_protein_cell = 343.2
#mass_per_protein_cell = mean_n_aa_per_protein_cell * mass_aa
#mean_mass_
mean_mass_per_protein_cell = 40000 # daltons (g/mol)
#lysis_rate = 1.14
lysis_rate_per_phage_per_time = 1e-7

phi_r_max = 0.55
#phi_r_0 = 0.2
#phi_p = 0.1

# converstion of mass fraction of phage protein in cell to number of phage
#mass_protein_cell = 156*1e-15 # (g)
mean_length_protein = 343.2

# death rate per-unit time
delta_phage = 1e-5


v = 3 # nutritional efficiency
# convert OD600 to cells 
od600_to_cell_density = 609375000
n_cell_0 = od600_to_cell_density*0.01
n_phage_0 = 1e5

n_cell = 1e8


#phage_conversion_factor = mass_protein_cell * (mean_mass_per_protein_cell**-1) * mean_n_aa_per_protein_cell * (n_aa_per_phage**-1)* lysis_rate


def calculate_lambda_cell_and_phage(phi_r_0, phi_p):

    # units of inverse time
    lambda_cell = (phi_r_max - phi_p - phi_r_0)*v * (od600_to_cell_density**-1)

    # phage per unit cell
    Y = n_proteins_per_cell*mean_mass_per_protein_cell*(mean_mass_aa**-1)*(n_aa_per_phage**-1)

    # units of inverse time and inverse cells
    lambda_phage = phi_p*Y*lysis_rate_per_phage_per_time
    #phage_birth_per_cell_per_time = n_proteins_per_cell*mean_mass_per_protein_cell*(mean_mass_aa**-1)*(n_aa_per_phage**-1)*lysis_rate_per_phage_per_time

    return lambda_cell, lambda_phage, Y


    




def derivative(X, t, lambda_cell, phage_birth_per_cell_per_time, lysis_rate_per_phage_per_time, death_phage):
    
    cell, phage = X

    deriv_cell = cell * od600_to_cell_density * (lambda_cell - lysis_rate_per_phage_per_time * phage)
    deriv_phage = phage * (phage_birth_per_cell_per_time*cell  - death_phage)
    
    return numpy.array([deriv_cell, deriv_phage])


    Nt = 1000
    tmax = 24
    t = numpy.linspace(0.,tmax, Nt)

    X0 = [n_cell_0, n_phage_0]
    #res = integrate.odeint(derivative, X0, t, args = (lambda_cell, phage_birth_per_cell_per_time, lysis_rate_per_phage_per_time, death_phage))
    #cell_t, phage_t = res.T

    #print(cell_t)


    #print(cell_t)


def lotka_volterra_V(n_cell, n_phage, phi_r_0, phi_p):

    lambda_cell, lambda_phage, Y = calculate_lambda_cell_and_phage(phi_r_0, phi_p)

    V = lambda_phage*n_cell - delta_phage*numpy.log(n_cell) + lambda_phage*(Y**-1) * n_phage - lambda_cell*numpy.log(n_phage)

    return V





def solve_n_cell_max(n_cell_0, n_phage_0, phi_r_0, phi_p):

    lambda_cell, lambda_phage, Y = calculate_lambda_cell_and_phage(phi_r_0, phi_p)

    V = lotka_volterra_V(n_cell_0, n_phage_0, phi_r_0, phi_p)

    #print(lambda_cell, lambda_phage, Y )

    #lambda_phage*(delta_phage**-1)

    lambert_w_term = delta_phage*(lambda_phage**-1) * numpy.exp( (V - lambda_cell*( numpy.log(lambda_cell*Y*(lambda_phage**-1)) - 1)) * (lambda_phage**-1))

    #(lambda_cell/death_phage) *
    #print(lambert_w_term)

    W = special.lambertw(lambert_w_term, k=0)

    W_real = numpy.real(W)

    n_cell_max = W_real*lambda_phage * (delta_phage**-1)

    print(n_cell_max)



# requirements for cell growth
# phi_r_max < phi_p + phi_r_0  


phi_r_0 = 0.2
phi_p = 0.1


if (phi_r_max <= phi_r_0 + phi_p):
    sys.exit('Cell growth not possible!')


solve_n_cell_max(1e2, 1e5, 0.2, 0.3)


# 


