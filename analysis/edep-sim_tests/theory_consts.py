# ln term consts
m_e = 0.511 # electron mass, MeV
M = 0.511 # incident particle mass, MeV
I = 188 * 1e-6 # MeV
 
z = -1
K = 0.307075 # MeV mol^-1 cm^2
Z = 18 # atomic number
A = 39.948 # atomic mass, g/mol
rho = 1.3973 # g/cm^3, liquid
#rho = 1.663e-3 # gas

steps = 2000

# delta correction

# gasous argon
#x0 = 1.7635
#x1 = 4.4855
# from the sternheimer density effect paper
#x0 = 1.96
#x1 = 4.0
#a = 3.89/10
#k = 2.80
#C = -11.92

# liquid argon. Check these!
x0 = 0.2
x1 = 3.0
C = -5.2146
a = 0.19559
k = 3.0 

# L1 correction
b = 1.8
F = 20 # approximate for 1 to 565 keV

Fx = [3/10 * 1/10, 2/10, 3/10,4/10,5/10,6/10,7/10,9/10,2/10 * 10, 4/10 * 10, 7/10 * 10, 8/10 * 10, 9/10 *10, 14]
Fy = [2/10 * 100, 9/10 * 10, 7/10 * 10, 6/10 * 10, 5/10 * 10, 4/10 * 10, 3/10 * 10, 2/10 * 10, 4/10, 5/10 * 1/10, 8/10 * 1/100, 5/10 * 1/100, 3/10 * 1/100, 1/10**3]