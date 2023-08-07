import numpy as np

# Declare variables
i, j, n, nend, kount1, kount2 = 0, 0, 0, 0, 0, 0
col = 301
injw = 301
obsw = 301 + 36
pi = 3.14159

# Stress period 1 (S1)
ColdS1 = np.zeros(col)
CnewS1 = np.zeros(col)
r = np.zeros(col)
vr = np.zeros(col)
Cr = np.zeros((col, 2))
a1, a2, a3, a4 = 0.278393, 0.230389, 0.000972, 0.078108


# Stress period 2 (S2)
ColdS2 = np.zeros(col + col - 1)
CnewS2 = np.zeros(col + col - 1)
x = np.zeros(col + col - 1)
Cx = np.zeros((col + col - 1, 3))


# Breakthrough curve
Ct = np.zeros((col - 1, 3))


Q, b, theta, d, ti, Ci, Ca = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
vx, alpha, Ret, k = 0.0, 0.0, 0.0, 0.0
delr, delx, deltS1, deltS2, time, vol, comS1, matS2, vcalc, Mi, Mo, Mb = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
Pe, rmin, AdvStab1, AdvStab2, DispStab1, DispStab2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

# Get basic input parameters (stress period 1)
# Q = float(input("Pumping rate (gpm)="))
# ti = float(input("Pumping time (hours)="))
# b = float(input("Aquifer thickness (ft)="))
# Ci = float(input("Solute concentration injection fluid (mg/L)="))
# Ca = float(input("Solute concentration aquifer fluid (mg/L)="))

Q = .25
ti = 7
b= 5
Ci = 484
Ca = .5

# Get basic input parameters (stress period 2)
# vx = float(input("Groundwater velocity (ft/day)="))
# alpha = float(input("Aquifer dispersivity (ft)="))
# Ret = float(input("Solute retention (-)="))

vx = 1.05
alpha = .15
Ret = 1


# k = float(input("Solute reaction rate (1/days)="))

# Define
Q = Q * 192.5  # gpm to ft^3/day
ti = ti * 0.04166667  # hours to days
k = 0.0  # Turn off k for now, this needs work...
theta = 0.25  # Porosity fixed, for now...
d = 1.0 / 12.0  # Grid spacing of one inch

delr = d  # Grid space
delx = d  # Grid space
deltS1 = 0.00001  # Time

Pe = ((Q / (2.0 * pi * delr * b * theta)) * delr) / ((Q / (2.0 * pi * delr * b * theta)) * alpha)  # Peclet needs to be less than 4 for stability
rmin = ((Q * ti) / (pi * b * theta)) ** (1.0 / 2.0)  # Radius of influence needs to be less than 10 feet for no boundary effects

AdvStab1 = Ret / ((Q / (2.0 * pi * delr * b * theta)) / delr)  # See MT3DMS Guide Pg. 54
AdvStab2 = Ret / (vx / delx)  # See MT3DMS Guide Pg. 54

DispStab1 = (0.5 * Ret) / ((alpha * (Q / (2.0 * pi * delr * b * theta))) / (delr ** 2.0))  # See MT3DMS Guide Pg. 54
DispStab2 = (0.5 * Ret) / ((alpha * vx) / (delx ** 2.0))  # See MT3DMS Guide Pg. 54


# ChemRxnS1 = 1./?
# ChemRxnS2 = 1./?

# WARNINGS :(
if Pe >= 4.0:
    print("\n WARNING: Peclet # >= 4, numerical dispersion/artificial oscillation\n")

if rmin >= 10.0:
    print("\n WARNING: Minimum radius of influence >= 10 feet, boundary condition effect\n")

if rmin < 0.5:
    print("\n WARNING: Minimum radius of influence < 0.5 feet, boundary condition effect\n")

if deltS1 > AdvStab1:
    print("\n WARNING: Advection term unstable for injection phase\n")

if deltS2 > AdvStab2:
    print("\n WARNING: Advection term unstable for drift phase\n")

if deltS1 > DispStab1:
    print("\n WARNING: Dispersion term unstable for injection phase\n")

if deltS2 > DispStab2:
    print("\n WARNING: Dispersion term unstable for drift phase\n")

# Fill arrays with zero placeholders
ColdS1 = [0.0] * col
CnewS1 = [0.0] * col
r = [0.0] * col

Cr = [[0.0] * 2 for _ in range(col)]

ColdS2 = [0.0] * (col + col - 1)
CnewS2 = [0.0] * (col + col - 1)
x = [0.0] * (col + col - 1)

Cx = [[0.0] * 2 for _ in range(col + col - 1)]

Ct = [[0.0] * 2 for _ in range(col - 1)]

# Set boundary conditions at inlet (left side) and outlet (right side)
ColdS1[0] = Ci  # Inlet
CnewS1[0] = Ci  # Inlet
ColdS1[col-1] = Ca  # Outlet
CnewS1[col-1] = Ca  # Outlet



# Set initial condition within boundaries
for j in range(1, col-1):
    ColdS1[j] = Ca

# Set radius
r[0] = 1.0 / 12.0  # Must NOT divide by zero, see below for vr
for j in range(1, col):
    r[j] = delr + r[j-1]

# Set velocity as a function of radius, i.e., v(r)
for j in range(col):
    vr[j] = Q / (2.0 * pi * r[j] * b * theta)


# Compute Stress Condition 1

nend = int(ti / deltS1)

time = 0.0
vol = 0.0
rmin = 0.0

for n in range(1, nend+1):
    for j in range(1, col-1):
        A = alpha * ((vr[j+1] + vr[j]) / 2.0 * (ColdS1[j+1] - ColdS1[j]) / delr - (vr[j] + vr[j-1]) / 2.0 * (ColdS1[j] - ColdS1[j-1]) / delr) / delr - \
            (vr[j] + vr[j-1]) / 2.0 * (ColdS1[j] - ColdS1[j-1]) / delr + k * ColdS1[j]
        CnewS1[j] = A * deltS1 / Ret + ColdS1[j]

    for j in range(1, col-1):
        ColdS1[j] = CnewS1[j]

    time += deltS1
    vol = Q * time
    rmin = (vol / (pi * b * theta))**(1./2.)

# Radial
for i in range(col):
    Cr[i][0] = r[i]
    Cr[i][1] = CnewS1[i]
    

# # Horizontal
# Cx = []
# for i in range(1, col+col):
#     Cx.append([-25.0 + (1./12.) * (i-1), 0])

# # Right (+) side numerical
# for i in range(1, col+1):
#     Cx[i+col-2][1] = CnewS1[i-1]

# # Left (-) side numerical
# for i in range(2, col+1):
#     Cx[302-i][1] = CnewS1[i-1]

print(Cx)


# Horizontal
for i in range(col+col-1):
    Cx[i][0] = -25.0 + (1./12.) * i


# Right (+) side numerical
for i in range(col):
    Cx[i+col-1][1] = CnewS1[i]
    

# Left (-) side numerical
for i in range(2, col):
    Cx[301-i][1] = CnewS1[i]
    
print("New Print")
print(Cx)

x,y = map(list,zip(*Cx))