from django.shortcuts import render,HttpResponse
from django.shortcuts import render

# Create your views here.

def home(request):
    # return HttpResponse("This is my homepage(/)")
    # context = {'name':'Hasan', 'project':'SWID'}
    data = {}
    if request.method=="POST":
        # print("This is post")
        PR = request.POST ['PR']
        PT = request.POST['PT']
        AT = request.POST['AT']
        SCIF = request.POST ['SCIF']
        SCAF = request.POST['SCAF']
        GWV = request.POST['GWV']
        AD = request.POST ['AD']
        SR = request.POST['SR']
        
        #data = [PR,PT]

        #print(PR,PT,AT,SCIF, SCAF,GWV,AD,SR)

        #labels = ['pumping_rate', 'pumping_time']
        #values = [PR, PT]

        #data = {'labels': labels, 'data': values}
        import numpy as np
        # Declare variables
        i, j, n, nend, kount1, kount2 = 0, 0, 0, 0, 0, 0
        col = 301
        injw = 301
        obsw = 301 + 36
        pi = 3.14159
        # Stress period 1 (S1)
        ColdS1 = np.zeros(int(col))
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
        
        Q = float(PR)
        ti = float(PT)
        b= float(AT)   
        Ci = float(SCIF)
        Ca = float(SCAF)

        # Get basic input parameters (stress period 2)
        # vx = float(input("Groundwater velocity (ft/day)="))
        # alpha = float(input("Aquifer dispersivity (ft)="))
        # Ret = float(input("Solute retention (-)="))
        vx = float(GWV)
        alpha = float(AD)
        Ret = float(SR)

        print(Q,ti,b)
        # k = float(input("Solute reaction rate (1/days)="))
        # Define
        Q = 192.5* Q  # gpm to ft^3/day
        ti = 0.04166667* ti   # hours to days
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
        WARNING = []

        # WARNINGS :(
        if Pe >= 4.0:
            print("\n WARNING: Peclet # >= 4, numerical dispersion/artificial oscillation\n")
            txt='WARNING: Peclet # >= 4, numerical dispersion/artificial oscillation'
            WARNING.append(txt)
            

        if rmin >= 10.0:
            print("\n WARNING: Minimum radius of influence >= 10 feet, boundary condition effect\n")
            txt1='WARNING: Minimum radius of influence >= 10 feet, boundary condition effect'
            WARNING.append(txt1)

        if rmin < 0.5:
            print("\n WARNING: Minimum radius of influence < 0.5 feet, boundary condition effect\n")
            txt2='WARNING: Minimum radius of influence < 0.5 feet, boundary condition effect'
            WARNING.append(txt2)

        if deltS1 > AdvStab1:
            print("\n WARNING: Advection term unstable for injection phase\n")
            txt3='WARNING: Advection term unstable for injection phase'
            WARNING.append(txt3)

        if deltS2 > AdvStab2:
            print("\n WARNING: Advection term unstable for drift phase\n")
            txt4='WARNING:  Advection term unstable for drift phase'
            WARNING.append(txt4)

        if deltS1 > DispStab1:
            print("\n WARNING: Dispersion term unstable for injection phase\n")
            txt5='WARNING: Dispersion term unstable for injection phase'
            WARNING.append(txt5)

        if deltS2 > DispStab2:
            print("\n WARNING: Dispersion term unstable for drift phase\n")
            txt6='WARNING: Dispersion term unstable for drift phase'
            WARNING.append(txt6)
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

        # Horizontal    
        for i in range(col+col-1):
            Cx[i][0] = -25.0 + (1./12.) * i


        # Right (+) side numerical
        for i in range(col):
            Cx[i+col-1][1] = CnewS1[i]
    

        # Left (-) side numerical
        for i in range(2, col):
            Cx[301-i][1] = CnewS1[i]
    

        x,y = map(list,zip(*Cx))


        x_data = x
        y_data = y

        # Analytical Solution 

        pi = 3.1416
        Ax = []
        for i in range(1, col+col):
            A = Cx[i-1][0]**2 - ((Q*ti)/(3.1416*b*theta))
            AA = ((((((Q*ti)/(3.1416*b*theta))**(3./2.))*16.*alpha))/3.)**(1./2.)
            
            if A/AA >= 0:
                AAA = 1. - (1./(((1.+a1*(A/AA)+a2*(A/AA)**2.+a3*(A/AA)**3.+a4*(A/AA)**4.))**4.))
                AAA = 1. - AAA
            else:
                AAA = -1.*(1. - (1./(((1.+a1*(-A/AA)+a2*(-A/AA)**2.+a3*(-A/AA)**3.+a4*(-A/AA)**4.))**4.)))
                AAA = 1. - AAA
            
            Ax.append((Ci/2.)*AAA)
        

        sumA = 0.
        sumAA = 0.

        for i in range(col-1):
            if Cr[i+1][1] <= Ca + Ca / 10.:
                break
            A = (Cr[i+1][0] - Cr[i][0]) * (1./2.) * (Cr[i+1][1] + Cr[i][1])  # Zero moment
            AA = (Cr[i+1][0] - Cr[i][0]) * (1./2.) * (Cr[i+1][1] * Cr[i+1][0] + Cr[i][1] * Cr[i][0])  # 1st moment
            sumA += A
            sumAA += AA

        comS1 = sumAA / sumA  # Center of mass

        Mi = Q * ti * Ci * 28.3168 / 1000.  # Mass in grams
        Mo = sumA * b * pi * rmin * theta * 28.3168 / 1000.  # Mass in grams
        Mb1 = (1. + (Mo - Mi) / Mi) * 100.
        print('Mass Balance=',Mb1)
    


        # Set boundary conditions at inlet (left side) and outlet (right side)
        ColdS2[0] = Ca  # Inlet
        CnewS2[0] = Ca  # Inlet
        ColdS2[col+col-2] = Ca  # Outlet
        CnewS2[col+col-2] = Ca  # Outlet

        # Set initial condition within boundaries
        for i in range(1, col+col-2):
            ColdS2[i] = y[i]

        nend = 30000
        kount1 = 0
        kount2 = 100
        time = 0.
        i = 0
        deltS2 = 0.001

        t_t = []  # Initialize Ct as an empty list
        Ca_Ca =[] 
        for n in range(1, nend+1):
            for j in range(1, col+col-3):
                A = vx * alpha * ((ColdS2[j+1] - 2. * ColdS2[j] + ColdS2[j-1]) / (delx**2.)) - vx * ((ColdS2[j] - ColdS2[j-1]) / (delx)) + k * ColdS2[j]
                CnewS2[j] = A * deltS2 / Ret + ColdS2[j]

            for j in range(1, col+col-3):
                ColdS2[j] = CnewS2[j]

            time += deltS2
            kount1 += 1



            #Output S2
            if kount1 == kount2:
                kount1 = 0

                i += 1
                t_t.append(time)
                #print(time, CnewS2[injw])
                Ca_Ca.append(CnewS2[injw])
        
        sumA = 0.0
        sumAA = 0.0

        for i in range(col-1):
            if  Ca_Ca[i+1]<= Ca + Ca/10.0:
                break
            
            A = (t_t[i+1] - t_t[i]) * (1./2.) * (Ca_Ca[i+1] + Ca_Ca[i])
            AA = (t_t[i+1] - t_t[i]) * (1./2.) * (Ca_Ca[i+1] * t_t[i+1] + Ca_Ca[i]* t_t[i])
            
            sumA += A
            sumAA += AA

        Mi = Mo
        Mo = sumA * vx * b * pi * rmin * theta * 28.3168 / 1000.0
        Mb2 = (1. + (Mo - Mi) / Mi) * 100.0
        print('Mass Balance=',Mb2)

    
        data={'pumping_rate':PR,'pumping_time':PT, 'x_data': x_data,'y_data': y_data, 't_t':t_t, 'Ca_Ca':Ca_Ca, 'Ax':Ax, 'WARNING':WARNING,'Mb1':Mb1, 'Mb2':Mb2 }

    return render(request,'home.html',data)


def about(request):
    # return HttpResponse("This is my Aboutpage(/about)")
    return render(request,'about.html')

def contact(request):
    #return HttpResponse("This is my contactpage(/contact)")
    return render(request,'contact.html')
