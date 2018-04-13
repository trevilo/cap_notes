import numpy as np

def formD2(h, Nx):
    # construct spatial operator (just 2nd derivative)
    D2 = np.zeros((Nx,Nx))

    # interior
    D2 = (     np.diag(np.ones(Nx-1), k=1) +
          -2.0*np.diag(np.ones(Nx  ), k=0) +
               np.diag(np.ones(Nx-1), k=-1))

    # boundary (for Neumann)
    D2[0, 0] = -2.0
    D2[0, 1] =  2.0

    D2[-1, -1] = -2.0
    D2[-1, -2] =  2.0

    # divide by dx^2
    D2 = D2/(h*h)
    return D2



def solve(t0, tf, xL, xR, Nt, Nx, alfa, beta, eta0, plot=False):

    # space and time grids
    xx,dx = np.linspace(xL, xR, num=Nx, endpoint=True, retstep=True)
    tt,dt = np.linspace(t0, tf, num=Nt, endpoint=True, retstep=True)

    # full space-time solution
    U = np.zeros((Nx,Nt))

    # set IC
    U[:,0] = eta0

    # form operators
    D2 = formD2(dx, Nx)
    Amat = np.identity(Nx) - 0.5*dt*D2
    Bmat = np.identity(Nx) + 0.5*dt*D2

    cvec = np.zeros(Nx)

    # march in time
    for jj in range(1,Nt):
        cvec[0]  = -0.5*dt*( 2.0*alfa[jj-1]/dx+2.0*alfa[jj]/dx);
        cvec[-1] = -0.5*dt*(-2.0*beta[jj-1]/dx-2.0*beta[jj]/dx);

        Bu = np.matmul(Bmat, U[:,jj-1])
        
        rvec = Bu + cvec

        U[:,jj] = np.linalg.solve(Amat, rvec)


    # plot if requested
    if plot:
        import matplotlib
        matplotlib.Use('Agg')
        import matplotlib.pyplot as plt

        plt.plot(xx, U[:,0], 'b--', linewidth=2, label="Initial")
        plt.plot(xx, U[:,-1], '--', linewidth=2, label="Final")
        plt.savefig('soln.pdf', bbox_inches='tight')

    # return full solution
    return xx, tt, U


#################################################################
# If standalone invocation, run a case
#################################################################
if __name__ == '__main__':
    import sys
    from scipy.signal import square


    # Set up step current case
    
    # Parameters
    kappa = 0.0195174;   #--- sec/m
    sigma = 52.1;        #--- sec/m
    gamma = kappa/sigma;
    L = 50e-6;           #--- m
    C = 0.03134;         #--- F/m2
    a = 4.19956e7/C;     #--- m
    V0 = 1.25;           #--- volt
    Iunscaled = 200;     #--- Amp/m^2
    Ls = 25e-6;          #--- m
    kappa_s = 0.0311627; #--- sec/m
    Time = 4;            #--- sec
    
    Nt = 1024;           #--- number of time steps
    Nx = 256;            #--- number of space steps


    ## Domain
    t = np.linspace(0,Time , Nt);  #--- sec
    x = np.linspace(0,L , Nx);

    # convert to nondimesional values
    tau = (kappa*sigma/(kappa+sigma))*(1.0/(a*C*L*L)) * t;
    xi = x / L;
    Nxi = Nx
    Ntau = Nt

    ## Applied Current
    Iamp = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa)));

    # square wave case
    #I = Iamp*square(6.0*np.pi*tau);

    # sine case
    I = Iamp*np.sin(5.0*np.pi*tau);


    ##========================================================
    ## Solution of HF model
    #
    # Neumann BCs:
    #-- alpha at xi=0
    alpha = -I*(gamma/(1+gamma));
    #-- beta at xi=1
    beta  =  I*(1/(1+gamma));

    # IC:
    eta0 = np.zeros(Nx)

    # solve
    xg, tg, U = solve(tau[0], tau[-1], xi[0], xi[-1], Ntau, Nxi, alpha, beta, eta0, plot=False)

    # plot
    import matplotlib.pyplot as plt
    plt.contourf(xg, tg, U.transpose(), 64, cmap='inferno')
    plt.colorbar()
    plt.xlabel(r'$\xi$', fontsize=22)
    plt.ylabel(r'$\tau$', fontsize=22)

    plt.contour(xg, tg, U.transpose(), 64, cmap='inferno')
    plt.savefig('full_soln.pdf', bbox_inches='tight')
    
