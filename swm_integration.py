## TIME INTEGRATION
def time_integration(u,v,eta,k,omega):

    tic = tictoc.time()     # measure time
    global dt
    dt = param['dt']        # for convenience
    t = param['t0']         # initial time
    feedback_ini(u,v,eta,k,omega,t)   # output

    ## RUNGE KUTTA 4th ORDER

    # RUNGE KUTTA coefficients
    rk_a = np.array([1/6.,1/3.,1/3.,1/6.])
    rk_b = np.array([0.5,0.5,1.])

    # trigger deep copy through [:]
    u0,v0,eta0,k0,omega0 = u.copy(),v.copy(),eta.copy(),k.copy(),omega.copy()
    u1,v1,eta1,k1,omega1 = u.copy(),v.copy(),eta.copy(),k.copy(),omega.copy()

    global i
    for i in range(param['Nt']):

        # trigger deep copy through [:]
        u1[:],v1[:],eta1[:],k1[:],omega1[:] = u,v,eta,k,omega

        for rki in range(4):
            du,dv,deta,dk,domega = rhs(u1,v1,eta1,k1,omega1)

            if rki < 3: # RHS update for the next RK-step
                u1 = u + rk_b[rki]*dt*du
                v1 = v + rk_b[rki]*dt*dv
                eta1 = eta + rk_b[rki]*dt*deta
                k1 = k + rk_b[rki]*dt*dk
                omega1 = omega + rk_b[rki]*dt*domega

            # Summing all the RHS on the go
            u0 += rk_a[rki]*dt*du
            v0 += rk_a[rki]*dt*dv
            eta0 += rk_a[rki]*dt*deta
            k0 += rk_a[rki]*dt*dk
            omega0 += rk_a[rki]*dt*domega

        # deep copy through [:]
        u[:],v[:],eta[:],k[:],omega[:] = u0,v0,eta0,k0,omega0

        t += dt
        feedback(u,v,eta,k,omega,t,tic)

    print(('Integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))
    output_txt(('\nTime integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))

    # finalising output
    if param['output']:
        output_nc_fin()         # finalise nc file
        output_txt_fin()        # finalise info txt file

    return u,v,eta,k,omega

### TIME STEPPING SCHEMES
def RK4(u,v,eta):
    """ Computes the right-hand side using RUNGE KUTTA 4th order scheme.
    u,v,h are coupled in every of the 4 sub-time steps of the RK4 scheme."""
    k1 = rhs(u,v,eta)
    k2 = rhs(u + dt/2.*k1[0],v + dt/2.*k1[1],eta + dt/2.*k1[2])
    k3 = rhs(u + dt/2.*k2[0],v + dt/2.*k2[1],eta + dt/2.*k2[2])
    k4 = rhs(u + dt*k3[0],v + dt*k3[1],eta + dt*k3[2])

    du = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.
    dv = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.
    deta = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6.

    return du,dv,deta

## FEEDBACK ON INTEGRATION
def feedback_ini(u,v,eta,k,omega,t):
    """ This function is only evaluated right before the integration loops starts.
    Initializes all files ready to allow output."""
    if param['output']:
        output_nc_ini()
        output_nc(u,v,eta,k,omega,t)  # store initial conditions
        output_param()      # store the param dictionnary as .npy
        output_param_txt()  # store the param dictionnary easily readable as .txt


        # Store information in txt file
        output_txt('Integrating %.1f days with dt=%.2f min in %i time steps' % (param['Ndays'],dt/60.,param['Nt']))
        output_txt('Time integration scheme is '+param['scheme']+' with CFL = %.2f' % param['cfl'])
        output_txt('')
        output_txt('Starting shallow water model on '+tictoc.asctime())
        print(('Starting shallow water model run %i on ' % param['run_id'])+tictoc.asctime())
    else:
        print('Starting shallow water model on '+tictoc.asctime())

def feedback(u,v,eta,k,omega,t,tic):
    """ Executed on every time step. Check whether it is time to produce netcdf
    output or give feedback on the progress. """

    if (i+1) % param['output_n'] == 0:  # check whether it is time to produce netcdf output.
        if param['output']:     # storing u,v,nc,k,omega as netCDF4
            output_nc(u,v,eta,k,omega,t)

    # feedback on progress every 5% step.
    if ((i+1)/param['Nt']*100 % 5) < (i/param['Nt']*100 % 5):
        progress = str(int((i+1)/param['Nt']*100.))+'%'
        print(progress, end='\r')
        if i > 100:
            output_txt(progress,'\n')

    # estimate total time for integration after 100 time steps.
    if i == 100:
        duration_est(tic)
