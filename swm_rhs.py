## RIGHT HAND SIDE OF THE EQUATIONS
def rhs(u,v,eta,k,omega):
    """ Set of equations:

    u_t = qhv - p_x + Fx + Mx(u,v) - bottom_friction
    v_t = -qhu - p_y + My(u,v)  - bottom_friction
    eta_t = -(uh)_x - (vh)_y

    with p = .5*(u**2 + v**2) + gh, the bernoulli potential
    and q = (v_x - u_y + f)/h the potential vorticity

    using the enstrophy and energy conserving scheme (Arakawa and Lamb, 1981) and
    a biharmonic lateral mixing term based on Shchepetkin and O'Brien (1996).
    """

    #TODO param[nu_B] is large, applying the biharmonic creates tiny values (as dx^4 is large)
    #TODO think about a way to avoid possibly involved rounding errors especially with single precision
    #TODO might be efficiently only possible for dx=dy
    #UPDATE does not seem to be a major issue.

    h = eta + param['H']

    h_u = ITu.dot(h)    # h on u-grid
    h_v = ITv.dot(h)    # h on v-grid
    h_q = ITq.dot(h)    # h on q-grid

    U = u*h_u    # volume fluxes: U on u-grid
    V = v*h_v    # and V on v-grid

    KE = IuT.dot(u**2) + IvT.dot(v**2)  # kinetic energy without .5-factor

    q = (f_q + Gvx.dot(v) - Guy.dot(u)) / h_q       # potential vorticity q
    p = .5*KE + param['g']*h            # Bernoulli potential p

    ## BOTTOM FRICTION: quadratic drag
    sqrtKE = np.sqrt(KE)
    bfric_u = param['c_D']*ITu.dot(sqrtKE)*u/h_u
    bfric_v = param['c_D']*ITv.dot(sqrtKE)*v/h_v

    ## ADVECTION
    # Sadourny, 1975 enstrophy conserving scheme.
    # adv_u = Iqu.dot(q)*Ivu.dot(V)
    # adv_v = -Iqv.dot(q)*Iuv.dot(U)

    # Arakawa and Lamb, 1981
    adv_u, adv_v = ALadvection(q,U,V)

    ## LATERAL MIXING OPERATOR
    # simple bi-harmonic mixing
    # bidiff_u = param['nu_B']*LLu.dot(u)
    # bidiff_v = param['nu_B']*LLv.dot(v)

    ## K-OMEGA MODEL
    # constants
    beta_star = 9/100.
    c_lim = 7/8.
    nu_ref = 135.
    sig_star = 3/5.
    sig_do = 1/8.
    sig = 1/2.
    alpha = 13/25.
    beta0 = 0.0708

    # precompute
    dudx = Gux.dot(u)
    dudy = G2uy.dot(u)
    dvdx = G2vx.dot(v)
    dvdy = Gvy.dot(v)

    dkdx = GTx.dot(k)
    dkdy = GTy.dot(k)
    domdx = GTx.dot(omega)
    domdy = GTx.dot(omega)

    # Sij = (S11,S12,S22), S12 = S21, Sij is symmetric do not store twice
    Sij = (dudx,0.5*(dudy+dvdx),dvdy)
    Sij_square = Sij[0]**2 + 2*IqT.dot(Sij[1]**2) + Sij[2]**2

    omega_bar = np.maximum(omega,c_lim*np.sqrt(2/beta_star*Sij_square))
    nu_T = k/omega_bar
    tij = (2*nu_T*Sij[0]-2/3.*k,2*nu_T*Sij[1],2*nu_T*Sij[2]-2/3.*k)

    # momentum diffusion term - harmonic
    diff_u = (GTx.dot(h*tij[0]) + Gqy.dot(h*tij[1])) / h_u
    diff_v = (Gqx.dot(h*tij[1]) - GTy.dot(h*tij[2])) / h_v

    # turbulent kinetic energy k, advection adv, stress term, dissipation diss, diffusion diff
    k_adv = -(IuT.dot(u*dkdx) + IvT.dot(v*dkdy))
    k_stress = tij[0]*dudx + IqT.dot(tij[1]*(dudy + dvdx)) + tij[2]*dvdy
    k_diss = -beta_star*k*omega

    k_visc = nu_ref + sig_star*k/omega  # the viscosity that appears within the diffusion of k
    k_diff = GTx.dot(ITu.dot(k_visc)*dkdx) + GTy.dot(ITv.dot(k_visc)*dkdy)

    # specific dissipation rate omega, naming as above
    om_adv = -(IuT.dot(u*domdx) + IvT.dot(v*domdy))
    om_stress = alpha*omega/k*k_stress
    om_diss = -beta*omega**2
    om_k = (sig_do*(IuT.dot(dkdx*domdx) + IvT.dot(dkdy*domdy)).clip(0,1e30)/omega

    om_visc = nu_ref + sig*k/omega
    om_diff = GTx.dot(ITu.dot(om_visc)*domdx) + GTy.dot(ITv.dot(om_visc)*domdy)

    ## RIGHT-HAND SIDE: ADD TERMS
    rhs_u = adv_u - GTx.dot(p) + Fx/h_u + diff_u - bfric_u
    rhs_v = adv_v - GTy.dot(p) + diff_v - bfric_v
    rhs_eta = -(Gux.dot(U) + Gvy.dot(V))
    rhs_k = k_adv + k_stress + k_diss + k_diff
    rhs_omega = om_adv + om_stress + om_diss + om_k + om_diff

    return rhs_u, rhs_v, rhs_eta, rhs_k, rhs_omega

def ALadvection(q,U,V):
    """ Arakawa and Lamb,1981 advection terms. See interpolation.py for further information. """

    AL1q = AL1.dot(q)
    AL2q = AL2.dot(q)

    adv_u = Seur.dot(ALeur.dot(q)*U) + Seul.dot(ALeul.dot(q)*U) +\
            Sau.dot(AL1q[indx_au]*V) + Sbu.dot(AL2q[indx_bu]*V) +\
            Scu.dot(AL2q[indx_cu]*V) + Sdu.dot(AL1q[indx_du]*V)

    adv_v = Spvu.dot(ALpvu.dot(q)*V) + Spvd.dot(ALpvd.dot(q)*V) -\
            Sav.dot(AL1q[indx_av]*U) - Sbv.dot(AL2q[indx_bv]*U) -\
            Scv.dot(AL2q[indx_cv]*U) - Sdv.dot(AL1q[indx_dv]*U)

    return adv_u, adv_v
