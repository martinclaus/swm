## RIGHT HAND SIDE OF THE EQUATIONS

import numpy as np

from global_var import param, gvar

def rhs(u,v,eta):
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

    h = eta + param['H']

    h_u = gvar['ITu'].dot(h)    # h on u-grid
    h_v = gvar['ITv'].dot(h)    # h on v-grid
    h_q = gvar['ITq'].dot(h)    # h on q-grid

    U = u*h_u    # volume fluxes: U on u-grid
    V = v*h_v    # and V on v-grid

    KE = gvar['IuT'].dot(u**2) + gvar['IvT'].dot(v**2)  # kinetic energy without .5-factor

    q = (gvar['f_q'] + gvar['Gvx'].dot(v) - gvar['Guy'].dot(u)) / h_q       # potential vorticity q
    p = .5*KE + param['g']*h            # Bernoulli potential p

    ## BOTTOM FRICTION: quadratic drag
    sqrtKE = np.sqrt(KE)
    bfric_u = param['c_D']*gvar['ITu'].dot(sqrtKE)*u/h_u
    bfric_v = param['c_D']*gvar['ITv'].dot(sqrtKE)*v/h_v

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

    # symmetric stress tensor S = (S11, S12, S12, -S11), store only S11, S12
    S = (gvar['Gux'].dot(u) - gvar['Gvy'].dot(v),gvar['G2vx'].dot(v) + gvar['G2uy'].dot(u))
    hS = (h*S[0],h_q*S[1])

    diff_u = (gvar['GTx']*hS[0] + gvar['Gqy']*hS[1]) / h_u
    diff_v = (gvar['Gqx']*hS[1] - gvar['GTy']*hS[0]) / h_v

    # biharmonic stress tensor R = (R11, R12, R12, -R11), store only R11, R12
    R = (gvar['Gux'].dot(diff_u) - gvar['Gvy'].dot(diff_v), gvar['G2vx'].dot(diff_v) + gvar['G2uy'].dot(diff_u))
    nuhR = (param['nu_B']*h*R[0],param['nu_B']*h_q*R[1])

    bidiff_u = (gvar['GTx'].dot(nuhR[0]) + gvar['Gqy'].dot(nuhR[1])) / h_u
    bidiff_v = (gvar['Gqx'].dot(nuhR[1]) - gvar['GTy'].dot(nuhR[0])) / h_v

    ## RIGHT-HAND SIDE: ADD TERMS
    rhs_u = adv_u - gvar['GTx'].dot(p) + gvar['Fx']/h_u - bidiff_u - bfric_u
    rhs_v = adv_v - gvar['GTy'].dot(p) - bidiff_v - bfric_v
    rhs_eta = -(gvar['Gux'].dot(U) + gvar['Gvy'].dot(V))

    return rhs_u, rhs_v, rhs_eta

def ALadvection(q,U,V):
    """ Arakawa and Lamb,1981 advection terms. See interpolation.py for further information. """

    AL1q = gvar['AL1'].dot(q)
    AL2q = gvar['AL2'].dot(q)

    adv_u = (gvar['Seur'].dot(gvar['ALeur'].dot(q)*U)
             + gvar['Seul'].dot(gvar['ALeul'].dot(q)*U)
             + gvar['Sau'].dot(AL1q[gvar['indx_au']]*V)
             + gvar['Sbu'].dot(AL2q[gvar['indx_bu']]*V)
             + gvar['Scu'].dot(AL2q[gvar['indx_cu']]*V)
             + gvar['Sdu'].dot(AL1q[gvar['indx_du']]*V))

    adv_v = (gvar['Spvu'].dot(gvar['ALpvu'].dot(q)*V)
             + gvar['Spvd'].dot(gvar['ALpvd'].dot(q)*V)
             - gvar['Sav'].dot(AL1q[gvar['indx_av']]*U)
             - gvar['Sbv'].dot(AL2q[gvar['indx_bv']]*U)
             - gvar['Scv'].dot(AL2q[gvar['indx_cv']]*U)
             - gvar['Sdv'].dot(AL1q[gvar['indx_dv']]*U))

    return adv_u, adv_v
