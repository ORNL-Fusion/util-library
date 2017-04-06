function psi = calc_psi_mpex(coil,current,r,z)
psi = r.*afield_circular_coilset(coil,current,r,z);