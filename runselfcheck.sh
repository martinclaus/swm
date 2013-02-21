
cd selfcheck
time .././model
cd output
set -x
cdo diff reference_000000000001_eta_out.nc new_000000000001_eta_out.nc
cdo diff reference_000000000001_psi_out.nc new_000000000001_psi_out.nc
cdo diff reference_000000000001_u_out.nc new_000000000001_u_out.nc
cdo diff reference_000000000001_v_out.nc new_000000000001_v_out.nc
cdo diff reference_000000000001_u_mean_out.nc new_000000000001_u_mean_out.nc
cdo diff reference_000000000001_v_mean_out.nc new_000000000001_v_mean_out.nc

