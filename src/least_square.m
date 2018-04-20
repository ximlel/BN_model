function d_U_sr = least_square(d_UL,d_UR,d_phi_sL,d_phi_sR,CA_LL,CA_Lsr,CA_RR,CA_Rsr)
%LEAST_SQUARE  least square method to compute \partial_x U in star region
%              for solid or gas phase
b=[CA_LL(2:4,2:4)*d_UL+(CA_LL(2:4,1)-CA_Lsr(2:4,1))*d_phi_sL;CA_RR(2:4,2:4)*d_UR-(CA_RR(2:4,1)-CA_Rsr(2:4,1))*d_phi_sR];
A=[CA_Lsr(2:4,2:4);CA_Rsr(2:4,2:4)];
d_U_sr=A\b;
end
