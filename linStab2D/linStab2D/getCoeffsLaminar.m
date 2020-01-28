% Continuity
A0_11 = -(U + R .* dUdr + 1i .* m .* V + R .* dWdz) ./ R;
A0_12 = -(RHO + R .* dRHOdr) ./ R;
A0_13 = -1i ./ R .* RHO .* m;
A0_14 = -dRHOdz;
A0_15 = Z;
Ar_11 = -U;
Ar_12 = -RHO;
Ar_13 = Z;
Ar_14 = Z;
Ar_15 = Z;
Arr_11 = Z;
Arr_12 = Z;
Arr_13 = Z;
Arr_14 = Z;
Arr_15 = Z;
Az_11 = -W;
Az_12 = Z;
Az_13 = Z;
Az_14 = -RHO;
Az_15 = Z;
Azz_11 = Z;
Azz_12 = Z;
Azz_13 = Z;
Azz_14 = Z;
Azz_15 = Z;
Arz_11 = Z;
Arz_12 = Z;
Arz_13 = Z;
Arz_14 = Z;
Arz_15 = Z;
B_11 = 1i*I;
B_12 = Z;
B_13 = Z;
B_14 = Z;
B_15 = Z;
% Radial momentum
A0_21 = -(2 .* dUdr .* U .* R .* c2 + W .* dUdz .* R .* c2 + dWdz .* U .* R .* c2 + 1i .* m .* U .* V .* c2 + dTdr .* R - V .^ 2 .* c2 + U .^ 2 .* c2) ./ R ./ c2;
A0_22 = -((3 .* MU .* m .^ 2 + 6 .* RHO .* U .* Re .* R + 6 .* dRHOdr .* U .* Re .* R .^ 2 + 6 .* dUdr .* RHO .* Re .* R .^ 2 + 3 .* W .* dRHOdz .* Re .* R .^ 2 + 3 .* dWdz .* RHO .* Re .* R .^ 2 + 2 .* dmudT .* dTdr .* R + 3 .* 1i .* RHO .* m .* V .* Re .* R + 4 .* MU) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_23 = -((2 .* 1i .* dmudT .* dTdr .* m .* R - 6 .* RHO .* V .* Re .* R + 7 .* 1i .* MU .* m + 3 .* 1i .* RHO .* U .* m .* Re .* R) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_24 = -RHO .* dUdz - dRHOdz .* U;
A0_25 = -((2 .* R .^ 2 .* c2 .* dWdz .* d2mudT2 .* dTdr + 2 .* R .* c2 .* U .* d2mudT2 .* dTdr - 4 .* R .^ 2 .* c2 .* dUdr .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c2 .* dWdr .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c2 .* dUdz .* d2mudT2 .* dTdz + 3 .* 1i .* dmudT .* m .* V .* c2 - 3 .* 1i .* dmudT .* m .* dVdr .* R .* c2 + 3 .* Re .* R .^ 2 .* dRHOdr) ./ Re ./ R .^ 2 ./ c2) ./ 0.3e1;
Ar_21 = -(T + U .^ 2 .* c2) ./ c2;
Ar_22 = 0.2e1 ./ 0.3e1 .* (2 .* dmudT .* dTdr .* R - 3 .* RHO .* U .* Re .* R + 2 .* MU) ./ Re ./ R;
Ar_23 = 1i ./ Re ./ R .* MU .* m ./ 0.3e1;
Ar_24 = 0.1e1 ./ Re .* dmudT .* dTdz;
Ar_25 = -((2 .* dWdz .* dmudT .* R .* c2 + 2 .* U .* dmudT .* c2 - 4 .* dUdr .* dmudT .* R .* c2 + 3 .* RHO .* Re .* R) ./ Re ./ R ./ c2) ./ 0.3e1;
Arr_21 = Z;
Arr_22 = 0.4e1 ./ 0.3e1 .* MU ./ Re;
Arr_23 = Z;
Arr_24 = Z;
Arr_25 = Z;
Az_21 = -W .* U;
Az_22 = -(W .* RHO .* Re - dmudT .* dTdz) ./ Re;
Az_23 = Z;
Az_24 = -((2 .* dmudT .* dTdr + 3 .* RHO .* U .* Re) ./ Re) ./ 0.3e1;
Az_25 = dmudT .* (dWdr + dUdz) ./ Re;
Azz_21 = Z;
Azz_22 = MU ./ Re;
Azz_23 = Z;
Azz_24 = Z;
Azz_25 = Z;
Arz_21 = Z;
Arz_22 = Z;
Arz_23 = Z;
Arz_24 = MU ./ Re ./ 0.3e1;
Arz_25 = Z;
B_21 = 1i .* U;
B_22 = 1i .* RHO;
B_23 = Z;
B_24 = Z;
B_25 = Z;
% Azimuthal momentum
A0_31 = -(1i .* m .* V .^ 2 .* c2 + 2 .* V .* U .* c2 + W .* dVdz .* R .* c2 + dWdz .* V .* R .* c2 + dVdr .* U .* R .* c2 + V .* dUdr .* R .* c2 + 1i .* m .* T) ./ R ./ c2;
A0_32 = ((-6 .* RHO .* V .* Re .* R - 3 .* V .* dRHOdr .* Re .* R .^ 2 - 3 .* dVdr .* RHO .* Re .* R .^ 2 + 7 .* 1i .* MU .* m + 3 .* 1i .* dmudT .* dTdr .* m .* R) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_33 = -((3 .* MU + 3 .* dmudT .* dTdr .* R + 6 .* RHO .* U .* Re .* R + 3 .* W .* dRHOdz .* Re .* R .^ 2 + 3 .* dWdz .* RHO .* Re .* R .^ 2 + 3 .* dRHOdr .* U .* Re .* R .^ 2 + 3 .* RHO .* dUdr .* Re .* R .^ 2 + 4 .* MU .* m .^ 2 + 6 .* 1i .* RHO .* V .* m .* Re .* R) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_34 = -(RHO .* dVdz .* Re .* R + dRHOdz .* V .* Re .* R - 1i .* dmudT .* dTdz .* m) ./ Re ./ R;
A0_35 = -((3 .* 1i .* RHO .* m .* Re .* R + 3 .* R .* c2 .* V .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c2 .* dVdr .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c2 .* dVdz .* d2mudT2 .* dTdz - 4 .* 1i .* dmudT .* m .* U .* c2 + 2 .* 1i .* dmudT .* m .* dUdr .* R .* c2 + 2 .* 1i .* dmudT .* m .* dWdz .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.3e1;
Ar_31 = -V .* U;
Ar_32 = -((3 .* RHO .* V .* Re .* R - 1i .* MU .* m) ./ Re ./ R) ./ 0.3e1;
Ar_33 = (MU + dmudT .* dTdr .* R - RHO .* U .* Re .* R) ./ Re ./ R;
Ar_34 = Z;
Ar_35 = -dmudT .* (V - dVdr .* R) ./ Re ./ R;
Arr_31 = Z;
Arr_32 = Z;
Arr_33 = MU ./ Re;
Arr_34 = Z;
Arr_35 = Z;
Az_31 = -W .* V;
Az_32 = Z;
Az_33 = (dmudT .* dTdz - W .* RHO .* Re) ./ Re;
Az_34 = -((3 .* RHO .* V .* Re .* R - 1i .* MU .* m) ./ Re ./ R) ./ 0.3e1;
Az_35 = 0.1e1 ./ Re .* dVdz .* dmudT;
Azz_31 = Z;
Azz_32 = Z;
Azz_33 = MU ./ Re;
Azz_34 = Z;
Azz_35 = Z;
Arz_31 = Z;
Arz_32 = Z;
Arz_33 = Z;
Arz_34 = Z;
Arz_35 = Z;
B_31 = 1i .* V;
B_32 = Z;
B_33 = 1i .* RHO;
B_34 = Z;
B_35 = Z;
% Axial momentum
A0_41 = -(dTdz .* R + W .* U .* c2 + W .* dUdr .* R .* c2 + 2 .* dWdz .* W .* R .* c2 + dWdr .* U .* R .* c2 + 1i .* m .* V .* W .* c2) ./ R ./ c2;
A0_42 = -((2 .* dmudT .* dTdz + 3 .* W .* RHO .* Re + 3 .* dWdr .* RHO .* R .* Re + 3 .* W .* dRHOdr .* R .* Re) ./ R ./ Re) ./ 0.3e1;
A0_43 = -(1i .* m .* (2 .* dmudT .* dTdz + 3 .* W .* RHO .* Re) ./ R ./ Re) ./ 0.3e1;
A0_44 = -(MU .* m .^ 2 + RHO .* U .* R .* Re + RHO .* dUdr .* R .^ 2 .* Re + 2 .* dWdz .* RHO .* R .^ 2 .* Re + dRHOdr .* U .* R .^ 2 .* Re + 2 .* dRHOdz .* W .* R .^ 2 .* Re + 1i .* RHO .* V .* m .* R .* Re) ./ R .^ 2 ./ Re;
A0_45 = -((3 .* dRHOdz .* R .* Re + 2 .* dUdr .* d2mudT2 .* dTdz .* R .* c2 - 4 .* dWdz .* d2mudT2 .* dTdz .* R .* c2 - 3 .* dWdr .* d2mudT2 .* dTdr .* R .* c2 - 3 .* dUdz .* d2mudT2 .* dTdr .* R .* c2 + 2 .* U .* d2mudT2 .* dTdz .* c2 - 3 .* 1i .* dmudT .* m .* dVdz .* c2) ./ R ./ Re ./ c2) ./ 0.3e1;
Ar_41 = -W .* U;
Ar_42 = -((2 .* dmudT .* dTdz + 3 .* W .* RHO .* Re) ./ Re) ./ 0.3e1;
Ar_43 = Z;
Ar_44 = -(-MU - dmudT .* dTdr .* R + RHO .* U .* R .* Re) ./ R ./ Re;
Ar_45 = dmudT .* (dWdr + dUdz) ./ Re;
Arr_41 = Z;
Arr_42 = Z;
Arr_43 = Z;
Arr_44 = MU ./ Re;
Arr_45 = Z;
Az_41 = -(T + W .^ 2 .* c2) ./ c2;
Az_42 = ((MU + 3 .* dmudT .* dTdr .* R) ./ R ./ Re) ./ 0.3e1;
Az_43 = 1i ./ R ./ Re .* MU .* m ./ 0.3e1;
Az_44 = 0.2e1 ./ 0.3e1 .* (2 .* dmudT .* dTdz - 3 .* W .* RHO .* Re) ./ Re;
Az_45 = -((3 .* RHO .* R .* Re + 2 .* dUdr .* dmudT .* R .* c2 - 4 .* dWdz .* dmudT .* R .* c2 + 2 .* U .* dmudT .* c2) ./ R ./ Re ./ c2) ./ 0.3e1;
Azz_41 = Z;
Azz_42 = Z;
Azz_43 = Z;
Azz_44 = 0.4e1 ./ 0.3e1 .* MU ./ Re;
Azz_45 = Z;
Arz_41 = Z;
Arz_42 = MU ./ Re ./ 0.3e1;
Arz_43 = Z;
Arz_44 = Z;
Arz_45 = Z;
B_41 = 1i .* W;
B_42 = Z;
B_43 = Z;
B_44 = 1i .* RHO;
B_45 = Z;
% Energy
A0_51 = -((1i .* m .* V .^ 3 .* c2 + 2 .* cv .* T .* U .* c2 + V .^ 2 .* dWdz .* R .* c2 + 3 .* U .^ 2 .* dUdr .* R .* c2 + 3 .* W .^ 2 .* dWdz .* R .* c2 + V .^ 2 .* dUdr .* R .* c2 + U .^ 2 .* dWdz .* R .* c2 + W .^ 2 .* dUdr .* R .* c2 + 2 .* 1i .* V .* m .* T + 2 .* T .* U + U .^ 3 .* c2 + 2 .* 1i .* V .* m .* cv .* T .* c2 + 2 .* cv .* dTdz .* W .* R .* c2 + 2 .* dUdz .* U .* W .* R .* c2 + 2 .* cv .* T .* dWdz .* R .* c2 + 2 .* dVdz .* V .* W .* R .* c2 + 2 .* dVdr .* V .* U .* R .* c2 + 2 .* dWdr .* U .* W .* R .* c2 + 2 .* cv .* T .* dUdr .* R .* c2 + 2 .* cv .* dTdr .* U .* R .* c2 + 1i .* V .* m .* U .^ 2 .* c2 + 1i .* V .* m .* W .^ 2 .* c2 + 2 .* dTdz .* W .* R + 2 .* T .* dUdr .* R + 2 .* dTdr .* U .* R + V .^ 2 .* U .* c2 + W .^ 2 .* U .* c2 + 2 .* T .* dWdz .* R) ./ R ./ c2) ./ 0.2e1;
A0_52 = ((-6 .* T .* RHO .* Re .* R - 8 .* MU .* dWdz .* c2 .* R + 6 .* MU .* d2Udz2 .* R .^ 2 .* c2 - 6 .* T .* dRHOdr .* Re .* R .^ 2 - 6 .* dTdr .* RHO .* Re .* R .^ 2 + 8 .* MU .* d2Udr2 .* R .^ 2 .* c2 + 2 .* MU .* d2Wdrz .* R .^ 2 .* c2 - 6 .* U .* MU .* m .^ 2 .* c2 - 6 .* cv .* T .* RHO .* Re .* c2 .* R - 6 .* W .* dRHOdz .* U .* Re .* R .^ 2 .* c2 - 18 .* dUdr .* RHO .* U .* Re .* R .^ 2 .* c2 - 6 .* dWdz .* RHO .* U .* Re .* R .^ 2 .* c2 - 6 .* dUdz .* W .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dVdr .* V .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dWdr .* W .* RHO .* Re .* R .^ 2 .* c2 - 6 .* cv .* T .* dRHOdr .* Re .* R .^ 2 .* c2 - 6 .* cv .* dTdr .* RHO .* Re .* R .^ 2 .* c2 + 12 .* 1i .* dVdr .* MU .* m .* R .* c2 - 9 .* U .^ 2 .* RHO .* Re .* c2 .* R - 3 .* V .^ 2 .* RHO .* Re .* c2 .* R - 3 .* W .^ 2 .* RHO .* Re .* c2 .* R - 4 .* W .* dmudT .* dTdz .* c2 .* R - 8 .* U .* dmudT .* dTdr .* c2 .* R + 2 .* 1i .* V .* MU .* m .* c2 - 3 .* W .^ 2 .* dRHOdr .* Re .* R .^ 2 .* c2 - 9 .* dRHOdr .* U .^ 2 .* Re .* R .^ 2 .* c2 - 3 .* V .^ 2 .* dRHOdr .* Re .* R .^ 2 .* c2 + 6 .* dUdz .* dmudT .* dTdz .* R .^ 2 .* c2 + 6 .* dWdr .* dmudT .* dTdz .* R .^ 2 .* c2 + 8 .* dUdr .* dmudT .* dTdr .* R .^ 2 .* c2 - 4 .* dWdz .* dmudT .* dTdr .* R .^ 2 .* c2 + 6 .* 1i .* V .* dmudT .* dTdr .* m .* R .* c2 - 6 .* 1i .* U .* RHO .* V .* m .* Re .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.6e1;
A0_53 = -((-6 .* MU .* d2Vdr2 .* R .^ 2 .* c2 - 6 .* MU .* d2Vdz2 .* R .^ 2 .* c2 + 6 .* dVdr .* MU .* R .* c2 + 8 .* 1i .* dUdr .* MU .* m .* R .* c2 + 3 .* 1i .* m .* RHO .* W .^ 2 .* Re .* R .* c2 + 3 .* 1i .* m .* RHO .* U .^ 2 .* Re .* R .* c2 + 6 .* 1i .* m .* RHO .* cv .* T .* Re .* R .* c2 + 8 .* 1i .* m .* MU .* dWdz .* R .* c2 - 2 .* 1i .* U .* MU .* m .* c2 + 6 .* 1i .* m .* RHO .* T .* Re .* R + 4 .* 1i .* W .* dmudT .* dTdz .* m .* R .* c2 + 8 .* V .* MU .* m .^ 2 .* c2 + 12 .* V .* dmudT .* dTdr .* R .* c2 - 6 .* dVdz .* dmudT .* dTdz .* R .^ 2 .* c2 - 6 .* dVdr .* dmudT .* dTdr .* R .^ 2 .* c2 + 6 .* U .* RHO .* V .* Re .* R .* c2 + 6 .* dWdz .* RHO .* V .* Re .* R .^ 2 .* c2 + 6 .* U .* dRHOdr .* V .* Re .* R .^ 2 .* c2 + 6 .* W .* dRHOdz .* V .* Re .* R .^ 2 .* c2 + 6 .* dUdr .* RHO .* V .* Re .* R .^ 2 .* c2 + 6 .* dVdz .* W .* RHO .* Re .* R .^ 2 .* c2 + 6 .* dVdr .* U .* RHO .* Re .* R .^ 2 .* c2 + 4 .* 1i .* U .* dmudT .* dTdr .* m .* R .* c2 + 9 .* 1i .* m .* RHO .* V .^ 2 .* Re .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.6e1;
A0_54 = ((6 .* MU .* d2Wdr2 .* R .^ 2 .* c2 - 6 .* T .* dRHOdz .* Re .* R .^ 2 - 6 .* dTdz .* RHO .* Re .* R .^ 2 + 2 .* MU .* d2Udrz .* R .^ 2 .* c2 + 8 .* MU .* d2Wdz2 .* R .^ 2 .* c2 + 6 .* MU .* dWdr .* R .* c2 + 2 .* MU .* dUdz .* R .* c2 - 6 .* W .* MU .* m .^ 2 .* c2 - 18 .* dWdz .* RHO .* W .* Re .* R .^ 2 .* c2 - 6 .* dUdr .* RHO .* W .* Re .* R .^ 2 .* c2 - 6 .* U .* dRHOdr .* W .* Re .* R .^ 2 .* c2 - 6 .* cv .* dTdz .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dUdz .* U .* RHO .* Re .* R .^ 2 .* c2 - 6 .* cv .* T .* dRHOdz .* Re .* R .^ 2 .* c2 - 6 .* dVdz .* V .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dWdr .* U .* RHO .* Re .* R .^ 2 .* c2 - 6 .* U .* RHO .* W .* Re .* R .* c2 + 12 .* 1i .* m .* MU .* dVdz .* R .* c2 - 9 .* dRHOdz .* W .^ 2 .* Re .* R .^ 2 .* c2 - 3 .* V .^ 2 .* dRHOdz .* Re .* R .^ 2 .* c2 - 3 .* U .^ 2 .* dRHOdz .* Re .* R .^ 2 .* c2 + 6 .* dUdz .* dmudT .* dTdr .* R .^ 2 .* c2 + 6 .* dWdr .* dmudT .* dTdr .* R .^ 2 .* c2 - 4 .* dUdr .* dmudT .* dTdz .* R .^ 2 .* c2 + 8 .* dWdz .* dmudT .* dTdz .* R .^ 2 .* c2 - 4 .* U .* dmudT .* dTdz .* R .* c2 + 6 .* 1i .* V .* dmudT .* dTdz .* m .* R .* c2 - 6 .* 1i .* W .* RHO .* V .* m .* Re .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.6e1;
A0_55 = -((2 .* R .* c1 .* c2 .* U .^ 2 .* d2mudT2 .* dTdr + 3 .* R .* c1 .* c2 .* V .^ 2 .* d2mudT2 .* dTdr - 3 .* Re .* R .^ 2 .* c2 .* d2mudT2 .* dTdr .^ 2 - 3 .* Re .* R .^ 2 .* c2 .* d2mudT2 .* dTdz .^ 2 + 3 .* Re .* R .* c1 .* U .* RHO + 3 .* Re .* R .^ 2 .* c1 .* dUdr .* RHO + 3 .* Re .* R .^ 2 .* c1 .* U .* dRHOdr + 3 .* Re .* R .^ 2 .* c1 .* dWdz .* RHO + 3 .* Re .* R .^ 2 .* c1 .* W .* dRHOdz + 3 .* Re .* R .* c1 .* c2 .* U .* RHO .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* W .* dRHOdz .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* dWdz .* RHO .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* dUdr .* RHO .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* U .* dRHOdr .* cv - 3 .* R .^ 2 .* c1 .* c2 .* dUdz .* W .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c1 .* c2 .* dWdr .* W .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c1 .* c2 .* dUdz .* U .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c1 .* c2 .* dWdr .* U .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c1 .* c2 .* dVdz .* V .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c1 .* c2 .* dVdr .* V .* d2mudT2 .* dTdr + 2 .* R .^ 2 .* c1 .* c2 .* dUdr .* W .* d2mudT2 .* dTdz - 4 .* R .^ 2 .* c1 .* c2 .* dWdz .* W .* d2mudT2 .* dTdz + 2 .* R .* c1 .* c2 .* U .* W .* d2mudT2 .* dTdz - 4 .* R .^ 2 .* c1 .* c2 .* dUdr .* U .* d2mudT2 .* dTdr + 2 .* R .^ 2 .* c1 .* c2 .* dWdz .* U .* d2mudT2 .* dTdr + 3 .* 1i .* V .* RHO .* m .* Re .* R .* c1 - 1i .* U .* dmudT .* m .* V .* c1 .* c2 + 3 .* MU .* m .^ 2 .* Re .* c2 + 3 .* 1i .* V .* RHO .* cv .* m .* Re .* R .* c1 .* c2 - 3 .* 1i .* U .* dmudT .* m .* dVdr .* R .* c1 .* c2 - 3 .* 1i .* W .* dmudT .* m .* dVdz .* R .* c1 .* c2 + 2 .* 1i .* V .* dmudT .* m .* dUdr .* R .* c1 .* c2 + 2 .* 1i .* V .* dmudT .* m .* dWdz .* R .* c1 .* c2) ./ Re ./ R .^ 2 ./ c1 ./ c2) ./ 0.3e1;
Ar_51 = -(U .* (V .^ 2 .* c2 + W .^ 2 .* c2 + 2 .* T + U .^ 2 .* c2 + 2 .* cv .* T .* c2) ./ c2) ./ 0.2e1;
Ar_52 = ((16 .* dUdr .* MU .* c2 .* R - 9 .* U .^ 2 .* RHO .* Re .* c2 .* R - 3 .* V .^ 2 .* RHO .* Re .* c2 .* R - 3 .* W .^ 2 .* RHO .* Re .* c2 .* R - 6 .* T .* RHO .* Re .* R - 8 .* MU .* dWdz .* c2 .* R + 2 .* 1i .* V .* MU .* m .* c2 - 6 .* cv .* T .* RHO .* Re .* c2 .* R - 4 .* W .* dmudT .* dTdz .* c2 .* R + 8 .* U .* dmudT .* dTdr .* c2 .* R) ./ Re ./ c2 ./ R) ./ 0.6e1;
Ar_53 = ((6 .* dVdr .* MU .* R - 3 .* V .* MU + 1i .* U .* MU .* m + 3 .* V .* dmudT .* dTdr .* R - 3 .* U .* RHO .* V .* Re .* R) ./ Re ./ R) ./ 0.3e1;
Ar_54 = (2 .* dWdr .* MU .* R + W .* MU + 2 .* MU .* dUdz .* R + W .* dmudT .* dTdr .* R + U .* dmudT .* dTdz .* R - W .* RHO .* U .* Re .* R) ./ Re ./ R;
Ar_55 = ((6 .* dmudT .* dTdr .* Re .* R .* c2 - 3 .* Re .* R .* c1 .* U .* RHO - 2 .* U .^ 2 .* dmudT .* c1 .* c2 - 3 .* V .^ 2 .* dmudT .* c1 .* c2 + 3 .* MU .* Re .* c2 + 3 .* dUdz .* W .* dmudT .* R .* c1 .* c2 + 3 .* dWdr .* W .* dmudT .* R .* c1 .* c2 - 3 .* Re .* R .* c1 .* c2 .* U .* RHO .* cv + 3 .* dVdr .* V .* dmudT .* R .* c1 .* c2 + 4 .* dUdr .* U .* dmudT .* R .* c1 .* c2 - 2 .* dWdz .* U .* dmudT .* R .* c1 .* c2) ./ Re ./ R ./ c1 ./ c2) ./ 0.3e1;
Arr_51 = Z;
Arr_52 = 0.4e1 ./ 0.3e1 ./ Re .* U .* MU;
Arr_53 = 0.1e1 ./ Re .* V .* MU;
Arr_54 = 0.1e1 ./ Re .* W .* MU;
Arr_55 = 0.1e1 ./ c1 .* MU;
Az_51 = -(W .* (V .^ 2 .* c2 + W .^ 2 .* c2 + 2 .* T + U .^ 2 .* c2 + 2 .* cv .* T .* c2) ./ c2) ./ 0.2e1;
Az_52 = ((6 .* MU .* dUdz .* R + W .* MU + 6 .* dWdr .* MU .* R + 3 .* W .* dmudT .* dTdr .* R + 3 .* U .* dmudT .* dTdz .* R - 3 .* W .* RHO .* U .* Re .* R) ./ Re ./ R) ./ 0.3e1;
Az_53 = ((6 .* dVdz .* MU .* R + 1i .* W .* MU .* m + 3 .* V .* dmudT .* dTdz .* R - 3 .* W .* RHO .* V .* Re .* R) ./ Re ./ R) ./ 0.3e1;
Az_54 = ((16 .* MU .* dWdz .* c2 .* R - 3 .* V .^ 2 .* RHO .* Re .* c2 .* R - 9 .* W .^ 2 .* RHO .* Re .* c2 .* R - 3 .* U .^ 2 .* RHO .* Re .* c2 .* R - 8 .* U .* MU .* c2 - 8 .* dUdr .* MU .* c2 .* R - 6 .* T .* RHO .* Re .* R + 2 .* 1i .* V .* MU .* m .* c2 - 6 .* cv .* T .* RHO .* Re .* c2 .* R + 8 .* W .* dmudT .* dTdz .* c2 .* R - 4 .* U .* dmudT .* dTdr .* c2 .* R) ./ Re ./ c2 ./ R) ./ 0.6e1;
Az_55 = -((-6 .* dmudT .* dTdz .* Re .* R .* c2 + 3 .* W .* RHO .* Re .* R .* c1 - 3 .* dUdz .* U .* dmudT .* R .* c1 .* c2 - 3 .* dWdr .* U .* dmudT .* R .* c1 .* c2 - 3 .* dVdz .* V .* dmudT .* R .* c1 .* c2 + 3 .* W .* RHO .* cv .* Re .* R .* c1 .* c2 + 2 .* dUdr .* W .* dmudT .* R .* c1 .* c2 - 4 .* dWdz .* W .* dmudT .* R .* c1 .* c2 + 2 .* U .* W .* dmudT .* c1 .* c2) ./ Re ./ R ./ c1 ./ c2) ./ 0.3e1;
Azz_51 = Z;
Azz_52 = 0.1e1 ./ Re .* U .* MU;
Azz_53 = 0.1e1 ./ Re .* V .* MU;
Azz_54 = 0.4e1 ./ 0.3e1 ./ Re .* W .* MU;
Azz_55 = 0.1e1 ./ c1 .* MU;
Arz_51 = Z;
Arz_52 = 0.1e1 ./ Re .* W .* MU ./ 0.3e1;
Arz_53 = Z;
Arz_54 = 0.1e1 ./ Re .* U .* MU ./ 0.3e1;
Arz_55 = Z;
B_51 = (1i .* (V .^ 2 + W .^ 2 + U .^ 2 + 2 .* cv .* T)) ./ 0.2e1;
B_52 = 1i .* U .* RHO;
B_53 = 1i .* V .* RHO;
B_54 = 1i .* W .* RHO;
B_55 = 1i .* RHO .* cv;