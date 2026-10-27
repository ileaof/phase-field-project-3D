module variaveis

double precision::               Tm, delta_H, E0, W0, M0, D,D_l,D_s, Cp_l,Cp_s,  ANI, DELTA_EPSLON, DELTA_T, T, D_sol, D_liq, D_soluto_C,dif_conc
double precision, allocatable :: phase_field_old(:,:),phase_field_new(:,:),XC(:,:), YC(:,:), TNEW(:,:),TOLD(:,:), DPP(:,:),conc_soluto_new(:,:), conc_soluto_old(:,:)
double precision, allocatable :: TETA_0(:,:,:)
double PRECISION, allocatable :: T_(:,:)
double precision              :: XX, YY,XL1,YL1, h, h_E, h_W, h_N, h_S, CSP, CSN, CSS, CSE, CSW, CLP, CLN, CLS, CLE, CLW, Cl_eq, Cs_eq, Vm, R
double precision::               dx, dy, PHIy_E, PHIx_E, TETA_E , EPSLON_E_LINHA, EPSLON_E, PHIy_W, Me, C0,InvDX,InvDY
double precision::               PHIx_W, TETA_W, EPSLON_W_LINHA, EPSLON_W, PHIy_N, PHIx_N, TETA_N, EPSLON_N_LINHA,Hphi
double precision::               EPSLON_N, PHIy_S, PHIx_S, TETA_S, EPSLON_S_LINHA, EPSLON_S, g_linha, h_linha, PHI, small, AA, BB, CC, DD, EE, FF, GG                
double precision::               PP,PN,PS,PE,PW,PNW,PSW,PNE,PSE,PPNEW, AAA,ruidos, FONTE, TNEW1, T0, TE, TP, TW, TS, TN, CPP, CPS, CPN, CPW, CPE
double precision::               GAMA_P, GAMA_S, GAMA_N, GAMA_W, GAMA_E, GAMAn, GAMAw, GAMAe, GAMAs, c_new, Thin_Interface, Lambda
integer ::                        i,j, Nx, Ny, it_max,it,E_R,NYF,n
integer ::                       WriteFileD_timeCounter
real            ::               ran, Ke,time, Sigma, Beta, ajuste,teta_1, ran_1, teta_2,ran_2 
real            ::               teta_3,ran_3,ran_0,av,D1,D2,TX

double precision             ::               stoptime   
character(len=10) ::             text
character(len=14) ::             text2
character(len=38) ::             text3
                                   


end module variaveis




  




