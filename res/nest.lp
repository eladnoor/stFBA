\* NEST *\
Minimize
OBJ: __dummy
Subject To
conc_adp_c_lower: dGf_adp_c >= -1924.96219836
conc_adp_c_upper: dGf_adp_c <= -1924.51047164
conc_amp_c_lower: dGf_amp_c >= -1061.55242456
conc_amp_c_upper: dGf_amp_c <= -1060.54783146
conc_atp_c_lower: dGf_atp_c >= -2779.7709367
conc_atp_c_upper: dGf_atp_c <= -2779.50989069
conc_co2_c_lower: dGf_co2_c >= -562.039676069
conc_co2_c_upper: dGf_co2_c <= -539.21989069
conc_co2_e_lower: dGf_co2_e >= -552.900010579
conc_co2_e_upper: dGf_co2_e <= -551.895417482
conc_h2o_c_lower: dGf_h2o_c >= -254.329739513
conc_h2o_c_upper: dGf_h2o_c <= -254.280183331
conc_h2o_e_lower: dGf_h2o_e >= -254.329739513
conc_h2o_e_upper: dGf_h2o_e <= -254.280183331
conc_h_c_lower: dGf_h_c >= -68.4593511389
conc_h_c_upper: dGf_h_c <= -17.1148365347
conc_h_e_lower: dGf_h_e >= -39.9595228928
conc_h_e_upper: dGf_h_e <= -39.9346199143
conc_nh4_e_lower: dGf_nh4_e >= -90.7447936683
conc_nh4_e_upper: dGf_nh4_e <= -90.6952374864
conc_o2_c_lower: dGf_o2_c >= -34.2296760694
conc_o2_c_upper: dGf_o2_c <= -11.4098906898
conc_o2_e_lower: dGf_o2_e >= -34.2296760694
conc_o2_e_upper: dGf_o2_e <= -22.8197823796
conc_pi_c_lower: dGf_pi_c >= -1130.32967607
conc_pi_c_upper: dGf_pi_c <= -1107.50989069
conc_pi_e_lower: dGf_pi_e >= -1107.53479367
conc_pi_e_upper: dGf_pi_e <= -1107.48523749
dGr_ACONTa: dGf_acon_C_c - dGf_cit_c + dGf_h2o_c <= -1e-06
dGr_ACONTb: - dGf_acon_C_c - dGf_h2o_c + dGf_icit_c <= -1e-06
dGr_AKGDH: - dGf_akg_c + dGf_co2_c - dGf_coa_c - dGf_nad_c + dGf_nadh_c
 + dGf_succoa_c <= -1e-06
dGr_ATPM: dGf_adp_c - dGf_atp_c - dGf_h2o_c + dGf_h_c + dGf_pi_c <= -1e-06
dGr_ATPS4r: - dGf_adp_c + dGf_atp_c + dGf_h2o_c + 3 dGf_h_c - 4 dGf_h_e
 - dGf_pi_c <= -1e-06
dGr_CO2t: dGf_co2_c - dGf_co2_e >= 1e-06
dGr_CS: - dGf_accoa_c + dGf_cit_c + dGf_coa_c - dGf_h2o_c + dGf_h_c
 - dGf_oaa_c <= -1e-06
dGr_CYTBD: dGf_h2o_c - 2 dGf_h_c + 2 dGf_h_e - 0.5 dGf_o2_c + dGf_q8_c
 - dGf_q8h2_c <= -1e-06
dGr_ENO: - dGf_2pg_c + dGf_h2o_c + dGf_pep_c <= -1e-06
dGr_FBA: dGf_dhap_c - dGf_fdp_c + dGf_g3p_c <= -1e-06
dGr_FUM: - dGf_fum_c - dGf_h2o_c + dGf_mal__L_c <= -1e-06
dGr_G6PDH2r: dGf_6pgl_c - dGf_g6p_c + dGf_h_c - dGf_nadp_c + dGf_nadph_c
 <= -1e-06
dGr_GAPD: dGf_13dpg_c - dGf_g3p_c + dGf_h_c - dGf_nad_c + dGf_nadh_c
 - dGf_pi_c <= -1e-06
dGr_GLCpts: dGf_g6p_c - dGf_glc__D_e - dGf_pep_c + dGf_pyr_c <= -1e-06
dGr_GLNS: dGf_adp_c - dGf_atp_c + dGf_gln__L_c - dGf_glu__L_c + dGf_h_c
 - dGf_nh4_c + dGf_pi_c <= -1e-06
dGr_GLUDy: dGf_akg_c - dGf_glu__L_c - dGf_h2o_c + dGf_h_c - dGf_nadp_c
 + dGf_nadph_c + dGf_nh4_c >= 1e-06
dGr_GND: - dGf_6pgc_c + dGf_co2_c - dGf_nadp_c + dGf_nadph_c + dGf_ru5p__D_c
 <= -1e-06
dGr_H2Ot: dGf_h2o_c - dGf_h2o_e >= 1e-06
dGr_ICDHyr: dGf_akg_c + dGf_co2_c - dGf_icit_c - dGf_nadp_c + dGf_nadph_c
 <= -1e-06
dGr_MDH: dGf_h_c - dGf_mal__L_c - dGf_nad_c + dGf_nadh_c + dGf_oaa_c <= -1e-06
dGr_NADH16: - 4 dGf_h_c + 3 dGf_h_e + dGf_nad_c - dGf_nadh_c - dGf_q8_c
 + dGf_q8h2_c <= -1e-06
dGr_NH4t: dGf_nh4_c - dGf_nh4_e <= -1e-06
dGr_O2t: dGf_o2_c - dGf_o2_e <= -1e-06
dGr_PDH: dGf_accoa_c + dGf_co2_c - dGf_coa_c - dGf_nad_c + dGf_nadh_c
 - dGf_pyr_c <= -1e-06
dGr_PFK: dGf_adp_c - dGf_atp_c - dGf_f6p_c + dGf_fdp_c + dGf_h_c <= -1e-06
dGr_PGI: dGf_f6p_c - dGf_g6p_c <= -1e-06
dGr_PGK: dGf_13dpg_c - dGf_3pg_c + dGf_adp_c - dGf_atp_c >= 1e-06
dGr_PGL: dGf_6pgc_c - dGf_6pgl_c - dGf_h2o_c + dGf_h_c <= -1e-06
dGr_PGM: - dGf_2pg_c + dGf_3pg_c >= 1e-06
dGr_PIt2r: dGf_h_c - dGf_h_e + dGf_pi_c - dGf_pi_e <= -1e-06
dGr_PPC: - dGf_co2_c - dGf_h2o_c + dGf_h_c + dGf_oaa_c - dGf_pep_c + dGf_pi_c
 <= -1e-06
dGr_PYK: - dGf_adp_c + dGf_atp_c - dGf_h_c - dGf_pep_c + dGf_pyr_c <= -1e-06
dGr_RPE: - dGf_ru5p__D_c + dGf_xu5p__D_c <= -1e-06
dGr_RPI: - dGf_r5p_c + dGf_ru5p__D_c >= 1e-06
dGr_SUCDi: dGf_fum_c - dGf_q8_c + dGf_q8h2_c - dGf_succ_c <= -1e-06
dGr_SUCOAS: dGf_adp_c - dGf_atp_c - dGf_coa_c + dGf_pi_c - dGf_succ_c
 + dGf_succoa_c >= 1e-06
dGr_TALA: dGf_e4p_c + dGf_f6p_c - dGf_g3p_c - dGf_s7p_c <= -1e-06
dGr_TKT1: dGf_g3p_c - dGf_r5p_c + dGf_s7p_c - dGf_xu5p__D_c <= -1e-06
dGr_TKT2: - dGf_e4p_c + dGf_f6p_c + dGf_g3p_c - dGf_xu5p__D_c <= -1e-06
dGr_TPI: - dGf_dhap_c + dGf_g3p_c <= -1e-06
Bounds
__dummy = 0
-1000000 <= dGf_13dpg_c <= 1000000
-1000000 <= dGf_2pg_c <= 1000000
-1000000 <= dGf_3pg_c <= 1000000
-1000000 <= dGf_6pgc_c <= 1000000
-1000000 <= dGf_6pgl_c <= 1000000
-1000000 <= dGf_accoa_c <= 1000000
-1000000 <= dGf_acon_C_c <= 1000000
-1000000 <= dGf_adp_c <= 1000000
-1000000 <= dGf_akg_c <= 1000000
-1000000 <= dGf_amp_c <= 1000000
-1000000 <= dGf_atp_c <= 1000000
-1000000 <= dGf_cit_c <= 1000000
-1000000 <= dGf_co2_c <= 1000000
-1000000 <= dGf_co2_e <= 1000000
-1000000 <= dGf_coa_c <= 1000000
-1000000 <= dGf_dhap_c <= 1000000
-1000000 <= dGf_e4p_c <= 1000000
-1000000 <= dGf_f6p_c <= 1000000
-1000000 <= dGf_fdp_c <= 1000000
-1000000 <= dGf_fum_c <= 1000000
-1000000 <= dGf_g3p_c <= 1000000
-1000000 <= dGf_g6p_c <= 1000000
-1000000 <= dGf_glc__D_e <= 1000000
-1000000 <= dGf_gln__L_c <= 1000000
-1000000 <= dGf_glu__L_c <= 1000000
-1000000 <= dGf_h2o_c <= 1000000
-1000000 <= dGf_h2o_e <= 1000000
-1000000 <= dGf_h_c <= 1000000
-1000000 <= dGf_h_e <= 1000000
-1000000 <= dGf_icit_c <= 1000000
-1000000 <= dGf_mal__L_c <= 1000000
-1000000 <= dGf_nad_c <= 1000000
-1000000 <= dGf_nadh_c <= 1000000
-1000000 <= dGf_nadp_c <= 1000000
-1000000 <= dGf_nadph_c <= 1000000
-1000000 <= dGf_nh4_c <= 1000000
-1000000 <= dGf_nh4_e <= 1000000
-1000000 <= dGf_o2_c <= 1000000
-1000000 <= dGf_o2_e <= 1000000
-1000000 <= dGf_oaa_c <= 1000000
-1000000 <= dGf_pep_c <= 1000000
-1000000 <= dGf_pi_c <= 1000000
-1000000 <= dGf_pi_e <= 1000000
-1000000 <= dGf_pyr_c <= 1000000
-1000000 <= dGf_q8_c <= 1000000
-1000000 <= dGf_q8h2_c <= 1000000
-1000000 <= dGf_r5p_c <= 1000000
-1000000 <= dGf_ru5p__D_c <= 1000000
-1000000 <= dGf_s7p_c <= 1000000
-1000000 <= dGf_succ_c <= 1000000
-1000000 <= dGf_succoa_c <= 1000000
-1000000 <= dGf_xu5p__D_c <= 1000000
End
