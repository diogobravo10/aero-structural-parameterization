from nastran.model import nastran_model
import numpy as np

import _prmtr_aero_utils as utils


######### SWEEP ANGLE PARAMETERIZATION ########## -> affects tan(angle)
# k = 1 sweptback, k = 0 no sweep, k = -1 sweptforward
# - affects x4
def sweep_angle(aero_bdf, sa, n_wc = 0.25):
    """
    Function to apply sweep angle parameterization to the wing CAERO1 cards.
    """

    caero_W_ids = [
        [6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
        [6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
    ]    

    tan_sa_original, sa_original = utils.calculate_tan_sa(aero_bdf, n_wc, [6203001, 6209001]) # tan sweep angle from original BWB model (k=100)
    ws_original, ar_original, b_original = utils.calculate_ws(aero_bdf, [6203001, 6209001])
    tan_di_original, di_original = utils.calculate_tan_di(aero_bdf, [6203001, 6209001])

    # sa2k - Convert sweep angle (sa) to parameterization value (k)
    if sa == 100:
        k = 1  # Original BWB model sweep angle
    else:    
        k = utils.sa2k(sa, tan_sa_original)

    # Parameterization and Reposition of Wing Panels
    for idx_W in range(2):
        for idx, id in enumerate(caero_W_ids[idx_W]):
            utils.SA_parameterization(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[id], n_wc, k) # - affects x4

            # Break out of the loop once last panel-winglet is repositioned by TE (do not change sweep of winglet)
            if id == caero_W_ids[idx_W][-2]:
                utils.reposition_x2x3(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]]) # reposition winglet by TE (do not change sweep of winglet)
                break
            else:              
                utils.reposition_x1x4(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])

    tan_sa_parameterized, sa_parameterized = utils.calculate_tan_sa(aero_bdf, n_wc, [6203001, 6209001])
    ws_parameterized, ar_parameterized, b_parameterized = utils.calculate_ws(aero_bdf, [6203001, 6209001])
    tan_di_parameterized, di_parameterized = utils.calculate_tan_di(aero_bdf, [6203001, 6209001])

    # print(f"Sweep angle of original BWB model: {sa_original:.2f} degrees")
    # print(f"Targeted Sweep Angle {sa:.2f} degrees")
    # print(f"Parameterized Sweep Angle {sa_parameterized:.2f} degrees")
    # print(f"Parameterization value k: {k:.2f}")
    
    # print(f"Wing Span of original BWB model: {b_original:.2f}")
    # print(f"Parameterized Wing Span: {b_parameterized:.2f}")
    # print(f"Aspect Ratio of original BWB model: {ar_original:.2f}")
    # print(f"Parameterized Aspect Ratio: {ar_parameterized:.2f}") 
    # print(f"Dihedral Angle of original BWB model: {di_original:.2f}")
    # print(f"Parameterized Dihedral Angle: {di_parameterized:.2f}\n")

    return sa_parameterized
    
########## DIHEDRAL PARAMETERIZATION ##########
# k = -1 no chord, k = 0 original chord, k = 1 double span
# - affects 
# mantain LE sweep angle
def dihedral_angle(aero_bdf, di, n_wc = 0.25):
    """
    Function to apply dihedral angle parameterization to the wing CAERO1 cards.
    """

    caero_W_ids = [
        [6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
        [6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
    ]    

    tan_di_original, di_original = utils.calculate_tan_di(aero_bdf, [6203001, 6209001]) # tan sweep angle from original BWB model (k=100)
    tan_sa_original, sa_original = utils.calculate_tan_sa(aero_bdf, n_wc, [6203001, 6209001]) # tan sweep angle from original BWB model (k=100)
    ws_original, ar_original, b_original = utils.calculate_ws(aero_bdf, [6203001, 6209001])

    # di2k - Convert sweep angle (sa) to parameterization value (k)
    if di == 100:
        k = 1  # Original BWB model sweep angle
    else:    
        k = utils.di2k(di, tan_di_original)

    for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
        for idx, id in enumerate(caero_W_ids[idx_W]):
                       
            utils.DI_parameterization(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[id], k)
            utils.reposition_y1y4(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])

            # Break out of the loop once last panel-winglet is repositioned by TE (do not change sweep of winglet)
            if id == caero_W_ids[idx_W][-2]:
                break
            
    tan_di_parameterized, di_parameterized = utils.calculate_tan_di(aero_bdf, [6203001, 6209001])  
    tan_sa_parameterized, sa_parameterized = utils.calculate_tan_sa(aero_bdf, n_wc, [6203001, 6209001])
    ws_parameterized, ar_parameterized, b_parameterized = utils.calculate_ws(aero_bdf, [6203001, 6209001])

    # print(f"Dihedral angle of original BWB model: {di_original:.2f} degrees")
    # print(f"Targeted Dihedral Angle {di:.2f} degrees")
    # print(f"Parameterized dihedral angle {di_parameterized:.2f} degrees")
    # print(f"Parameterization value k: {k:.2f}")

    # print(f"Wing Span of original BWB model: {b_original:.2f}")
    # print(f"Parameterized Wing Span: {b_parameterized:.2f}")
    # print(f"Aspect Ratio of original BWB model: {ar_original:.2f}")
    # print(f"Parameterized Aspect Ratio: {ar_parameterized:.2f}") 
    # print(f"Sweep Angle of original BWB model: {sa_original:.2f}")
    # print(f"Parameterized Sweep Angle: {sa_parameterized:.2f}\n")

    return di_parameterized


# ########## ASPECT RATIO (SPAN) PARAMETERIZATION ########## -> projection of span in xz plane!
# # k = -1 no span, k = 0 original span, k = 1 double span
# # - affects z4 and x4
# # to mantain sweep angle -> (x1-x4)_original/(x1-x4)_new = (z1-z4)_original/(z1-z4)_new 
def aspect_ratio(aero_bdf, ar, n_wc = 0.25):
    """
    Function to apply wing span parameterization to the wing CAERO1 cards.
    """

    caero_W_ids = [
        [6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
        [6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
    ]    

    ws_original, ar_original, b_original = utils.calculate_ws(aero_bdf, [6203001, 6209001]) # tan sweep angle from original BWB model (k=100)
    tan_sa_original, sa_original = utils.calculate_tan_sa(aero_bdf, n_wc, [6203001, 6209001]) 
    tan_di_original, di_original = utils.calculate_tan_di(aero_bdf, [6203001, 6209001])
    
    # ar2k - Convert aspect ratio (ar) to parameterization value (k)
    if ar == 100:
        k = 1  # Original BWB model saspect ratio
    else:    
        k = utils.ar2k(aero_bdf, ar, ws_original, [6203001, 6209001])


    for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
        for idx, id in enumerate(caero_W_ids[idx_W]):
                       
            utils.WS_parameterization(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[id], k)
            utils.reposition_z1z4(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])
            
            utils.SA_parameterization(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[id], n_wc, k)

            utils.DI_parameterization(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[id], k)
            utils.reposition_y1y4(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])
            
            # Break out of the loop once last panel-winglet is repositioned by TE (do not change sweep of winglet)
            if id == caero_W_ids[idx_W][-2]:
                utils.reposition_x2x3(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]]) # reposition winglet by TE (do not change sweep of winglet)
                break
            else:              
                utils.reposition_x1x4(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])

    ws_parameterized, ar_parameterized, b_parameterized = utils.calculate_ws(aero_bdf, [6203001, 6209001])
    tan_sa_parameterized, sa_parameterized = utils.calculate_tan_sa(aero_bdf, n_wc, [6203001, 6209001])
    tan_di_parameterized, di_parameterized = utils.calculate_tan_di(aero_bdf, [6203001, 6209001])


    # print(f"Aspect Ratio of original BWB model: {ar_original:.2f}")
    # print(f"Targeted Aspect Ratio: {ar:.2f}")
    # print(f"Parameterized Aspect Ratio: {ar_parameterized:.2f}")  
    # print(f"Parameterization value k: {k:.2f}")

    # print(f"Wing Span of original BWB model: {b_original:.2f}")
    # print(f"Parameterized Wing Span: {b_parameterized:.2f}")
    # print(f"Sweep Angle of original BWB model: {sa_original:.2f}")
    # print(f"Parameterized Sweep Angle: {sa_parameterized:.2f}")
    # print(f"Dihedral Angle of original BWB model: {di_original:.2f}")
    # print(f"Parameterized Dihedral Angle: {di_parameterized:.2f}\n")

    return ar_parameterized

# ########## Fuselage Length Parameterization ##########
# # k = 1 BWB, k = 0 flyingwing
# # - affects x1, x12
def fuselage_length(aero_bdf, k_fl, n_wc = 0.25):

    caero_F_ids = [
        [6202001, 6201001, 6200001],  # LHS Fuselage (centrebody -> fuselage wing root)
        [6002001, 6001001, 6000001]   # RHS Fuselage (centrebody -> fuselage wing root)
    ]
  
    caero_WR = [6203001, 6003001]  # Wing root panels, LHS Wing, RHS Wing

    for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
        for idx, id in enumerate(caero_F_ids[idx_W]):  #parameterize LHS Wing, RHS Wing
        
            caero = aero_bdf.caero1_cards[id]

            utils.RL_parameterization(caero, aero_bdf.caero1_cards[caero_WR[idx_W]], 1, k_fl/100)
            utils.NL_parameterization(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_WR[idx_W]], 0, k_fl/100)     

            # Break out of the loop once central panel is parameterized
            if idx + 1 == len(caero_F_ids[idx_W]):
                break

            # Reposition of Next Panel x4 x43 (Continuity)
            utils.reposition_x4x43(aero_bdf.caero1_cards[id], aero_bdf.caero1_cards[caero_F_ids[idx_W][idx + 1]])

    return k_fl


# bdf = nastran_model.from_file("wing-dlm-coarser.bdf")

# # Panels for parameterization
# caero_WR = [6203001, 6003001]  # Wing root panels, LHS Wing, RHS Wing

# caero_FR = [6202001, 6002001]  # Fuselage root panels, LHS Wing, RHS Wing

# caero_F_ids = [
#     [6202001, 6201001, 6200001],  # LHS Fuselage (centrebody -> fuselage wing root)
#     [6002001, 6001001, 6000001]   # RHS Fuselage (centrebody -> fuselage wing root)
# ]

# caero_W_ids = [
#     [6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
#     [6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
# ]


# ########## Fuselage Length Parameterization ##########
# # k = 0 BWB, k = 1 flyingwing
# # - affects x1, x12
# for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
#     for idx, id in enumerate(caero_F_ids[idx_W]):  #parameterize LHS Wing, RHS Wing
        
#         caero = bdf.caero1_cards[id]

#         RL_parameterization(caero, bdf.caero1_cards[caero_WR[idx_W]], 1, k_FL/per)
#         NL_parameterization(bdf.caero1_cards[id], bdf.caero1_cards[caero_WR[idx_W]], 0, k_FL/per)     

#         # Break out of the loop once central panel is parameterized
#         if idx + 1 == len(caero_F_ids[idx_W]):
#             break

#         # Reposition of Next Panel x4 x43 (Continuity)
#         reposition_x4x43(bdf.caero1_cards[id], bdf.caero1_cards[caero_F_ids[idx_W][idx + 1]])
        

# ########## ASPECT RATIO (SPAN) PARAMETERIZATION ########## -> projection of span in xz plane!
# # k = -1 no span, k = 0 original span, k = 1 double span
# # - affects z4 and x4
# # to mantain sweep angle -> (x1-x4)_original/(x1-x4)_new = (z1-z4)_original/(z1-z4)_new 
# for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
#         for idx, id in enumerate(caero_W_ids[idx_W]):
                       
#             WS_parameterization(bdf.caero1_cards[id], bdf.caero1_cards[id], k_WS/per)
#             reposition_z1z4(bdf.caero1_cards[id], bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])
            
#             SA_parameterization(bdf.caero1_cards[id], bdf.caero1_cards[id], n_wc, k_WS/per)
#             reposition_x1x4(bdf.caero1_cards[id], bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])

#             # Break out of the loop once last panel/winglet is repositioned (do not change span of winglet)
#             if idx + 1 == len(caero_W_ids[idx_W]) - 1:
#                 break

# # ########## ASPECT RATIO (Cavg) PARAMETERIZATION ##########
# # # k = -1 no chord, k = 0 original chord, k = 1 double span
# # # - affects x12 and x43
# # # mantain LE sweep angle

# caero1_WC = [
#     [6202001, 6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
#     [6002001, 6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
# ]

# for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
#         for idx, id in enumerate(caero1_WC[idx_W]):
                       
#             WC_parameterization(bdf.caero1_cards[id], bdf.caero1_cards[id], k_WC/per)
                        
#             # Break out of the loop once last panel/winglet is repositioned (do not change span of winglet)
#             if idx < len(caero1_WC[idx_W]) - 1:
#                 reposition_x12(bdf.caero1_cards[id], bdf.caero1_cards[caero1_WC[idx_W][idx + 1]])


# caero1_TR = [
#     [6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
#     [6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
# ]

# ########## TAPER RATIO (CHORD) PARAMETERIZATION ##########
# # k = -1 no chord, k = 0 original chord, k = 1 double span
# # - affects x12 and x43
# # mantain LE sweep angle

# for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
#         for idx, id in enumerate(caero1_TR[idx_W]):
                                             
#             if idx + 1 < len(caero1_TR[idx_W]):
#                 TR_parameterization(bdf.caero1_cards[id], bdf.caero1_cards[caero1_TR[idx_W][0]], bdf.caero1_cards[caero1_TR[idx_W][-2]], k_TR/per)
#                 reposition_x12(bdf.caero1_cards[id], bdf.caero1_cards[caero1_TR[idx_W][idx + 1]])




# ########## DIHEDRAL PARAMETERIZATION ##########
# # k = -1 no chord, k = 0 original chord, k = 1 double span
# # - affects x12 and x43
# # mantain LE sweep angle

# for idx_W in range(2):  #parameterize LHS Wing, RHS Wing
#         for idx, id in enumerate(caero_W_ids[idx_W]):
                       
#             # Break out of the loop once last panel/winglet is repositioned (do not change span of winglet)
#             if idx+1 < len(caero1_WC[idx_W]) - 1:
#                 DI_parameterization(bdf.caero1_cards[id], bdf.caero1_cards[id], k_DI/per)
#                 reposition_y1y4(bdf.caero1_cards[id], bdf.caero1_cards[caero_W_ids[idx_W][idx + 1]])




# ########## Wing cg Parameterization ########## -> xcoordinate of CG
# # k = -1 cg2TE, k = 0 original, k = 1 cg2LE
# # - affects x1
# for idx, id in enumerate(caero_FR):  #parameterize LHS Wing, RHS Wing
#     # Geometry of Wing root panel
#     X1w = bdf.caero1_cards[caero_WR[idx]].x1
#     X12w = bdf.caero1_cards[caero_WR[idx]].x12
#     X4w = bdf.caero1_cards[caero_WR[idx]].x4
#     Z1w = bdf.caero1_cards[caero_WR[idx]].z1
#     Z4w = bdf.caero1_cards[caero_WR[idx]].z4

#     X1f = bdf.caero1_cards[id].x1
#     X12f = bdf.caero1_cards[id].x12
#     X4f = bdf.caero1_cards[id].x4
#     Z1f = bdf.caero1_cards[id].z1
#     Z4f = bdf.caero1_cards[id].z4

#     if k_AC > 0:
#         m_fTE, b_fTE = line_boundary(X1f + X12f, Z1f, X1f + X12f, Z4f)
#         x3_lb, x3_ub = point_boundaries(bdf.caero1_cards[id].x4 + bdf.caero1_cards[id].x43, bdf.caero1_cards[id].z4, m_fTE, b_fTE)
#         bdf.caero1_cards[id].x4 = parameterization(k_AC, x3_lb, x3_ub) - bdf.caero1_cards[id].x43

#     if k_AC <= 0:
#         #m_fLE, b_fLE = line_boundary(X1f, Z1f, X4f, Z4f)
#         m_wLE, b_wLE = line_boundary(X1w, Z1w, X4w, Z4w)
#         x4_lb, x4_ub = point_boundaries(bdf.caero1_cards[id].x4, bdf.caero1_cards[id].z4, m_wLE, bdf.caero1_cards[id].x1 - m_wLE * bdf.caero1_cards[id].z1)
#         bdf.caero1_cards[id].x4 = parameterization(-k_AC, x4_lb, x4_ub)
    
#     dx = bdf.caero1_cards[caero_FR[idx]].x4 - bdf.caero1_cards[caero_WR[idx]].x1

#     # Reposition of Wing Panels x4 (Continuity)
#     for id_W in caero_W_ids[idx]:
#         bdf.caero1_cards[id_W].x1 = bdf.caero1_cards[id_W].x1 + dx
#         bdf.caero1_cards[id_W].x4 = bdf.caero1_cards[id_W].x4 + dx


# # Write to file
# bdf.to_file("wing-dlm-coarser-parameterized.bdf")

# # folder = "models/"
# # filename = f"wing-{k_FL}_{k_SA}_{k_WS}_{k_AC}.bdf"
# # filepath = folder + filename
# # bdf.to_file(filepath)  # Save the BDF file