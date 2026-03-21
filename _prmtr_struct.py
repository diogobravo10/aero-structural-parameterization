from nastran.model import nastran_model
import nastran.cards as nastran_cards
import numpy as np
from _prmtr_aero_utils import parameterization as prmtr

def grid_update(bdf, grid_x1, le_gwt, te_gwt, k):

    if k > 0 and k <= 1: #positive moves stick aft (2TE)
        dx = prmtr(k, grid_x1, te_gwt) - grid_x1

    elif k < 0 and k >= -1: #negative moves stick forward (2LE)
        dx = prmtr(-k, grid_x1, le_gwt) - grid_x1

    else:
        return

    for grid_id in bdf.grid_cards:
        grid_card = bdf.grid_cards[grid_id]
        grid_card.x1 += dx


def geo_parameterization(geo_bdf, mass_wing_bdf, fuse_mass_bdf, aero_bdf, k_geo, k_mass_wing, k_mass_fuse):

    grid_wt = geo_bdf.grid_cards[119]
    caero_wt = aero_bdf.caero1_cards[6209001]

    eta = (grid_wt.x3 - caero_wt.z1)/(caero_wt.z4 - caero_wt.z1)
    le_gwt = caero_wt.x1 + eta*(caero_wt.x4 - caero_wt.x1)
    te_gwt = (caero_wt.x1 + caero_wt.x12) + eta*(caero_wt.x4+caero_wt.x43 - (caero_wt.x1+caero_wt.x12))

    grid_update(mass_wing_bdf, grid_wt.x1, le_gwt, te_gwt, k_mass_wing)
    grid_update(fuse_mass_bdf, grid_wt.x1, le_gwt, te_gwt, k_mass_fuse)
    grid_update(geo_bdf, grid_wt.x1, le_gwt, te_gwt, k_geo)  


def prop_parameterization(mat_bdf, k_E, k_G):

    # for pbar_id in out_bdf.pbar_cards:
    #     pbar_card = out_bdf.pbar_cards[pbar_id]
    #     pbar_card.a = pbar_card.a * k_A/per
    #     pbar_card.i1 = pbar_card.i1 * k_I1/per
    #     pbar_card.i2 = pbar_card.i2 * k_I2/per
    #     pbar_card.j = pbar_card.j * k_J/per
    #     pbar_card.k1 = pbar_card.k1 * k_K1/per
    #     pbar_card.k2 = pbar_card.k2 * k_K2/per

    for mat_id in mat_bdf.mat1_cards:
        mat_card = mat_bdf.mat1_cards[mat_id]
        mat_card.e *= k_E 
        mat_card.g *= k_G


def inertia_parameterization(fuse_bdf, k_IF):

    for mass_id in fuse_bdf.conm2_card:
        mass_card = fuse_bdf.conm2_card[mass_id]
        mass_card.i11 *= k_IF
        mass_card.i21 *= k_IF
        mass_card.i22 *= k_IF
        mass_card.i31 *= k_IF
        mass_card.i32 *= k_IF
        mass_card.i33 *= k_IF