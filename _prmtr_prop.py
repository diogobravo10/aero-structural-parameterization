from nastran.model import nastran_model
import nastran.cards as nastran_cards
import numpy as np


def parameterization(mat_bdf, k_E, k_G):

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

