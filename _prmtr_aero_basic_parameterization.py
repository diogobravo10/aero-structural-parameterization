from _prmtr_aero_utils import *

def SA_parameterization(caero_card, caero_ref, n, k):
    
    x4n_original = caero_card.x4 + n*caero_card.x43

    x1n_ref = caero_ref.x1 + n*caero_ref.x12

    m_ref, b_ref = line_boundary(x1n_ref, caero_ref.z4, x1n_ref, caero_ref.z1)

    x4n_ub, x4n_lb = point_boundaries(x4n_original, caero_ref.z4, m_ref, b_ref)  #upper bound (1) is the original, and lower bound (0) is calculated from the reference

    x4n_new = parameterization(k, x4n_lb, x4n_ub)
    
    caero_card.x4 = x4n_new - n*caero_ref.x43

def WS_parameterization(caero_card, caero_ref, k):

    z4_original = caero_card.z4
    z1_ref = caero_ref.z1

    z4_ub, z4_lb = point_boundaries(z4_original, 0, 0, z1_ref)

    caero_card.z4 = parameterization(k, z4_lb, z4_ub)

def DI_parameterization(caero_card, caero_ref, k):

    y4_original = caero_card.y4
    y1_ref = caero_ref.y1

    y4_ub, y4_lb = point_boundaries(y4_original, 0, 0, y1_ref) 

    caero_card.y4 = parameterization(k, y4_lb, y4_ub)


def WC_parameterization(caero_card, caero_ref, k):

    x43_original = caero_card.x43
    x43_ref = 0

    x43_ub, x43_lb = point_boundaries(x43_original, 0, 0, x43_ref)

    caero_card.x43 = parameterization(k, x43_lb, x43_ub)

def TR_parameterization(caero_card, caero_ref1, caero_ref2, k):

    x43_original = caero_card.x43
    d = caero_card.z4 - caero_ref1.z1
    cr = caero_ref1.x12
    s = caero_ref2.z4 - caero_ref1.z1

    delta_ref = 0
    x43_ub, x43_lb = point_boundaries(x43_original, d, -cr/s, cr)

    caero_card.x43 = parameterization(k, x43_lb, x43_ub)




def NL_parameterization(caero_card, caero_ref, n, k): #nose_length

    x1n_original = caero_card.x1 + n*caero_card.x12

    x1n_ref = caero_ref.x1 + n*caero_ref.x12
    x4n_ref = caero_ref.x4 + n*caero_ref.x43

    m_ref, b_ref = line_boundary(x1n_ref, caero_ref.z1, x4n_ref, caero_ref.z4)
    x1n_ub, x1n_lb = point_boundaries(x1n_original, caero_card.z1, m_ref, b_ref)
    x1n_new = parameterization(k, x1n_lb, x1n_ub)
    


    d = caero_card.x1 - x1n_new
    caero_card.x12 = caero_card.x12 + d
    caero_card.x1 = x1n_new

def RL_parameterization(caero_card, caero_ref, n, k): #rear_length

    x1n_original = caero_card.x1 + n*caero_card.x12

    x1n_ref = caero_ref.x1 + n*caero_ref.x12
    x4n_ref = caero_ref.x4 + n*caero_ref.x43

    m_ref, b_ref = line_boundary(x1n_ref, caero_ref.z1, x4n_ref, caero_ref.z4)
    x1n_ub, x1n_lb = point_boundaries(x1n_original, caero_card.z1, m_ref, b_ref)
    x1n_new = parameterization(k, x1n_lb, x1n_ub)




    caero_card.x12 = (x1n_new - caero_card.x1)/n




# def BCL_parameterization(bdf, id, caero_ref, n, k):
#     bdf.grid_cards[id] = cards.grid_card()
#     bdf.cbeam_cards[id] = cards.cbeam_cards()






# from basic_functions import *

# def SA_parameterization(caero_card, caero_ref, n, k):
    
#     x1n_ref = caero_ref.x1 + n*caero_ref.x12
#     x4n_ref = caero_ref.x4 + n*caero_ref.x43

#     x4n_lb, x4n_ub = point_boundaries(x4n_ref, 0, 0, x1n_ref)

#     x4n = parameterization(k, x4n_lb, x4n_ub)

#     caero_card.x4 = x4n - n*caero_ref.x43

# def WS_parameterization(caero_card, k):

#     z4_lb, z4_ub = point_boundaries(caero_card.z4, 0, 0, caero_card.z1)

#     caero_card.z4 = parameterization(k, z4_lb, z4_ub)

# def NL_parameterization(caero_card, caero_ref, n, k): #nose_length

#     m_ref, b_ref = line_boundary(caero_ref.x1, caero_ref.z1, caero_ref.x4, caero_ref.z4)

#     x1_lb, x1_ub = point_boundaries(caero_card.x1, caero_card.z1, m_ref, b_ref)


#     x1_new = parameterization(k, x1_lb, x1_ub)
#     d = caero_card.x1 - x1_new

#     caero_card.x12 = caero_card.x12 + d
#     caero_card.x1 = x1_new

# def RL_parameterization(caero_card, caero_ref, n, k): #rear_length

#     x1n = caero_card.x1 + caero_card.x12
#     x1n_ref = caero_ref.x1 + caero_ref.x12

#     x1n_lb, x1n_ub = point_boundaries(x1n, 0, 0, x1n_ref)

#     caero_card.x12 = parameterization(k, x1n_lb, x1n_ub) - caero_card.x1 