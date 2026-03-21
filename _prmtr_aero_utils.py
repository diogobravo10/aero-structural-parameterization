import numpy as np

########### BASIC FUNCTIONS ###########
def parameterization(k, lb, ub):
    return lb + (ub - lb) * k

def point_boundaries(y, x, m, b):
    return y, m * x + b     # ub, lb

def line_boundary(x1, y1, x4, y4): #projection
    m = (x1 - x4) / (y1 - y4)
    b = x1 - m * y1
    return m, b


########### PARAMETERIZATION FUNCTIONS ###########
def SA_parameterization(caero_card, caero_ref, n, k):
    
    x4n_original = caero_card.x4 + n*caero_card.x43

    x1n_ref = caero_ref.x1 + n*caero_ref.x12

    m_ref, b_ref = line_boundary(x1n_ref, caero_ref.z4, x1n_ref, caero_ref.z1)

    x4n_ub, x4n_lb = point_boundaries(x4n_original, caero_ref.z4, m_ref, b_ref)  #upper bound (1) is the original, and lower bound (0) is calculated from the reference

    x4n_new = parameterization(k, x4n_lb, x4n_ub)
    
    caero_card.x4 = x4n_new - n*caero_ref.x43

def DI_parameterization(caero_card, caero_ref, k):

    y4_original = caero_card.y4
    y1_ref = caero_ref.y1

    y4_ub, y4_lb = point_boundaries(y4_original, 0, 0, y1_ref) 

    caero_card.y4 = parameterization(k, y4_lb, y4_ub)

def WS_parameterization(caero_card, caero_ref, k):

    z4_original = caero_card.z4
    z1_ref = caero_ref.z1

    z4_ub, z4_lb = point_boundaries(z4_original, 0, 0, z1_ref)

    caero_card.z4 = parameterization(k, z4_lb, z4_ub)

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


########### CALCULATION FUNCTIONS ###########
def calculate_tan_sa(aero_bdf, n_wc, panels):
    """
    Calculate the sweep angle of the model wing_root to winglet.
    This function should return the sweep angle at n_wc in degrees.
    """
    
    z1 = aero_bdf.caero1_cards[panels[0]].z1
    z4 = aero_bdf.caero1_cards[panels[1]].z4

    x1 = aero_bdf.caero1_cards[panels[0]].x1 + n_wc * aero_bdf.caero1_cards[panels[0]].x12
    x4 = aero_bdf.caero1_cards[panels[1]].x4 + n_wc * aero_bdf.caero1_cards[panels[1]].x43

    tan_sa = (x4-x1)/(z4-z1)  # tan(sweep angle) = (x4 - x1)/(z4 - z1)
    
    return tan_sa, np.degrees(np.arctan(tan_sa)) 

def calculate_tan_di(aero_bdf, panels):
    """
    Calculate the dihedral angle of the model wing_root to winglet.
    This function should return the dihedral angle at n_wc in degrees.
    """
    
    z1 = aero_bdf.caero1_cards[panels[0]].z1
    z4 = aero_bdf.caero1_cards[panels[1]].z4

    y1 = aero_bdf.caero1_cards[panels[0]].y1
    y4 = aero_bdf.caero1_cards[panels[1]].y4

    tan_di = (y4-y1)/(z4-z1)  # tan(dihedral angle) = (y4 - y1)/(z4 - z1)
    
    return tan_di,  np.degrees(np.arctan(tan_di))

def calculate_ws(aero_bdf, panels):
    """
    Calculate the aspect ratio of the model.
    This function should return the aspect ratio.
    """
    
    z_r = aero_bdf.caero1_cards[panels[0]].z1    
    z_t = aero_bdf.caero1_cards[panels[1]].z4

    ws = z_t - z_r 
    b = 2 * z_t

    c_r = aero_bdf.caero1_cards[panels[0]].x12
    c_t = aero_bdf.caero1_cards[panels[1]].x43

    c_avg = (c_r + c_t) / 2

    ar = b / c_avg
    
    return ws, ar, b


########### CONVERSION FUNCTIONS ###########
def sa2k(sa, tan_sa_original):
    """
    Convert sweep angle to parameterization value.
    """
    
    tan_sa_0 = 0    # tan no sweep angle (k=0)
    tan_sa = np.tan(np.radians(sa))  # Convert sweep angle from degrees to radians and calculate tan(sweep angle)

    k = (tan_sa - tan_sa_0) / (tan_sa_original - tan_sa_0)

    return k

def di2k(di, tan_di_original):
    """
    Convert sweep angle to parameterization value.
    """
    
    tan_di_0 = 0    # tan no sweep angle (k=0)
    tan_di = np.tan(np.radians(di))  # Convert sweep angle from degrees to radians and calculate tan(sweep angle)

    k = (tan_di - tan_di_0) / (tan_di_original - tan_di_0)

    return k

def ar2k(aero_bdf, ar, ws_original, panels):
    """
    Convert aspect ratio to parameterization value.
    """
    
    fs = aero_bdf.caero1_cards[panels[0]].z1

    c_r = aero_bdf.caero1_cards[panels[0]].x12
    c_t = aero_bdf.caero1_cards[panels[1]].x43

    cavg = (c_r + c_t) / 2

    ws = ar * cavg /2 - fs # AR = 2* (ws + fs) / cavg  
    ws_0 = 0    # no aspect ratio (ar=0, k=0)

    k = (ws - ws_0) / (ws_original - ws_0)

    return k


########### REPOSITION FUNCTIONS ###########
def reposition_x1x4(caero_card, caero_card_next):
    dx = caero_card.x4 - caero_card_next.x1
    caero_card_next.x1 = caero_card_next.x1 + dx
    caero_card_next.x4 = caero_card_next.x4 + dx

def reposition_x2x3(caero_card, caero_card_next):
    dx = (caero_card.x4 + caero_card.x43) - (caero_card_next.x1 + caero_card_next.x12)
    caero_card_next.x1 = caero_card_next.x1 + dx
    caero_card_next.x4 = caero_card_next.x4 + dx

def reposition_y1y4(caero_card, caero_card_next):
    dx = caero_card.y4 - caero_card_next.y1
    caero_card_next.y1 = caero_card_next.y1 + dx
    caero_card_next.y4 = caero_card_next.y4 + dx

def reposition_z1z4(caero_card, caero_card_next):
    dz = caero_card.z4 - caero_card_next.z1
    caero_card_next.z1 = caero_card_next.z1 + dz
    caero_card_next.z4 = caero_card_next.z4 + dz

def reposition_x4x43(caero_card, caero_card_next):
    dx = caero_card.x1 - caero_card_next.x4
    caero_card_next.x4 = caero_card_next.x4 + dx
    caero_card_next.x43 = caero_card.x12