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

def reposition_x12(caero_card, caero_card_next):
    caero_card_next.x12 = caero_card.x43