import nastran.cards as nastran_cards

class normalized_grid_point:
    def __init__(self, grid_id, wing_side, caero1_id, chordwise_location=0.0, delta_thickness=0.0, spanwise_location=0.0):
        self.grid_id = grid_id
        self.wing_side = wing_side
        self.caero1_id = caero1_id
        self.chordwise_location = chordwise_location
        self.delta_thickness = delta_thickness
        self.spanwise_location = spanwise_location

def mapping(aero_bdf, aero_parameterized_bdf, stick_geo_mapped_bdf, stick_mass_mapped_bdf, stick_spline_mapped_bdf, stick_winglet_mapped_bdf, winglet_trim_mass_mapped_bdf):
  
    caero_W_ids = [
        [6200001, 6201001, 6202001, 6203001, 6204001, 6205001, 6206001, 6207001, 6208001, 6209001, 6210001],  # LHS Wing (wing root -> winglet)
        [6000001, 6001001, 6002001, 6003001, 6004001, 6005001, 6006001, 6007001, 6008001, 6009001, 6010001]   # RHS Wing (wing root -> winglet)
    ]

    map_unmap_grid_points(aero_bdf, aero_parameterized_bdf, stick_geo_mapped_bdf, caero_W_ids)
    map_unmap_grid_points(aero_bdf, aero_parameterized_bdf, stick_mass_mapped_bdf, caero_W_ids)
    map_unmap_grid_points(aero_bdf, aero_parameterized_bdf, stick_spline_mapped_bdf, caero_W_ids)    
    map_unmap_grid_points(aero_bdf, aero_parameterized_bdf, stick_winglet_mapped_bdf, caero_W_ids)
    map_unmap_grid_points(aero_bdf, aero_parameterized_bdf, winglet_trim_mass_mapped_bdf, caero_W_ids)


def read_CORD2R_create_tmp_grids(cord2r_ids, bdf):
    
    cord2r_grid_ids = []

    for id in cord2r_ids:
        cord2r_card = bdf.cord2r_cards[id]
        A1, A2, A3 = cord2r_card.a1, cord2r_card.a2, cord2r_card.a3
        B1, B2, B3 = cord2r_card.b1, cord2r_card.b2, cord2r_card.b3
        C1, C2, C3 = cord2r_card.c1, cord2r_card.c2, cord2r_card.c3

        # Create GRID cards for points A, B, and C
        grid_A_id = id * 10 + 1
        grid_B_id = id * 10 + 2
        grid_C_id = id * 10 + 3

        bdf.grid_cards[grid_A_id] = nastran_cards.grid_card(None, A1, A2, A3, id)
        bdf.grid_cards[grid_B_id] = nastran_cards.grid_card(None, B1, B2, B3, id)
        bdf.grid_cards[grid_C_id] = nastran_cards.grid_card(None, C1, C2, C3, id)

        cord2r_grid_ids.append(grid_A_id)
        cord2r_grid_ids.append(grid_B_id) 
        cord2r_grid_ids.append(grid_C_id)
    return cord2r_grid_ids

def update_CORD2R_pop_tmp_grids(cord2r_ids, cord2r_grid_ids, bdf):
    
    for id in cord2r_ids:

        grid_A_id = id * 10 + 1
        grid_B_id = id * 10 + 2
        grid_C_id = id * 10 + 3

        bdf.cord2r_cards[id].a1, bdf.cord2r_cards[id].a2, bdf.cord2r_cards[id].a3 = bdf.grid_cards[grid_A_id].x1, bdf.grid_cards[grid_A_id].x2, bdf.grid_cards[grid_A_id].x3
        bdf.cord2r_cards[id].b1, bdf.cord2r_cards[id].b2, bdf.cord2r_cards[id].b3 = bdf.grid_cards[grid_B_id].x1, bdf.grid_cards[grid_B_id].x2, bdf.grid_cards[grid_B_id].x3
        bdf.cord2r_cards[id].c1, bdf.cord2r_cards[id].c2, bdf.cord2r_cards[id].c3 = bdf.grid_cards[grid_C_id].x1, bdf.grid_cards[grid_C_id].x2, bdf.grid_cards[grid_C_id].x3
        bdf.grid_cards.pop(id * 10 + 1)
        bdf.grid_cards.pop(id * 10 + 2)
        bdf.grid_cards.pop(id * 10 + 3)


def map_unmap_grid_points(wing_dlm_coarser_bdf, wing_dlm_coarser_modified_bdf, bdf, caero_W_ids):
    
    grid_point_ids = list(bdf.grid_cards.keys())
    
    cord2r_ids = list(bdf.cord2r_cards.keys())
    
    if cord2r_ids:
        cord2r_grid_ids = read_CORD2R_create_tmp_grids(cord2r_ids, bdf)
        grid_point_ids.extend(cord2r_grid_ids)
    
    # Create a dictionary to hold normalized grid points
    normalized_grid_point_dict = {}

    # mapping structure to aerodynamics
    for i, grid_id in enumerate(grid_point_ids):
        x_gp = bdf.grid_cards[grid_id].x1
        y_gp = bdf.grid_cards[grid_id].x2
        z_gp = bdf.grid_cards[grid_id].x3
        found_panel = None
        wing_side = None
        chordwise_location = None
        delta_thickness = None
        spanwise_location = None

        # Determine wing_side from grid_id
        wing_side = 0 if bdf.grid_cards[grid_id].x3 > 0 else 1

        for caero_id in caero_W_ids[wing_side]:
            caero = wing_dlm_coarser_bdf.caero1_cards[caero_id]
            winglet_ls = wing_dlm_coarser_bdf.caero1_cards[caero_W_ids[0][-1]]
            winglet_rs = wing_dlm_coarser_bdf.caero1_cards[caero_W_ids[1][-1]]

            z1, z4 = wing_dlm_coarser_bdf.caero1_cards[caero_id].z1, wing_dlm_coarser_bdf.caero1_cards[caero_id].z4

            if (wing_side == 0 and z1 <= z_gp < z4) or (wing_side == 1 and z4 < z_gp <= z1): #point within aero panels 
                found_panel = caero_id
            elif (wing_side == 0 and z_gp >= winglet_ls.z4): #point outboard left winglet aero panel
                found_panel = caero_W_ids[0][-1]
            elif (wing_side == 1 and z_gp <= winglet_rs.z4): #point outboard rigth winglet aero panel
                found_panel = caero_W_ids[1][-1]

            if found_panel is not None:
                caero = wing_dlm_coarser_bdf.caero1_cards[found_panel]
                x1, x4 = caero.x1, caero.x4
                x12, x43 = caero.x12, caero.x43
                y1, y4 = caero.y1, caero.y4
                z1, z4 = caero.z1, caero.z4

                spanwise_location = (z_gp - z1) / (z4 - z1)
                delta_thickness = y_gp - (y1 + spanwise_location * (y4 - y1))
                x_le = x1 + spanwise_location * (x4 - x1)
                x_te = (x1 + x12) + spanwise_location * ((x4 + x43) - (x1 + x12))
                chord = x_te - x_le
                chordwise_location = (x_gp - x_le) / chord
                break


        normalized_grid_point_dict[grid_id] = normalized_grid_point(
            grid_id,
            wing_side,
            found_panel,
            chordwise_location,
            delta_thickness,
            spanwise_location
        )

    # mapping aerodynamics to structure
    for grid_id, norm_gp in normalized_grid_point_dict.items():
        caero_id = norm_gp.caero1_id
        wing_side = norm_gp.wing_side
        modified_caero = wing_dlm_coarser_modified_bdf.caero1_cards[caero_id]
        x1, x4 = modified_caero.x1, modified_caero.x4
        x12, x43 = modified_caero.x12, modified_caero.x43
        y1, y4 = modified_caero.y1, modified_caero.y4
        z1, z4 = modified_caero.z1, modified_caero.z4

        spanwise_location = norm_gp.spanwise_location
        chordwise_location = norm_gp.chordwise_location
        delta_thickness = norm_gp.delta_thickness

        z_gp = z1 + spanwise_location * (z4 - z1)
        y_gp = y1 + spanwise_location * (y4 - y1) + delta_thickness
        x_le = x1 + spanwise_location * (x4 - x1)
        x_te = (x1 + x12) + spanwise_location * ((x4 + x43) - (x1 + x12))
        chord = x_te - x_le
        x_gp = x_le + chordwise_location * chord

        bdf.grid_cards[grid_id].x1 = x_gp
        bdf.grid_cards[grid_id].x2 = y_gp
        bdf.grid_cards[grid_id].x3 = z_gp

    if cord2r_ids:
        cord2r_grid_ids = update_CORD2R_pop_tmp_grids(cord2r_ids, cord2r_grid_ids, bdf)


def CORD2R(G1_card, G2_card):
    A1, A2, A3 = G1_card.x1, G1_card.x2, G1_card.x3
    B1, B2, B3 = G2_card.x1, G2_card.x2, G2_card.x3
    C1, C2, C3 = G2_card.x1 + 1, G2_card.x2, G2_card.x3
    return nastran_cards.cord2r_card(None, A1, A2, A3, B1, B2, B3, C1, C2, C3)
