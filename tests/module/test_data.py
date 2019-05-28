# -*- coding: utf-8 -*-


# Quaternary ligand binding to aromatic residues in the active-site gorge of acetylcholinesterase.
# Harel, M., Schalk, I., Ehret-Sabatier, L., Bouet, F., Goeldner, M., Hirth, C., Axelsen,
# P.H., Silman, I., Sussman, J.L. (1993) Proc.Natl.Acad.Sci.USA 90: 9031-9035
reference_1acj = {
    'pdb_id': '1acj',
    'ligand': 999, #'THA'
    'rings': [(5195, 5196, 5198, 5201, 5200, 5199), (5195, 5196, 5194, 5197, 5190, 5191),
             (5190, 5191, 5192, 5193, 5188, 5189)],
    'neighbours': [72, 80, 81, 84, 85, 117, 118, 119, 121, 122, 130, 199, 200, 201, 330, 334,
                  432, 436, 439, 440, 441, 442, 444],
    'hydrophobic': [84, 330, 334, 432, 439, 442],
    'ps_stacking': [84, 330],
    'ts_stacking': [],
    'hbond': [],
    'water_bridge': [{999, 634, 84}],
    'salt_bridge': [],
    'pi_cation': [],
    'halogen': []
}


# Crystallographic investigation of the role of aspartate 95 in the modulation of the redox potentials
# of Desulfovibrio vulgaris flavodoxin.
# McCarthy, A.A., Walsh, M.A., Verma, C.S., O'Connell, D.P., Reinhold, M., Yalloway, G.N., D'Arcy, D.,
# Higgins, T.M., Voordouw, G., Mayhew, S.G. (2002) Biochemistry 41: 10950-10962
reference_1aku = {
    'pdb_id': '1aku',
    'ligand': 150, # FMN
    'rings': [(1355, 1356, 1358, 1359, 1360, 1371), (1360, 1371, 1363, 1372, 1370, 1362),
             (1362, 1370, 1364, 1369, 1365, 1367)],
    'neighbours': [9, 10, 11, 12, 13, 14, 15, 16, 17, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                  68, 93, 94, 95, 96, 98, 99, 100, 101, 102, 103, 104, 126, 128, 130],
    'hydrophobic': [60],
    'ps_stacking': [],
    'ts_stacking': [98],
    'hbond': [59, 60],
    'water_bridge': [{94, 171, 150}, {63, 185, 150}, {100, 185, 150}],
    'salt_bridge': [],
    'pi_cation': [],
    'halogen': []
}


# Crystal structures of Paracoccus denitrificans aromatic amino acid aminotransferase: a substrate recognition
# site constructed by rearrangement of hydrogen bond network.
# Okamoto, A., Nakai, Y., Hayashi, H., Hirotsu, K., Kagamiyama, H. (1998) J.Mol.Biol. 280: 443-461
reference_1ay8 = {
    'pdb_id': '1ay8',
    'chain': 7,
    'ligand': 413, # PLP (one of two. chain 7)
    'rings': [(3666, 3667, 3669, 3671, 3673, 3675)],
    'neighbours': [70, 106, 107, 108, 109, 110, 111, 112, 139, 140, 143, 189, 192, 194, 222, 224,
                  225, 254, 255, 256, 257, 258, 263, 266, 267, 268, 295, 296, 297, 360],
    'hydrophobic': [194, 140, 224],
    'ps_stacking': [140],
    'ts_stacking': [],
    'hbond': [108, 109, 194, 257],
    'water_bridge': [],
    'salt_bridge': [258, 266],
    'pi_cation': [],
    'halogen': []
}


# Oxyanion-mediated inhibition of serine proteases.
# Presnell, S.R., Patil, G.S., Mura, C., Jude, K.M., Conley, J.M., Bertrand, J.A., Kam, C.M.,
# Powers, J.C., Williams, L.D. (1998) Biochemistry 37: 17068-17081
reference_1bju = {
    'pdb_id': '1bju',
    'ligand': 910, # GP6
    'rings': [(2341, 2340, 2339, 2338, 2337, 2342, 2341), (2327, 2329, 2331, 2333, 2334, 2335)],
    'neighbours': [57, 94, 96, 99, 102, 172, 183, 189, 190, 191, 192, 193, 194, 195, 213, 214, 215,
                  216, 217, 219, 220, 221, 224, 225, 226, 227, 228, 901], # 901 = SO4
    'hydrophobic': [99],
    'ps_stacking': [57],
    'ts_stacking': [],
    'hbond': [189, 190, 195, 219],
    'water_bridge': [{190, 513, 910}, {227, 513, 910}],
    'salt_bridge': [],
    'pi_cation': [],
    'halogen': []
}


# Interaction of a peptidomimetic aminimide inhibitor with elastase.
# Peisach, E., Casebier, D., Gallion, S.L., Furth, P., Petsko, G.A., Hogan Jr., J.C., Ringe, D.
# (1995) Science 269: 66-69
reference_1bma = {
    'pdb_id': '1bma',
    'ligand': 256, # 0QH
    'rings': [(2281, 2283, 2285, 2286, 2287, 2288), (2276, 2277, 2278, 2279, 2282, 2284)],
    'neighbours': [60, 98, 101, 102, 103, 104, 105, 108, 152, 153, 179, 182, 198, 199, 200,
                   201, 202, 203, 204, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 236,
                   237],
    'hydrophobic': [103, 104, 152, 182, 223, 224], # literature: 103, 223
    'ps_stacking': [],
    'ts_stacking': [],
    'hbond': [224, 200],
    'water_bridge': [],
    'salt_bridge': [],
    'pi_cation': [],
    'halogen': []
}


# Structure of acetylcholinesterase complexed with E2020 (Aricept): implications for the design of
# new anti-Alzheimer drugs. Kryger, G., Silman, I., Sussman, J.L.
# (1999) Structure Fold.Des. 7: 297-307
reference_1eve = {
    'pdb_id': '1eve',
    'ligand': 2001, # E20
    'rings': [(6349, 6350, 6351, 6352, 6353, 6354), (6352, 6353, 6355, 6356, 6357),
              (6359, 6360, 6361, 6362, 6363, 6364), (6366, 6367, 6368, 6369, 6370, 6371)],
    'neighbours': [70, 72, 80, 81, 84, 117, 118, 119, 121, 122, 130, 199, 200, 201, 279, 282, 286, 287,
                  288, 289, 290, 330, 331, 334, 440, 441, 444],
    'hydrophobic': [330, 331, 334, 279, 84],
    'ps_stacking': [84, 279],
    'ts_stacking': [],
    'hbond': [288, 70],
    'water_bridge': [],
    'salt_bridge': [],
    'pi_cation': [330],
    'halogen': []
}


# Crystal structures of the choline/acetylcholine substrate-binding protein ChoX from Sinorhizobium meliloti
# in the liganded and unliganded-closed states.
# Oswald, C., Smits, S.H., Hoing, M., Sohn-Bosser, L., Dupont, L., Le Rudulier, D., Schmitt, L., Bremer, E.
# (2008) J.Biol.Chem. 283: 32848-32859
reference_2reg = {
    'pdb_id': '2reg',
    'ligand': 1, # CHT
    'chain': 0,
    'rings': [],
    'neighbours': [43, 45, 46, 71, 88, 90, 93, 94, 118, 119, 152, 153, 156, 157, 158, 159, 183, 204, 205],
    'hydrophobic': [],
    'ps_stacking': [],
    'ts_stacking': [],
    'hbond': [119, 156],
    'water_bridge': [],
    'salt_bridge': [], # TODO: should be 45, not found because charged N1 + O6 have formal charge of 0
    'pi_cation': [43, 90, 119, 205],
    'halogen': []
}


# Crystal Structure of Poxvirus Thymidylate Kinase: An Unexpected Dimerization Has Implications for Antiviral Therapy
# Caillat, C., Topalis, D., Agrofoglio, L.A., Pochet, S., Balzarini, J., Deville-Bonne, D., Meyer, P.
# (2008) Proc.Natl.Acad.Sci.USA 105: 16900
reference_2w0s = {
    'pdb_id': '2w0s',
    'ligand': 1207, # BVP
    'chain': 1,
    'rings': [],
    'neighbours': [13, 14, 17, 18, 37, 38, 39, 41, 49, 52, 53, 61, 64, 65, 68, 69, 72, 92, 93, 94, 95, 97,
                   98, 99, 101, 102, 105, 116, 142, 143, 144, 145, 149, 150, 1205, 1206, 1208],
    'hydrophobic': [68, 101],
    'ps_stacking': [68],
    'ts_stacking': [],
    'hbond': [72, 101],
    'water_bridge': [],
    'salt_bridge': [], # TODO: should be 41 and 93. fix this.
    'pi_cation': [], # Res 68 if pication_amine_angle to 35 deg.
    'halogen': [65]
}
