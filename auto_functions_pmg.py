# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 20:08:32 2019

@author: HZX
"""

import os
import re
import copy
import logging
import numpy as np
import cv2
import yaml 

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry import site_symmetries
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

stream = logging.StreamHandler() 
stream.setLevel(logging.WARNING)

log_file = logging.FileHandler('info.log')
log_file.setLevel(logging.INFO)

logger.addHandler(stream)
logger.addHandler(log_file)


ModulePath =  os.environ.get('HOME')


def control_yaml(paras=None):
    yaml.warnings({'YAMLLoadWarning': False})
    setting_path = os.path.join(ModulePath, 'auto_run_vasp/SETTING.yaml')
    if os.path.exists(setting_path):
        with open(setting_path,'r', encoding='utf-8') as f:   
            cfg = f.read()
        datas = yaml.load_all(cfg)  
        for paras in datas:
            print(paras)
    else:
        print(setting_path, ' does not exist!')
    return paras


def read_poscar():
    elements = False
    if os.path.exists('POSCAR') and os.path.getsize('POSCAR') != 0:
        with open('POSCAR') as file:
            lines_poscar = file.readlines()
        elements = lines_poscar[5].split()
    else:
        logger.warning('The POSCAR is not Effective!')
    return elements


def structure_pmg():
    if os.path.exists('POSCAR'):
        structure = Structure.from_file('POSCAR')
    return structure


def read_ratios():
    if os.path.exists('POSCAR') and os.path.getsize('POSCAR') != 0:
        with open('POSCAR') as file:
            lines_poscar = file.readlines()
        ratios = lines_poscar[6].split()
    else:
        logger.warning('The POSCAR is not Effective!')
    return ratios


def copy_pot(path_pot, pot_file, element):
    path_element = path_pot + '/' + element
    path_element_sv = path_pot + '/' + element + '_sv'
    path_element_pv = path_pot + '/' + element + '_pv'
    if os.path.exists(path_element):
        os.system('cat ' + path_element + '/POTCAR >> '+ pot_file)
    elif os.path.exists(path_element_pv):
        os.system('cat ' + path_element_pv + '/POTCAR >> '+ pot_file)
    else:
        os.system('cat ' + path_element_sv + '/POTCAR >> '+ pot_file)


def build_potcar(elements, path_pot):
    if not os.path.exists('POTCAR') or os.path.getsize('POTCAR') == 0 :
        pot_file = os.path.join(os.getcwd(),'POTCAR')
        if elements is not False:
            for element in elements:
                copy_pot(path_pot, pot_file, element)
    else:
        logger.warning('The POTCAR is alreadly existed')


def line_mode(structure, intersection=60, symprec=0.5, angle_tolerance=1e-1, cutoff=None, kpath_cutoff=None):
    print('generating Kpath.')
    kpath = HighSymmKpath(structure=structure, symprec=symprec, angle_tolerance=angle_tolerance).kpath
    print(kpath)
    kpath_lines = []
    for path in kpath['path']:
        if not cutoff:
            n_path = len(path) - 1
        else:
            if cutoff > len(path) - 1:
                n_path = len(path) - 1
            else:
                n_path = cutoff               
        i = 0
        while i < n_path:
            j = i + 1
            diff = kpath['kpoints'][path[i]] - kpath['kpoints'][path[j]]
            if np.linalg.norm(diff) < 0.01:
                if j < n_path:
                    j += 1
            kpath_lines.append('{:>12.5f}{:>12.5f}{:>12.5f}'.format(
                *kpath['kpoints'][path[i]]))
            kpath_lines.append('{:>8s}{:s}\n'.format('! ', path[i]))
            kpath_lines.append('{:>12.5f}{:>12.5f}{:>12.5f}'.format(
                *kpath['kpoints'][path[j]]))
            kpath_lines.append('{:>8s}{:s}\n\n'.format('! ', path[j]))
            i = j

    if isinstance(kpath_cutoff, int):
        if kpath_cutoff%4 != 0:
            kpath_cutoff = kpath_cutoff - kpath_cutoff%4
        if kpath_cutoff > len(kpath_lines_win):
            kpath_cutoff = None

    with open('KPOINTS.band', 'w') as f:
        f.write('KPoints Along Paths\n')
        f.write('{}\n'.format(intersection))
        f.write('Line-mode\n')
        f.write('rec\n')
        f.writelines(kpath_lines[:kpath_cutoff])
#    print("Notice: overlook 'does not match' warning! It does not matter!\nbuild KPOINTS for band successfully...")


def build_pbs(pbs):
    # warning: this func can only be used in producing pbs in order lsf, band, wannier, post_wannier, otherwise will result in bugs!
    print('building pbs')
    pbs_file = 'job_lsf'
    if not os.path.exists(pbs_file):
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    elif not os.path.exists('job_band'):
        pbs_file = 'job_band'
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    elif not os.path.exists('job_wannier'):
        pbs_file = 'job_wannier'
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    elif not os.path.exists('job_post_wannier'):
        pbs_file = 'job_post_wannier'
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    else:
        logger.warning('The pbs is alreadly existed')


def produce_pbs(pbs, name):
    print('building pbs')
    pbs_file = 'job_lsf'
    if name == 'lsf':
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    elif name == 'band':
        pbs_file = 'job_band'
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    elif name == 'wannier':
        pbs_file = 'job_wannier'
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    elif name == 'post_wannier':
        pbs_file = 'job_post_wannier'
        with open(pbs_file,'w') as file:
            file.writelines(pbs)
    else:
        logger.warning('The pbs is alreadly existed')


def calc_NBANDS_pmg(ratios, pbs, noncolinear=True):
#    logging.warning("Notice!\nPBE poseudopotential 'does not match' error should be overlock!\nIt does not matter!")
    print('calculating NBANDS')
    NBANDS = 0
    if os.path.exists('POTCAR'):
        potcar = Potcar.from_file("POTCAR")
        zvals = {}
        for pot in potcar:
            symbol = pot.symbol  
            zval = pot.nelectrons
            zvals[symbol] = zval
        print('\nzval : ',zvals)
        for i,key in enumerate(zvals):
            print(key,ratios[i])
            NBANDS += int(ratios[i])*int(zvals[key])
        print('total num of eletrons : ', NBANDS)
    if noncolinear:
        NBANDS = NBANDS*2
    num_cores = int(pbs[6].split()[2])    
    if NBANDS%num_cores ==0:
        NBANDS = int(NBANDS)
    else:
        NBANDS = (int(NBANDS/num_cores)+1)*num_cores
    print('Set NBANDS: ', NBANDS, ' multiples of core: ', num_cores)
    return NBANDS


def set_LMAXMIX_pmg(structure):
    LMAXMIX = 2
    species = structure.species
    for specie in species:
        if specie.is_lanthanoid or specie.is_actinoid:
            LMAXMIX = int(6)
            break
        elif specie.is_transition_metal:
            LMAXMIX = int(4)
    return LMAXMIX


def init_incar_pmg(incar_relax, incar_scf, incar_band, incar_wannier):
    print('init_incar executes')
    if incar_relax:
        relax = Incar.from_string(incar_relax)
    if incar_scf:
        scf = Incar.from_string(incar_scf)
    if incar_band:
        band = Incar.from_string(incar_band)
    if incar_wannier:
        wannier = Incar.from_string(incar_wannier)
    logging.info('init_incar executes successfully')
    return relax, scf, band, wannier


def config_U_calc(structure, element_U):
    print('configure U')
    list_L = []
    list_U = []
    line_L= ' '
    line_U = ' '
    U_tag = False
    elements = structure.composition.elements
    for element in elements:
        if str(element) in element_U:
            list_L.append('2')
            list_U.append(element_U[str(element)])
            U_tag = True
        else:
            list_L.append('-1')
            list_U.append('0')
    if U_tag:
        for L in list_L:
            line_L = line_L + L + '  '
        for U in list_U:
            line_U = line_U + U + '  '
    return U_tag, line_L, line_U


def custom_incar_pmg(relax, scf, band, wannier, U_tag, **dicts):
    print('custom_incar executes')
    if relax and scf and band and wannier:
        if U_tag:
            for key, value in dicts.items():
                scf[key] = value
                band[key] = value
                wannier[key] = value
                if key in ['NBANDS','LDAU','LDAUTYPE','LDAUL', 'LDAUU', 'LMAXMIX']:
                    relax[key] = value
        else:
            for key, value in dicts.items():
                relax[key] = value
                scf[key] = value
                band[key] = value
                wannier[key] = value
    return relax, scf ,band, wannier
    

def build_incar_pmg(relax, scf, band, wannier):
    print('build_incar executes')
    if not os.path.exists('INCAR'):
        os.system('touch INCAR')
    if relax:
        Incar(relax).write_file('INCAR.relax')
    if scf:
        Incar(scf).write_file('INCAR.scf')
    if band:
        Incar(band).write_file('INCAR.band')
    if wannier:
        Incar(wannier).write_file('INCAR.wannier')


def set_projectors(elements, structure):
    print('generating projectors')
    dict_projectors = {}
    projectors = []
    stoich = structure.composition.get_el_amt_dict()
    num_wann = 0
    for element in elements:
        if Element(element).is_transition_metal and not Element(element).is_lanthanoid and not Element(element).is_actinoid:
            projector = element + ': s; p; d\n'
            num_wann_element = stoich[element]*9
        elif Element(element).is_lanthanoid or Element(element).is_actinoid:
            projector = element + ': f\n'
            num_wann_element = stoich[element]*7
        elif Element(element).is_alkali or Element(element).is_alkaline or Element(element) == Element('H'):
            projector = element + ': s\n'
            num_wann_element = stoich[element]*1
        else:
            projector = element + ': s; p\n'
            num_wann_element = stoich[element]*4
        dict_projectors[element] = projector
        projectors.append(projector)
        num_wann = num_wann + num_wann_element
    print('projectors:')
    for element in dict_projectors:
        print(dict_projectors[element].strip())
    return projectors,num_wann
            

def build_wannier_win_with_projectors(structure, wannier_win, projectors, num_wann, NBANDS, noncolinear=True,
                                        symprec=0.5,angle_tolerance=1e-1,cutoff=None,kpath_cutoff=None):
    print('generating Kpath for wannier.win.')
    kpath = HighSymmKpath(structure=structure,symprec=symprec,angle_tolerance=angle_tolerance).kpath
    if noncolinear:
        num_wann_line = 'num_wann = ' + str(int(num_wann*2)) + '\n'
    else:
        num_wann_line = 'num_wann = ' + str(int(num_wann)) + '\n'
    nbands_line = 'num_bands = ' + str(NBANDS) + '\n'
    kpath_lines_win = []
    for path in kpath['path']:
        if not cutoff:
            n_path = len(path) - 1
        else:
            if cutoff > len(path) - 1:
                n_path = len(path) - 1
            else:
                n_path = cutoff               
        i = 0
        while i < n_path:
            j = i + 1
            m_i = re.search(r'[0-9]',path[i])
            if m_i is not None:
                init = re.search(r'[A-Z]',path[i]).group() + re.search(r'[0-9]',path[i]).group()
            else:
                init = re.search(r'[A-Z]',path[i]).group()
            m_j= re.search(r'[0-9]',path[j])
            if m_j is not None:
                end = re.search(r'[A-Z]',path[j]).group() + re.search(r'[0-9]',path[j]).group()
            else:
                end = re.search(r'[A-Z]',path[j]).group()
            diff = kpath['kpoints'][path[i]] - kpath['kpoints'][path[j]]
            if np.linalg.norm(diff) < 0.01:
                if j < n_path:
                    j += 1
            kpath_lines_win.append('{:<3s}'.format(init))
            kpath_lines_win.append('{:<7.3f}{:<7.3f}{:<8.3f}'.format(
                *kpath['kpoints'][path[i]]))
            kpath_lines_win.append('{:<3s}'.format(end))
            kpath_lines_win.append('{:<7.3f}{:<7.3f}{:<7.3f}\n'.format(
                *kpath['kpoints'][path[j]]))
            i = j

    if isinstance(kpath_cutoff, int):
        if kpath_cutoff%4 != 0:
            kpath_cutoff = kpath_cutoff - kpath_cutoff%4
        if kpath_cutoff > len(kpath_lines_win):
            kpath_cutoff = None

    with open('wannier90.win', 'w') as f:
        f.write(num_wann_line)
        f.write(nbands_line)
        f.writelines(wannier_win)   
        f.write('#kpath = True\n')
        f.write('#kpath_task = curv+bands\n')
        f.write('#kpath_num_points = 1000\n')  
        f.write('bands_plot  =  true\n')
        f.write('begin kpoint_path\n')
        f.writelines(kpath_lines_win[:kpath_cutoff])
        f.write('end kpoint_path\n')
        f.write('bands_num_points 120\n')
        f.write('bands_plot_format gnuplot xmgrace\n')
        f.write('\nBegin Projections\n')
        f.writelines(projectors)
        f.write('End Projections\n')
#    print("Notice: overlook 'does not match' warning! It does not matter!\nbuild wannier90.win successfully...")


def build_wannier_win_with_scdm(structure, wannier_win, num_wann, NBANDS, 
                                symprec=0.5,angle_tolerance=1e-1,cutoff=None,kpath_cutoff=None):
    print('build wannier.win with scdm.')
    kpath = HighSymmKpath(structure=structure,symprec=symprec,angle_tolerance=angle_tolerance).kpath
    num_wann_line = 'num_wann = ' + str(int(num_wann*2)) + '\n'
    nbands_line = 'num_bands = ' + str(NBANDS) + '\n'
    kpath_lines_win = []
    for path in kpath['path']:
        if not cutoff:
            n_path = len(path) - 1
        else:
            if cutoff > len(path) - 1:
                n_path = len(path) - 1
            else:
                n_path = cutoff               
        i = 0
        while i < n_path:
            j = i + 1
            m_i = re.search(r'[0-9]',path[i])
            if m_i is not None:
                init = re.search(r'[A-Z]',path[i]).group() + re.search(r'[0-9]',path[i]).group()
            else:
                init = re.search(r'[A-Z]',path[i]).group()
            m_j= re.search(r'[0-9]',path[j])
            if m_j is not None:
                end = re.search(r'[A-Z]',path[j]).group() + re.search(r'[0-9]',path[j]).group()
            else:
                end = re.search(r'[A-Z]',path[j]).group()
            diff = kpath['kpoints'][path[i]] - kpath['kpoints'][path[j]]
            if np.linalg.norm(diff) < 0.01:
                if j < n_path:
                    j += 1
            kpath_lines_win.append('{:<3s}'.format(init))
            kpath_lines_win.append('{:<7.3f}{:<7.3f}{:<8.3f}'.format(
                *kpath['kpoints'][path[i]]))
            kpath_lines_win.append('{:<3s}'.format(end))
            kpath_lines_win.append('{:<7.3f}{:<7.3f}{:<7.3f}\n'.format(
                *kpath['kpoints'][path[j]]))
            i = j

    if isinstance(kpath_cutoff, int):
        if kpath_cutoff%4 != 0:
            if kpath_cutoff%2 == 0:
                kpath_cutoff = kpath_cutoff + 2
            else:
                kpath_cutoff = kpath_cutoff + 1
        if kpath_cutoff > len(kpath_lines_win):
            kpath_cutoff = None

    with open('wannier90.win', 'w') as f:
        f.write(num_wann_line)
        f.write(nbands_line)
        f.write('! Automatic generation of initial projections\n') 
        f.write('auto_projections = .true.\n')
        f.writelines(wannier_win)  
        f.write('#kpath = True\n')
        f.write('#kpath_task = curv+bands\n')
        f.write('#kpath_num_points = 1000\n')  
        f.write('bands_plot  =  true\n')
        f.write('begin kpoint_path\n')
        f.writelines(kpath_lines_win[:kpath_cutoff])
        f.write('end kpoint_path\n')
        f.write('bands_num_points 120\n')
        f.write('bands_plot_format gnuplot xmgrace\n')
#    print("Notice: overlook 'does not match' warning! It does not matter!\nbuild wannier90.win successfully...")


def calc_valences_pmg(structure):
    print('calcute the probable valences of elements in the composition')
    _structure = copy.deepcopy(structure)
    BV = BVAnalyzer(symm_tol=0.01)
    valence_of_elements = {}
    try:
        sites_valences = BV.get_oxi_state_decorated_structure(_structure)
        species = sites_valences.types_of_specie
        for i in species:
            valence_of_elements[i.element] = i.oxi_state
        print('\nprediction for valence of elements:')
        for key in valence_of_elements:
            print(key,':',valence_of_elements[key])
    except:
        _structure.add_oxidation_state_by_guess()
        sites_valences = _structure
    print('\nstructure of valences:\n',sites_valences,'\n')
    return sites_valences


def crystal_info(structure):
    analyer = SpacegroupAnalyzer(structure)
    crystal_type = analyer.get_crystal_system()
    print('crystal type : ', crystal_type)
    sym_str =analyer.get_symmetrized_structure()
    spg_info = sym_str.get_space_group_info()
    print('space group: ',spg_info[0],' space group number: ',spg_info[1],sym_str,'\n')
    data_symmetry = analyer.get_symmetry_dataset()
    for key,value in data_symmetry.items():
        print(key,':',value)
    return crystal_type


def is_senior_cry_sys(crystal_type):
    ### exclude triclinic and monoclinic crystal system ###
    high_sym_cry_system = ["orthorhombic", "tetragonal",
                      "trigonal", "hexagonal", "cubic"]
    if crystal_type in high_sym_cry_system:
        return True
    else:
        return False


def define_SAXIS(structure, rotation_matrices, high_sym):
    logging.info('calcute the SAXIS for structure')
    if high_sym:
        print('\ncrystal is senior crystal, initiate SAXIS by symmetry')
        analyer = SpacegroupAnalyzer(structure)
        convention_cell = analyer.get_conventional_standard_structure()
        if structure.num_sites != convention_cell.num_sites:
            wyckoffs = analyer.get_symmetry_dataset()['wyckoffs']
            sort_wyckoffs = copy.deepcopy(wyckoffs)
            sort_wyckoffs.sort()
            
            for i,wyck in enumerate(wyckoffs):
                if wyck == sort_wyckoffs[0]:
                    site_order = i
                    break
            
            site = structure[site_order].frac_coords
            operations = site_symmetries.get_site_symmetries(structure)[site_order]
        
            dict_mat = {}
            dict_p = {}
            vector_rot = []
            for i,ops in enumerate(operations): 
                matrix_rot = np.around(ops.rotation_matrix,3)
                det_rot = np.around(np.linalg.det(matrix_rot),2)
                dict_mat[i] = matrix_rot
                vec = cv2.Rodrigues(matrix_rot)[0]
                arc = np.linalg.norm(vec)
                angle = np.around(arc*180/np.pi,1)
                if arc - 0 < 0.01:
                    vector_rot.append(vec.swapaxes(0,1))
                    dict_p[i] = (vec.swapaxes(0,1),angle,det_rot,matrix_rot)
                else:
                    vector_rot.append(vec.swapaxes(0,1)/arc)
                    dict_p[i] = (vec.swapaxes(0,1)/arc,angle,det_rot,matrix_rot)
        
                
            # convention_cell = analyer.get_conventional_standard_structure()
            
            analyer_std = SpacegroupAnalyzer(convention_cell)
            wyckoffs_std = analyer_std.get_symmetry_dataset()['wyckoffs']
            sort_wyckoffs_std = copy.deepcopy(wyckoffs_std)
            sort_wyckoffs_std.sort()

            site_order_std = 0
            for i,wyck in enumerate(wyckoffs_std):
                if wyck == sort_wyckoffs_std[0] and np.all(convention_cell[i].frac_coords == site):
                    site_order_std = i
                    break   
            
            operations_std = site_symmetries.get_site_symmetries(convention_cell)[site_order_std]
            
            dict_mat_std = {}
            list_mat_std = []
            dict_std = {}
            vector_rot_std = [] 
            for i,ops in enumerate(operations_std): 
                matrix_rot = np.around(ops.rotation_matrix,3)
                det_rot = np.around(np.linalg.det(matrix_rot),2)
                dict_mat_std[i] = matrix_rot
                list_mat_std.append(matrix_rot.astype(int).tolist())
                vec=cv2.Rodrigues(matrix_rot)[0]
                arc=np.linalg.norm(vec)
                angle=np.around(arc*180/np.pi,1)
                if arc - 0 < 0.01:
                    vector_rot_std.append(vec.swapaxes(0,1))
                    dict_std[i]=(vec.swapaxes(0,1),angle,det_rot,matrix_rot)
                else:
                    vector_rot_std.append(vec.swapaxes(0,1)/arc)
                    dict_std[i]=(vec.swapaxes(0,1)/arc,angle,det_rot,matrix_rot)
            
            print('\nrotations')
            rots = {}            
            for i,mat in enumerate(list_mat_std):
                for rot in rotation_matrices:
                    if np.all((np.linalg.det(mat)*np.array(mat)).astype(int).tolist() == rotation_matrices[rot]):
                        rots[i] = rot
                        print("{:<6s} {:<38s} ---- matrix: {:<38s} vector: {}".
                            format(rot,str(rotation_matrices[rot]),str(mat),str(np.around(vector_rot_std[i],5))))
            print('\nrots:',rots, '   it represents all symmetry operations on the highest site.')        
            list_rot_order = []           
            for i,rot in rots.items():
                m = re.match(r'\d',rot)
                list_rot_order.append((i,m.group()))
                
            max_list = []
            max_order = list_rot_order[0]
            for i in list_rot_order:
                if i[1] > max_order[1]:
                    max_list = []
                    max_order = i
                    max_list.append(i)
                if i[1] == max_order[1]:
                    max_list.append(i)
            max_list = list(set(max_list))
            print('\nmax_list: ',max_list)

            order_sym = ''
            for tup in set(max_list):
                if np.all(dict_std[tup[0]][0] == [0,0,1]) \
                or np.all(dict_std[tup[0]][0] == [0,0,-1]):
                    order_sym = tup[0]
                    break
            if isinstance (order_sym,int):
                if np.all(np.isnan( vector_rot[order_sym][0]) == np.array([ True,  True,  True])) or np.all(vector_rot[order_sym][0] == np.array([0., 0., 0.])):
                    saxis = [0,0,1]
                else:   
                    saxis = vector_rot[order_sym][0]
            else:
                saxis = [0,0,1]
        else:
            print("cell is convetional, initiate SAXIS along 'z'")
            saxis = [0, 0, 1]
    else:
        print("crystal is inferior crystal, initiate SAXIS along 'z'")
        saxis = [0, 0, 1]

    print('\nSAIXS along direction : ',saxis,' \n')
    
    return saxis


def calc_MAGMOM_pmg(sites_valences, structure, saxis, noncolinear=True, ini_zero=False):
    print('calculating MAGMOM')
    if sites_valences: 
        Mag = CollinearMagneticStructureAnalyzer(sites_valences,overwrite_magmom_mode='replace_all')
    else:
        Mag = CollinearMagneticStructureAnalyzer(structure,overwrite_magmom_mode='replace_all')
    mag_species = Mag.magnetic_species_and_magmoms
    print('mag_species',mag_species)
          
    list_mag_species = []
    for _specie in mag_species:
        list_mag_species.append(_specie)
    _species = sites_valences.species
    species = list(set(_species))
    list_species = list(map(str,species))
    for specie in list_species:
        if specie not in list_mag_species:
            mag_species[specie] = 0.0
    
    print('\n',mag_species)
    
    if noncolinear:
        noncolinear_mag = {}
        for _sp,_mag in mag_species.items():
            sp = re.match(r'[A-Za-z]*',_sp).group()
            if not ini_zero:
                if _mag >= 4.95:
                    noncolinear_mag[sp] = [0, 0, _mag - 1.5]
                else:
                    noncolinear_mag[sp] = [0, 0, _mag]
            else:
                noncolinear_mag[sp] = [0, 0, 0]
    
    for _sp_ in noncolinear_mag:
        print(_sp_,' : ',noncolinear_mag[_sp_])
    
    Magmoments = []
    stoich = structure.composition.get_el_amt_dict()
    print(stoich)
    for sp in structure.species:
        Magmoments.append(noncolinear_mag[str(sp)])
    
    MAG_SAXIS = Magmom(moment=Magmoments,saxis=saxis)
    SAXIS = list(MAG_SAXIS.saxis.flatten())
    MAGMOM = list(MAG_SAXIS.moment.flatten())
    
    print('MAGMOM :', MAGMOM)
    
    return SAXIS, MAGMOM 


def read_magm_from_poscar():
    _MAGMOM = []
    if os.path.exists('POSCAR') and os.path.getsize('POSCAR') != 0:
        with open('POSCAR') as file:
            lines_poscar = file.readlines()
    for line in lines_poscar[8:]:
        _MAGMOM.append(list(map(float,line.strip().split()[4:])))
    MAGMOM = list(np.array(_MAGMOM).flatten())
    SAXIS = [ 0, 0, 1 ]
    return SAXIS, MAGMOM
    
    
    
    

        

            

