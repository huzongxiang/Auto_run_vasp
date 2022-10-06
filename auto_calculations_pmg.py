#! /work/phy-liuqh/softwares/python/anaconda3/envs/opencv/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 23:09:18 2019

@author: HZX
"""

import os
import warnings

try:
    import pymatgen as pmg
except ImportError:
    pmg = None
    raise ImportError("Module 'pymatgen' does not exist. Executing 'conda install --channel conda-forge pymatgen' to install the module please.")

from settings import Settings
from incar_pmg import INCAR
from pbs import PBS
from wannier import WANNIER
import auto_functions_pmg as af

def run_auto_calculations():

    if pmg is not None: 
        warnings.filterwarnings('ignore')
        
        print('----- AUTO RUN VASP -----\n')
        print('\n----- initiate varirable for calculation -----\n')
        # get control parameters from default yaml file, you can change it for your demends
        paras_control = af.control_yaml()
        noncolinear = paras_control['noncolinear']
        ini_magm = paras_control['ini_magm']  # initial magmoms from poscar that includes magmoms info
        ini_zero = paras_control['ini_zero'] # initial magmoms with zero
        add_U = paras_control['add_U']
        scdm = paras_control['scdm']
        symprec = paras_control['symprec']
        angle_tolerance = paras_control['angle_tolerance']
        kpath_cutoff = paras_control['kpath_cutoff']  # default 'None' for kpath cutoff, effective values must be int
        wannierwin = paras_control['wannierwin']

        auto_Settings = Settings()
        element_U = auto_Settings.element_U
        path_pot = auto_Settings.path_pot
        auto_INCAR = INCAR()
        auto_PBS = PBS()
        auto_WANNIER = WANNIER()
        pbs = auto_PBS.pbs
        pbs_band = auto_PBS.pbs_band
        pbs_wannier = auto_PBS.pbs_wannier
        pbs_post_wannier =auto_PBS.pbs_post_wannier
        incar_relax = auto_INCAR.INCAR_relax

        if noncolinear:
            incar_scf = auto_INCAR.INCAR_scf
            incar_band = auto_INCAR.INCAR_band
            incar_wannier = auto_INCAR.INCAR_wannier
            wannier_win = auto_WANNIER.wannier_win
        else:
            incar_scf = auto_INCAR.INCAR_scf_nonsoc
            incar_band = auto_INCAR.INCAR_band_nonsoc
            incar_wannier = auto_INCAR.INCAR_wannier_nonsoc         
            wannier_win = auto_WANNIER.wannier_win_nonsoc
        rotation_matrices = auto_Settings.rotation_matrices

        elements = af.read_poscar()
        ratios = af.read_ratios()
        structure = af.structure_pmg()    
        sites_valences = af.calc_valences_pmg(structure)  
        
        print('\n ----- symmetry -----\n')
        crystal_type = af.crystal_info(structure)
        senior_or_inferior = af.is_senior_cry_sys(crystal_type)
        
        print('\n\n----- build VASP input files -----\n')
        af.build_potcar(elements, path_pot)
        print('build POTCAR successfully...')
        NBANDS = af.calc_NBANDS_pmg(ratios, pbs, noncolinear)
        LMAXMIX = af.set_LMAXMIX_pmg(structure)
        af.build_pbs(pbs)
        af.build_pbs(pbs_band)
        af.build_pbs(pbs_wannier)
        af.build_pbs(pbs_post_wannier)
        print('build pbs successfully...')
        af.line_mode(structure)
        print('build KPOINTS for band calculation successfully...')
                
        relax, scf, band, wannier = af.init_incar_pmg(incar_relax, incar_scf, incar_band, incar_wannier)
        U_tag, line_L, line_U = af.config_U_calc(structure, element_U)
        if noncolinear:
            if not ini_magm:
                saxis = af.define_SAXIS(structure, rotation_matrices, senior_or_inferior)
                SAXIS, MAGMOM = af.calc_MAGMOM_pmg(sites_valences, structure, saxis, noncolinear, ini_zero)
            else:
                SAXIS, MAGMOM = af.read_magm_from_poscar()
            if add_U:
                if U_tag:
                    para_incar = {'NBANDS': NBANDS, 'MAGMOM': MAGMOM, 'SAXIS': SAXIS, 'LMAXMIX': LMAXMIX,
                            'LDAU':U_tag, 'LDAUTYPE': 2,'LDAUL': line_L,'LDAUU': line_U}
                else:
                    para_incar = {'LDAU': U_tag, 'NBANDS': NBANDS, 'MAGMOM': MAGMOM, 'SAXIS': SAXIS, 'LMAXMIX': LMAXMIX}
            else:
                U_tag = False
                para_incar = {'NBANDS': NBANDS, 'MAGMOM': MAGMOM, 'SAXIS': SAXIS, 'LMAXMIX': LMAXMIX}
            relax, scf, band, wannier = af.custom_incar_pmg(relax, scf, band, wannier, U_tag, **para_incar)
        else:
            if add_U:
                if U_tag:
                    para = {'NBANDS': NBANDS,'LDAU': U_tag,'LDAUTYPE': 2,'LDAUL': line_L, 'LDAUU': line_U, 'LMAXMIX': LMAXMIX}
                else:
                    para = {'LMAXMIX': LMAXMIX}
                relax, scf, band, wannier = af.custom_incar_pmg(relax, scf, band, wannier, U_tag, **para)
            else:
                U_tag = False
                para = {'NBANDS':NBANDS, 'LMAXMIX': LMAXMIX}
                relax, scf, band, wannier = af.custom_incar_pmg(relax, scf, band, wannier, U_tag, **para)
        af.build_incar_pmg(relax, scf, band, wannier)
        print('build INCAR successfully...\n')
        
        print('\n----- build wannier90 input files -----\n')
        if wannierwin:
            projectors, num_wann = af.set_projectors(elements, structure)
            if scdm:
                af.build_wannier_win_with_scdm(structure, wannier_win, num_wann, NBANDS, noncolinear, symprec, angle_tolerance, kpath_cutoff)
            else:
                af.build_wannier_win_with_projectors(structure, wannier_win, projectors, num_wann, NBANDS, noncolinear, symprec, angle_tolerance, kpath_cutoff)
        
        print('\n----- ALL DONE! FINISHED! ------\n\n')
    else:
        warnings.warn("Warning: module 'pymatgen' is not found.")

if __name__ == '__main__':
    run_auto_calculations()