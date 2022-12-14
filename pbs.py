# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 21:41:49 2019

@author: HZX
"""

class PBS():
    
    def __init__(self):
        self.pbs = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q short\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'vasp_exe=/work/phy-liuqh/softwares/wannier/vasp.5.4.1/bin/vasp_ncl\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $vasp_exe"\n',
         '#---------------------------------------------------------------------\n',
         '# run vasp\n',
         '$COMMAND_std > vasp.out\n']
        
        self.pbs_band = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q short\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'vasp_exe=/work/phy-liuqh/softwares/wannier/vasp.5.4.1/bin/vasp_ncl\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $vasp_exe"\n',
         '#---------------------------------------------------------------------\n',
         '# run vasp\n',
         'ln ../CHG CHG\n',
         'ln ../CHGCAR CHGCAR\n',
         'ln ../WAVECAR WAVECAR\n',
         'cp INCAR.band INCAR\n',
         'cp KPOINTS.band KPOINTS\n',
         '$COMMAND_std > vasp.out\n']
        
        self.pbs_wannier = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q short\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'vasp_exe=/work/phy-liuqh/softwares/wannier/vasp.5.4.1/bin/vasp_ncl\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $vasp_exe"\n',
         '#---------------------------------------------------------------------\n',
         '# run vasp\n',
         'ln ../CHG CHG\n',
         'ln ../CHGCAR CHGCAR\n',
         'ln ../WAVECAR WAVECAR\n',
         'cp INCAR.wannier INCAR\n',
         '$COMMAND_std > vasp.out\n']

        self.pbs_post_wannier_up = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q short\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'wannier_exe=/work/phy-liuqh/softwares/wannier/wannier90.3/wannier90.x\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $wannier_exe wannier90.up"\n',
         '#---------------------------------------------------------------------\n',
         '# run post wannier90\n',
         '$COMMAND_std  > wannier90.up.wout\n']

        self.pbs_post_wannier_dn = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q short\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'wannier_exe=/work/phy-liuqh/softwares/wannier/wannier90.3/wannier90.x\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $wannier_exe wannier90.dn"\n',
         '#---------------------------------------------------------------------\n',
         '# run post wannier90\n',
         '$COMMAND_std  > wannier90.dn.wout\n']

        self.pbs_post_wannier = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q short\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'wannier_exe=/work/phy-liuqh/softwares/wannier/wannier90.3/wannier90.x\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $wannier_exe wannier90"\n',
         '#---------------------------------------------------------------------\n',
         '# run post wannier90\n',
         '$COMMAND_std\n']

        self.pbs_wt = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q debug\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'wt_exe=/work/phy-liuqh/softwares/wanniertools/wannier_tools_ANC/bin/wt.x\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $wt_exe"\n',
         '#---------------------------------------------------------------------\n',
         '# run wt.x\n',
         'ln ../wannier90_hr.dat wannier90_hr.dat\n',
         '$COMMAND_std\n']

        self.pbs_wt_up = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q debug\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'wt_exe=/work/phy-liuqh/softwares/wanniertools/wannier_tools_ANC/bin/wt.x\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $wt_exe"\n',
         '#---------------------------------------------------------------------\n',
         '# run wt.x\n',
         'ln ../wannier90.up_hr.dat wannier90_hr.dat\n',
         '$COMMAND_std\n']

        self.pbs_wt_dn = ['#!/bin/bash\n',
         '#BSUB -J name\n',
         '#BSUB -q debug\n',
         '#BSUB -e %J.err\n',
         '#BSUB -o %J.out\n',
         '#BSUB -R "span[ptile=40]"\n',
         '#BSUB -n 40\n',
         'hostfile=`echo $LSB_DJOB_HOSTFILE`\n',
         'NP=`cat $hostfile | wc -l`\n',
         'cd $LS_SUBCWD\n',
         '#\n',
         'source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux\n',
         'source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh\n',
         '#\n',
         'wt_exe=/work/phy-liuqh/softwares/wanniertools/wannier_tools_ANC/bin/wt.x\n',
         'COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP $wt_exe"\n',
         '#---------------------------------------------------------------------\n',
         '# run wt.x\n',
         'ln ../wannier90.dn_hr.dat wannier90_hr.dat\n',
         '$COMMAND_std\n']