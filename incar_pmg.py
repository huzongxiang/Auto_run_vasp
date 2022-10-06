# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 21:43:25 2019

@author: HZX
"""

class INCAR():
    
    def __init__(self):
        self.INCAR_relax =  ' System    = Geo (VASP)\n \
        ####################################################################\n \
        Startparameter\n PREC      = Accurate  ! medium, high low\n \
        ISTART    = 0         ! job   : 0-new  1-cont  2-samecut\n \
        ICHARG    = 2         ! charge: 1-file 2-atom 10-const\n \
        ISPIN     = 2         ! spin polarized calculation?\n \
        #MAGMOM    = 4*4 20*0\n \
        #Selects the HSE06 hybrid function\n \
        #LHFCALC   = .TRUE. \n#HFSCREEN  = 0.2 \n \
        #AEXX      = 0.25  \n \
        #ALGO      = D \n \
        #TIME      = 0.4\n \
        #Electronic Relaxation\n \
        ENCUT     = 450.0     ! cut-off energy (eV)\n \
        ISMEAR    = 0         ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n SIGMA     = 0.1       ! broadening (eV) \n \
        NELM      = 500       ! maximum number of electronic SC steps\n#NELMIN    = 5         ! minimum number of electronic SC steps\n \
        #NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n \
        EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n \
        #IALGO     = 48        ! algorithm\n#LDIAG     = T         ! sub-space diagonalisation\n#LREAL     = Auto      ! real-space projection\n \
        #VOSKOWN = 1           !\n#IDIPOL   = 3\n#DIPOL    = 0.0 0.0 0.5\n#Mixing Parameters\n \
        #AMIX      = 0.05\n \
        #BMIX      = 0.01\n \
        #AMIX_MAG  = 0.05\n \
        #BMIX_MAG  = 0.0001\n \
        #MAXMIX    = 80\n \
        #Writing Parameters\n LWAVE     = F          !write WAVECAR\n \
        LCHARG    = F          !write CHGCAR\n \
        LVTOT     = F          !write LOCPOT\n \
        #RWIGS     = 0.77 1.17\n \
        #LORBIT    = 11\n \
        #LAECHG    = T\n \
        #LELF      = TRUE \n \
        #NBLOCK    = 1\n \
        #KBLOCK    = 50\n \
        KSPACING  = 0.3\n \
        KGAMMA    = TRUE\n \
        #LPARD     = T \n#IBAND     = 374 375 376 \n \
        #KPUSE     = 62 71 78   \n \
        #LSEPB     = T\n \
        #LSEPK     = T\n \
        ####################################################################\n \
        #Ionic relaxation\n ISIF      = 3          !stress and relaxation        \n \
        ISYM      = -1         !symmetry\n \
        EDIFFG    = -0.005     !stopping-criterion for ionic relaxation loop            \n \
        NSW       = 500        !number of steps for ionic relaxation loop\n \
        IBRION    = 2          !ionic relax: 0-MD 1-quasi-New 2-CG\n \
        #POTIM     = 0.2\n \
        #NFREE     = 3          ! initial steepest desc\n \
        ####################################################################\n \
        #LDA+U parameters\n#LDAU     = .TRUE.\n \
        #LDAUTYPE = 2\n#LDAUL    = -1    2    -1 \n \
        #LDAUU    =  0   2.50   0      / Na  Ru  O   \n \
        #LDAUJ    = 0.0  0.0 0.0 \n#LDAUPRINT= 2\n \
        #LMAXMIX  = 4\n\n \
        #ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n \
        LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n \
        NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n \
        #NCORE   = 4\n'

        self.INCAR_scf = ' System    = Scf (VASP)\n \
        ####################################################################\n \
        #Startparameter\n \
        PREC      = Accurate  ! medium, high low\n \
        ISTART    = 0         ! job   : 0-new  1-cont  2-samecut\n \
        ICHARG    = 2         ! charge: 1-file 2-atom 10-const\n \
        ISPIN     = 2         ! spin polarized calculation?\n \
        #MAGMOM    = 4*4 20*0\n \
        ####################################################################\n \
        #SOC parameter \n \
        LSORBIT    =.TRUE.\n \
        SAXIS     = 0 0 1 ! direction of the magnetic field \n \
        MAGMOM    = 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 \n \
        NBANDS = 40\n \
        GGA_COMPAT=.FALSE.\n \
        LORBMOM   =.TRUE.\n \
        #Selects the HSE06 hybrid function\n \
        #LHFCALC   = .TRUE. \n \
        #HFSCREEN  = 0.2 \n \
        #AEXX      = 0.25  \n \
        #ALGO      = D \n \
        #TIME      = 0.4\n \
        #Electronic Relaxation\n \
        ENCUT     = 450.0     ! cut-off energy (eV)\n \
        ISMEAR    = -5        ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n \
        SIGMA     = 0.1       ! broadening (eV) \n \
        NELM      = 1000       ! maximum number of electronic SC steps\n \
        #NELMIN    = 5         ! minimum number of electronic SC steps\n \
        #NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n \
        EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n \
        #IALGO     = 48        ! algorithm\n#LDIAG     = T         ! sub-space diagonalisation\n \
        #LREAL     = Auto      ! real-space projection\n \
        #VOSKOWN = 1           !\n \
        #IDIPOL   = 3\n \
        #DIPOL    = 0.0 0.0 0.5\n \
        #Mixing Parameters\n \
        #AMIX      = 0.05\n \
        #BMIX      = 0.01\n \
        #AMIX_MAG  = 0.05\n \
        #BMIX_MAG  = 0.0001\n \
        #MAXMIX    = 80\n \
        #Writing Parameters\n LWAVE     = T          !write WAVECAR\n \
        LCHARG    = T          !write CHGCAR\n \
        LVTOT     = F          !write LOCPOT\n \
        #RWIGS     = 0.77 1.17\n \
        LORBIT    = 11\n#LAECHG    = T\n \
        #LELF      = TRUE \n#NBLOCK    = 1\n \
        #KBLOCK    = 50\n \
        KSPACING  = 0.2\n \
        KGAMMA    = TRUE\n \
        #LPARD     = T \n \
        #IBAND     = 374 375 376 \n \
        #KPUSE     = 62 71 78   \n \
        #LSEPB     = T\n \
        #LSEPK     = T\n \
        ####################################################################\n \
        #Ionic relaxation\n \
        ISIF      = 3          !stress and relaxation        \n \
        #ISYM      = -1         !symmetry\n EDIFFG    = -0.01      !stopping-criterion for ionic relaxation loop            \n \
        NSW       = 0          !number of steps for ionic relaxation loop\n \
        IBRION    = -1         !ionic relax: 0-MD 1-quasi-New 2-CG\n \
        #POTIM     = 0.2\n \
        #NFREE     = 3          ! initial steepest desc\n \
        ####################################################################\n \
        #LDA+U parameters\n#LDAU     = .TRUE.\n \
        #LDAUTYPE = 2\n \
        #LDAUL    = -1    2    -1 \n \
        #LDAUU    =  0   2.50   0      / Na  Ru  O   \n \
        #LDAUJ    = 0.0  0.0 0.0 \n \
        #LDAUPRINT= 2\n \
        #LMAXMIX  = 4\n\n \
        #ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n \
        LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n \
        NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n \
        #NCORE   = 4\n'

        self.INCAR_band = ' System    = Band (VASP)\n \
        ####################################################################\n \
        #Startparameter\n \
        PREC      = Accurate  ! medium, high low\n \
        ISTART    = 1         ! job   : 0-new  1-cont  2-samecut\n \
        ICHARG    = 11        ! charge: 1-file 2-atom 10-const\n \
        ISPIN     = 2         ! spin polarized calculation?\n \
        #MAGMOM    = 4*4 20*0\n \
        ####################################################################\n \
        #SOC parameter \n \
        LSORBIT    =.TRUE.\n \
        SAXIS     = 0 0 1 ! direction of the magnetic field \n \
        MAGMOM    = 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 \n \
        NBANDS = 40\n \
        GGA_COMPAT=.FALSE.\n \
        LORBMOM   =.TRUE.\n \
        #Selects the HSE06 hybrid function\n \
        #LHFCALC   = .TRUE. \n \
        #HFSCREEN  = 0.2 \n \
        #AEXX      = 0.25  \n \
        #ALGO      = D \n \
        #TIME      = 0.4\n \
        #Electronic Relaxation\n \
        ENCUT     = 450.0     ! cut-off energy (eV)\n \
        ISMEAR    = 0         ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n \
        SIGMA     = 0.1       ! broadening (eV) \n \
        NELM      = 500       ! maximum number of electronic SC steps\n \
        #NELMIN    = 5         ! minimum number of electronic SC steps\n \
        #NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n \
        EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n \
        #IALGO     = 48        ! algorithm\n \
        #LDIAG     = T         ! sub-space diagonalisation\n \
        #LREAL     = Auto      ! real-space projection\n \
        #VOSKOWN = 1           !\n \
        #IDIPOL   = 3\n \
        #DIPOL    = 0.0 0.0 0.5\n \
        #Mixing Parameters\n \
        #AMIX      = 0.05\n \
        #BMIX      = 0.01\n \
        #AMIX_MAG  = 0.05\n \
        #BMIX_MAG  = 0.0001\n \
        #MAXMIX    = 80\n \
        #Writing Parameters\n \
        LWAVE     = T          !write WAVECAR\n \
        LCHARG    = F          !write CHGCAR\n \
        LVTOT     = F          !write LOCPOT\n \
        #RWIGS     = 0.77 1.17\n \
        LORBIT    = 11\n \
        #LAECHG    = T\n \
        #LELF      = TRUE \n \
        #NBLOCK    = 1\n \
        #KBLOCK    = 50\n \
        #KSPACING  = 0.3\n \
        #KGAMMA    = TRUE\n \
        #LPARD     = T \n \
        #IBAND     = 374 375 376 \n \
        #KPUSE     = 62 71 78   \n \
        #LSEPB     = T\n \
        #LSEPK     = T\n####################################################################\n \
        #Ionic relaxation\n \
        ISIF      = 3          !stress and relaxation        \n \
        #ISYM      = -1         !symmetry\n \
        EDIFFG    = -0.01      !stopping-criterion for ionic relaxation loop            \n \
        NSW       = 0          !number of steps for ionic relaxation loop\n \
        IBRION    = -1         !ionic relax: 0-MD 1-quasi-New 2-CG\n \
        #POTIM     = 0.2\n \
        #NFREE     = 3          ! initial steepest desc\n \
        ####################################################################\n \
        #LDA+U parameters\n# \
        LDAU     = .TRUE.\n \
        #LDAUTYPE = 2\n \
        #LDAUL    = -1    2    -1 \n \
        #LDAUU    =  0   2.50   0      / Na  Ru  O   \n \
        #LDAUJ    = 0.0  0.0 0.0 \n \
        #LDAUPRINT= 2\n \
        #LMAXMIX  = 4\n\n \
        #ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n \
        LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n \
        NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n \
        #NCORE   = 4\n'

        self.INCAR_wannier = ' System    = Wannier (VASP)\n####################################################################\n#Startparameter\n PREC      = Accurate  ! medium, high low\n ISTART    = 1         ! job   : 0-new  1-cont  2-samecut\n ICHARG    = 11        ! charge: 1-file 2-atom 10-const\n ISPIN     = 2         ! spin polarized calculation?\n#MAGMOM    = 4*4 20*0\n####################################################################\n#SOC parameter \nLSORBIT    =.TRUE.\nSAXIS     = 0 0 1 ! direction of the magnetic field \nMAGMOM    = 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 \nNBANDS = 40\nGGA_COMPAT=.FALSE.\nLORBMOM   =.TRUE.\n#Selects the HSE06 hybrid function\n#LHFCALC   = .TRUE. \n#HFSCREEN  = 0.2 \n#AEXX      = 0.25  \n#ALGO      = D \n#TIME      = 0.4\n#Electronic Relaxation\n ENCUT     = 450.0     ! cut-off energy (eV)\n ISMEAR    = -5        ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n SIGMA     = 0.1       ! broadening (eV) \n NELM      = 500       ! maximum number of electronic SC steps\n#NELMIN    = 5         ! minimum number of electronic SC steps\n#NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n#IALGO     = 48        ! algorithm\n#LDIAG     = T         ! sub-space diagonalisation\n#LREAL     = Auto      ! real-space projection\n#VOSKOWN = 1           !\n#IDIPOL   = 3\n#DIPOL    = 0.0 0.0 0.5\n#Mixing Parameters\n#AMIX      = 0.05\n#BMIX      = 0.01\n#AMIX_MAG  = 0.05\n#BMIX_MAG  = 0.0001\n#MAXMIX    = 80\n#Writing Parameters\n LWAVE     = F          !write WAVECAR\n LCHARG    = F          !write CHGCAR\n LVTOT     = F          !write LOCPOT\n#RWIGS     = 0.77 1.17\n LORBIT    = 11\n  LWANNIER90=.TRUE.\n#LAECHG    = T\n#LELF      = TRUE \n#NBLOCK    = 1\n#KBLOCK    = 50\n KSPACING  = 0.2\n KGAMMA    = TRUE\n#LPARD     = T \n#IBAND     = 374 375 376 \n#KPUSE     = 62 71 78   \n#LSEPB     = T\n#LSEPK     = T\n####################################################################\n#Ionic relaxation\n ISIF      = 3          !stress and relaxation        \n ISYM      = -1         !symmetry\n EDIFFG    = -0.01      !stopping-criterion for ionic relaxation loop            \n NSW       = 0          !number of steps for ionic relaxation loop\n IBRION    = -1         !ionic relax: 0-MD 1-quasi-New 2-CG\n#POTIM     = 0.2\n#NFREE     = 3          ! initial steepest desc\n####################################################################\n#LDA+U parameters\n#LDAU     = .TRUE.\n#LDAUTYPE = 2\n#LDAUL    = -1    2    -1 \n#LDAUU    =  0   2.50   0      / Na  Ru  O   \n#LDAUJ    = 0.0  0.0 0.0 \n#LDAUPRINT= 2\n #LMAXMIX  = 4\n\n#ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n#NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n#NCORE   = 4\n'

        self.INCAR_scf_nonsoc = ' System    = Scf (VASP)\n####################################################################\n#Startparameter\n PREC      = Accurate  ! medium, high low\n ISTART    = 0         ! job   : 0-new  1-cont  2-samecut\n ICHARG    = 2         ! charge: 1-file 2-atom 10-const\n ISPIN     = 2         ! spin polarized calculation?\n#MAGMOM    = 4*4 20*0\n####################################################################\n#SOC parameter \n#LSORBIT    =.TRUE.\n#SAXIS     = 0 0 1 ! direction of the magnetic field \n#MAGMOM    = 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 \n#NBANDS = 40\n#GGA_COMPAT=.FALSE.\n#LORBMOM   =.TRUE.\n#Selects the HSE06 hybrid function\n#LHFCALC   = .TRUE. \n#HFSCREEN  = 0.2 \n#AEXX      = 0.25  \n#ALGO      = D \n#TIME      = 0.4\n#Electronic Relaxation\n ENCUT     = 450.0     ! cut-off energy (eV)\n ISMEAR    = -5        ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n SIGMA     = 0.1       ! broadening (eV) \n NELM      = 1000       ! maximum number of electronic SC steps\n#NELMIN    = 5         ! minimum number of electronic SC steps\n#NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n#IALGO     = 48        ! algorithm\n#LDIAG     = T         ! sub-space diagonalisation\n#LREAL     = Auto      ! real-space projection\n#VOSKOWN = 1           !\n#IDIPOL   = 3\n#DIPOL    = 0.0 0.0 0.5\n#Mixing Parameters\n#AMIX      = 0.05\n#BMIX      = 0.01\n#AMIX_MAG  = 0.05\n#BMIX_MAG  = 0.0001\n#MAXMIX    = 80\n#Writing Parameters\n LWAVE     = T          !write WAVECAR\n LCHARG    = T          !write CHGCAR\n LVTOT     = F          !write LOCPOT\n#RWIGS     = 0.77 1.17\n LORBIT    = 11\n#LAECHG    = T\n#LELF      = TRUE \n#NBLOCK    = 1\n#KBLOCK    = 50\n KSPACING  = 0.2\n KGAMMA    = TRUE\n#LPARD     = T \n#IBAND     = 374 375 376 \n#KPUSE     = 62 71 78   \n#LSEPB     = T\n#LSEPK     = T\n####################################################################\n#Ionic relaxation\n ISIF      = 3          !stress and relaxation        \n ISYM      = -1         !symmetry\n EDIFFG    = -0.01      !stopping-criterion for ionic relaxation loop            \n NSW       = 0          !number of steps for ionic relaxation loop\n IBRION    = -1         !ionic relax: 0-MD 1-quasi-New 2-CG\n#POTIM     = 0.2\n#NFREE     = 3          ! initial steepest desc\n####################################################################\n#LDA+U parameters\n#LDAU     = .TRUE.\n#LDAUTYPE = 2\n#LDAUL    = -1    2    -1 \n#LDAUU    =  0   2.50   0      / Na  Ru  O   \n#LDAUJ    = 0.0  0.0 0.0 \n#LDAUPRINT= 2\n #LMAXMIX  = 4\n\n#ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n#NCORE   = 4\n'

        self.INCAR_band_nonsoc = ' System    = Band (VASP)\n####################################################################\n#Startparameter\n PREC      = Accurate  ! medium, high low\n ISTART    = 1         ! job   : 0-new  1-cont  2-samecut\n ICHARG    = 11        ! charge: 1-file 2-atom 10-const\n ISPIN     = 2         ! spin polarized calculation?\n#MAGMOM    = 4*4 20*0\n####################################################################\n#SOC parameter \n#LSORBIT    =.TRUE.\n#SAXIS     = 0 0 1 ! direction of the magnetic field \n#MAGMOM    = 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 \n#NBANDS = 40\n#GGA_COMPAT=.FALSE.\n#LORBMOM   =.TRUE.\n#Selects the HSE06 hybrid function\n#LHFCALC   = .TRUE. \n#HFSCREEN  = 0.2 \n#AEXX      = 0.25  \n#ALGO      = D \n#TIME      = 0.4\n#Electronic Relaxation\n ENCUT     = 450.0     ! cut-off energy (eV)\n ISMEAR    = 0         ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n SIGMA     = 0.1       ! broadening (eV) \n NELM      = 500       ! maximum number of electronic SC steps\n#NELMIN    = 5         ! minimum number of electronic SC steps\n#NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n#IALGO     = 48        ! algorithm\n#LDIAG     = T         ! sub-space diagonalisation\n#LREAL     = Auto      ! real-space projection\n#VOSKOWN = 1           !\n#IDIPOL   = 3\n#DIPOL    = 0.0 0.0 0.5\n#Mixing Parameters\n#AMIX      = 0.05\n#BMIX      = 0.01\n#AMIX_MAG  = 0.05\n#BMIX_MAG  = 0.0001\n#MAXMIX    = 80\n#Writing Parameters\n LWAVE     = F          !write WAVECAR\n LCHARG    = F          !write CHGCAR\n LVTOT     = F          !write LOCPOT\n#RWIGS     = 0.77 1.17\n LORBIT    = 11\n#LAECHG    = T\n#LELF      = TRUE \n#NBLOCK    = 1\n#KBLOCK    = 50\n#KSPACING  = 0.3\n#KGAMMA    = TRUE\n#LPARD     = T \n#IBAND     = 374 375 376 \n#KPUSE     = 62 71 78   \n#LSEPB     = T\n#LSEPK     = T\n####################################################################\n#Ionic relaxation\n ISIF      = 3          !stress and relaxation        \n ISYM      = -1         !symmetry\n EDIFFG    = -0.01      !stopping-criterion for ionic relaxation loop            \n NSW       = 0          !number of steps for ionic relaxation loop\n IBRION    = -1         !ionic relax: 0-MD 1-quasi-New 2-CG\n#POTIM     = 0.2\n#NFREE     = 3          ! initial steepest desc\n####################################################################\n#LDA+U parameters\n#LDAU     = .TRUE.\n#LDAUTYPE = 2\n#LDAUL    = -1    2    -1 \n#LDAUU    =  0   2.50   0      / Na  Ru  O   \n#LDAUJ    = 0.0  0.0 0.0 \n#LDAUPRINT= 2\n #LMAXMIX  = 4\n\n#ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n#NCORE   = 4\n'

        self.INCAR_wannier_nonsoc = ' System    = Wannier (VASP)\n####################################################################\n#Startparameter\n PREC      = Accurate  ! medium, high low\n ISTART    = 1         ! job   : 0-new  1-cont  2-samecut\n ICHARG    = 11        ! charge: 1-file 2-atom 10-const\n ISPIN     = 2         ! spin polarized calculation?\n#MAGMOM    = 4*4 20*0\n####################################################################\n#SOC parameter \n#LSORBIT    =.TRUE.\n#SAXIS     = 0 0 1 ! direction of the magnetic field \n#MAGMOM    = 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 \nNBANDS = 40\n#GGA_COMPAT=.FALSE.\n#LORBMOM   =.TRUE.\n#Selects the HSE06 hybrid function\n#LHFCALC   = .TRUE. \n#HFSCREEN  = 0.2 \n#AEXX      = 0.25  \n#ALGO      = D \n#TIME      = 0.4\n#Electronic Relaxation\n ENCUT     = 450.0     ! cut-off energy (eV)\n ISMEAR    = -5        ! 0:Gaussian -1:Fermi -5:Tetrahedron method\n SIGMA     = 0.1       ! broadening (eV) \n NELM      = 500       ! maximum number of electronic SC steps\n#NELMIN    = 5         ! minimum number of electronic SC steps\n#NELMDL    = -8        ! number of non-selfconsistent steps at the beginning\n EDIFF     = 1E-6      ! stopping-criterion for electronic SC-loop\n#IALGO     = 48        ! algorithm\n#LDIAG     = T         ! sub-space diagonalisation\n#LREAL     = Auto      ! real-space projection\n#VOSKOWN = 1           !\n#IDIPOL   = 3\n#DIPOL    = 0.0 0.0 0.5\n#Mixing Parameters\n#AMIX      = 0.05\n#BMIX      = 0.01\n#AMIX_MAG  = 0.05\n#BMIX_MAG  = 0.0001\n#MAXMIX    = 80\n#Writing Parameters\n LWAVE     = F          !write WAVECAR\n LCHARG    = F          !write CHGCAR\n LVTOT     = F          !write LOCPOT\n#RWIGS     = 0.77 1.17\n LORBIT    = 11\n  LWANNIER90=.TRUE.\n#LAECHG    = T\n#LELF      = TRUE \n#NBLOCK    = 1\n#KBLOCK    = 50\n KSPACING  = 0.2\n KGAMMA    = TRUE\n#LPARD     = T \n#IBAND     = 374 375 376 \n#KPUSE     = 62 71 78   \n#LSEPB     = T\n#LSEPK     = T\n####################################################################\n#Ionic relaxation\n ISIF      = 3          !stress and relaxation        \n ISYM      = -1         !symmetry\n EDIFFG    = -0.01      !stopping-criterion for ionic relaxation loop            \n NSW       = 0          !number of steps for ionic relaxation loop\n IBRION    = -1         !ionic relax: 0-MD 1-quasi-New 2-CG\n#POTIM     = 0.2\n#NFREE     = 3          ! initial steepest desc\n####################################################################\n#LDA+U parameters\n#LDAU     = .TRUE.\n#LDAUTYPE = 2\n#LDAUL    = -1    2    -1 \n#LDAUU    =  0   2.50   0      / Na  Ru  O   \n#LDAUJ    = 0.0  0.0 0.0 \n#LDAUPRINT= 2\n #LMAXMIX  = 4\n\n#ALGO    = Fast      ! electronic minimization method. Fast combines speed and robustness; VeryFast and Normal can also be used\n LPLANE  = T         ! T, reduce the communication during FFT. check balance among nodes!\n#NPAR    = 4         ! parallel over bands, equals an integer close to sqrt(total-cores-you-use) to speed up the calculation; the total cores must be divided by this number\n#NCORE   = 4\n'