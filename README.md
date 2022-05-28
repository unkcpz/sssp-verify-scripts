# sssp verify scripts

This repository store the pseudopotential files as input for SSSP verification and the scripts to run the verification.

## Command to run

Batch launching precheck or standard verification of pseudos in a folder:
```bash
╰$ for FILE in _sssp_pbe/<element>/*.upf; do
echo $FILE && python run_verify.py --mode <precheck/standard> --computer mr32 -- $FILE
done
```

## Computer to run

To make things easy and clean, to avoid caching issue that bands workchain of standard that may use precheck calcjob which may be cleaned and fail. 
I'll try to sure the precheck and standard are run on different machine/account.
The principles are 

- precheck verification first go to imxgesrv1
- if precheck on eiger (mr0 if not mentioned), then standard should use another account. mr0/mr32 or mr32/mr0

## Pseudopotential libraries 

### Lanthanides

Lanthanides pseudopotentials are triky since the lanthanides contain f state which is hard to precisely describe in DFT. 

### ATOMPAW 

The latest version of ATOMPAW is 4.2.0.0 which ... SCAN

The La and Ce run on Wentzcovich sinec these two can not generated with latest version.

### ONCVPSP

SG15 inputs, some (list them) may failed to generate

### PSLibrary

The PSLibrary pseudopotentials are gerenated by using `ld1.x`. 
I use inputs file from [PSLibrary](https://github.com/dalcorso/pslibrary) and re-generate pseudopotentials using docker image `pspgen/ld1:0.1.0`. 
`ld1.x` in the image is compiled from Quantum-Espresso v6.3-MAX with `gfortran7` and `lapack` the in package.

The libraries re-generated are `PBE` and `PBEsol` with inputs:

1. PAW/HIGH from https://github.com/dalcorso/pslibrary/blob/master/paw_ps_high.job
2. US/HIGH from https://github.com/dalcorso/pslibrary/blob/master/us_ps_high.job
3. PAW/LOW from
4. US/LOW from 
5. PAW/0.x from https://github.com/dalcorso/pslibrary/blob/master/paw_ps_collection.job 
6. US/0.x from https://github.com/dalcorso/pslibrary/blob/master/us_ps_collection.job

The 0.x PPs have been tested in some cases and no error has been reported so far. The PPs that were in pslibrary.0.1 (See [PSLibrary ChangeLog file](https://github.com/dalcorso/pslibrary/blob/master/ChangeLog)) have been tested extensively

#### Notes 

- Cs.paw.z_9.ld1.psl.v1.0.0-low.upf: manually set z=9 and modified z since -5 in UPF.
- Cs.us.z_9.ld1.psl.v1.0.0-low.upf: manually set z=9 and modified z since -5 in UPF.

## Name convention of pseudopotential UPF file

The filename (it is easy to be done by using `rename` command of linux and script `psp_fn_z_set.ipynb`) will be used to deduct the label name when submit the verification.
The naming conventios is: `<element>.<psp_type>.z_<enum>.<tool>.<library>.<version>.<extra>.upf`

    element: the element of pseudopotential.
    psp_type: nc for norm-conserving, us for ultra-solf and paw for paw pseudopotentials
    enum: is the number of valence electrons
    tool: is the tool used for pseudopotential generation
        oncvpsp4: ONCVPSP code version 4.0.1
        oncvpsp3: this is where the dojo from version 3.3.0
        ld1: used to generate PSlibrary
        uspp: Vanderbilt code used to generate GBRV library
        atompaw: ATOMPAW code for jth and wentzcovitch
    library: is the pseudopotential library, we will run on: 
        dojo for PSEUDO-DOJO norm-conserving pseudos
        jth for ATOMPAW paw pseudos from pseudo-dojo webpage
        sg15 for SG15 pseudopotestials re-generated with oncvpsp 4.0.1
        psl for PSLibrary pseudos
        gbrv for ultrasoft pseudos of GBRV library by Vanderbilt
        wentzcovitch for lanthanides by wentzcovith, comment: legacy is downloaded, new is generated with same inputs but on ATOMPAW 4.2.0.0
    version: is for version of pseudopotential. `v0` for no version specified.
    extra: in the extra information related to pseudos that can be used to distingish the small distinctien. 
         
# Acknowledgements

The SSSP is a verification effort, but it is very important to give credit to the different authors that have generated the pseudopotential libraries that are tested here, and to the original methodological developments that underlie the generation of these pseudopotential tables and datasets. Citations can e.g. be taken from this list (contact us if we need to add more), appropriately for the libraries, methods, and datasets used. Please make an effort to acknowledge these and to ensure reproducibility of your calculations by listing/citing all pseudopotentials used.

## VERIFICATION
- SSSP: G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet and N. Marzari, npj Computational Materials 4, 72 (2018).
WEB: http://materialscloud.org/sssp.
- K. Lejaeghere et al., Science 351 (6280), 1415 (2016).
DOI: 10.1126/science.aad3000, WEB: http://molmod.ugent.be/deltacodesdft. An open-access copy is available from Cottenier's page.
- https://github.com/aiidateam/aiida-sssp-workflow

## LIBRARIES
- GBRV: K. F. Garrity, J. W. Bennett, K. M. Rabe, and D. Vanderbilt, Comput. Mater. Sci. 81, 446 (2014).
DOI: 10.1016/j.commatsci.2013.08.053, WEB: http://www.physics.rutgers.edu/gbrv, LICENSE: GNU General Public License (version 3).
- SG15: M. Schlipf and F. Gygi, Comp. Phys. Comm. 196, 36 (2015).
DOI: 10.1016/j.cpc.2015.05.011, WEB: http://www.quantum-simulation.org/potentials/sg15_oncv, LICENSE: Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0).
- Goedecker: A. Willand, Y. O. Kvashnin, L. Genovese, A. Vázquez-Mayagoitia, A. K. Deb, A. Sadeghi, T. Deutsch, and S. Goedecker, J. Chem. Phys. 138, 104109 (2013).
DOI: 10.1063/1.4793260, WEB: http://bigdft.org/Wiki/index.php?title=New_Soft-Accurate_NLCC_pseudopotentials, LICENSE: Creative Commons Attribution 3.0 Unported License (CC BY 3.0).
- Pslibrary 0.3.1: E. Küçükbenli et al., arXiv:1404.3015.
WEB: http://theossrv1.epfl.ch/Main/Pseudopotentials, LICENSE: GNU General Public License (version 2 or later).
- Pslibrary 1.0.0: A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014).
DOI: 10.1016/j.commatsci.2014.07.043, WEB: http://www.quantum-espresso.org/pseudopotentials, LICENSE: GNU General Public License (version 2 or later).
- RE Wentzcovitch: M. Topsakal and R. M. Wentzcovitch, Comput. Mater. Sci. 95, 263 (2014).
DOI: 10.1016/j.commatsci.2014.07.030, WEB: http://www.mineralscloud.com/resources/repaw/index.shtml.
- Pseudo Dojo: M.J. van Setten, M. Giantomassi, E. Bousquet, M.J. Verstraete, D.R. Hamann, X. Gonze, G.-M. Rignanese, Comp. Phys. Comm. 226, 39 (2018).
DOI: 10.1016/j.cpc.2018.01.012, WEB: http://www.pseudo-dojo.org/.

## METHODS
- Ultrasoft pseudopotentials: D. Vanderbilt, Phys. Rev. B 41, 7892(R) (1990).
DOI: 10.1103/PhysRevB.41.7892, WEB: http://physics.rutgers.edu/~dhv/uspp.
- Projector-augmented wave (PAW) method: P. E. Blöchl, Phys. Rev. B 50, 17953 (1994).
DOI: 10.1103/PhysRevB.50.17953.
- Norm conserving, multiple projectors pseudopotentials: D. R. Hamann, Phys. Rev. B 88, 085117 (2013).
DOI: 10.1103/PhysRevB.88.085117, WEB: http://www.mat-simresearch.com.
- Separable dual-space Gaussian pseudopotentials: S. Goedecker, M. Teter, J. Hutter, Phys. Rev. B 54, 1703 (1996).
DOI: 10.1103/PhysRevB.54.1703, WEB: http://cp2k.web.psi.ch/potentials.

## OTHER SITES
- http://www.pseudo-dojo.org
- http://opium.sourceforge.net