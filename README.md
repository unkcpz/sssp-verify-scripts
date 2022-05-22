# sssp verify scripts

This repository store the pseudopotential files as input for SSSP verification and the scripts to run the verification.

## Pseudopotential libraries 

### Lanthanides

Lanthanides pseudopotentials are triky since the lanthanides contain f state which is hard to precisely describe in DFT. 

### ATOMPAW 

The latest version of ATOMPAW is 4.2.0.0 which ... SCAN

The La and Ce run on Wentzcovich sinec these two can not generated with latest version.

### ONCVPSP

SG15 inputs, some (list them) may failed to generate

## Name convention of pseudopotential UPF file

The filename (it is easy to be done by using `rename` command of linux and script `psp_fn_z_set.py`) will be used to deduct the label name when submit the verification.
The naming conventios is: `<element>.<psp_type>.z_<enum>.<tool>.<library>.<version>.<extra>.upf`

    element: the element of pseudopotential.
    psp_type: nc for norm-conserving, us for ultra-solf and paw for paw pseudopotentials
    enum: is the number of valence electrons
    tool: is the tool used for pseudopotential generation
        oncvpsp4: ONCVPSP code version 4.0.1
        oncvpsp3: this is where the dojo from version 3.3.0
        ldx1: used to generate PSlibrary
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
         
