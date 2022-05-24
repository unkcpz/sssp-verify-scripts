#!/usr/bin/env bash

## AFTER CREATE A NEW PROFILE< RUN THIS SCRIPYT TO CONFIGURE COMPUTER AND CODES

profile="$1"
WDIR="$HOME/Projects/WP-SSSP/sssp-verify-scripts/computer_and_code"

# localhost
verdi -p ${profile} computer setup --config $WDIR/localhost/theos38/computer-setup.yaml -n
verdi -p ${profile} computer configure local localhost --config $WDIR/localhost/theos38/computer-configure.yaml -n
verdi -p ${profile} code setup --config $WDIR/localhost/codes/pw-6.7.yaml
verdi -p ${profile} code setup --config $WDIR/localhost/codes/ph-6.7.yaml

# eiger mr0
verdi -p ${profile} computer setup --config $WDIR/eiger.cscs.ch/multicore-mr0/computer-setup.yaml -n
verdi -p ${profile} computer configure ssh eiger-mc-mr0 --config $WDIR/eiger.cscs.ch/multicore-mr0/computer-configure.yaml
verdi -p ${profile} code setup --config $WDIR/eiger.cscs.ch/codes/pw-7.0-multicore.yaml --computer eiger-mc-mr0
verdi -p ${profile} code setup --config $WDIR/eiger.cscs.ch/codes/ph-7.0-multicore.yaml --computer eiger-mc-mr0

# eiger mr32
verdi -p ${profile} computer setup --config $WDIR/eiger.cscs.ch/multicore-mr32/computer-setup.yaml -n
verdi -p ${profile} computer configure ssh eiger-mc-mr32 --config $WDIR/eiger.cscs.ch/multicore-mr32/computer-configure.yaml
verdi -p ${profile} code setup --config $WDIR/eiger.cscs.ch/codes/pw-7.0-multicore.yaml --computer eiger-mc-mr32
verdi -p ${profile} code setup --config $WDIR/eiger.cscs.ch/codes/ph-7.0-multicore.yaml --computer eiger-mc-mr32

# eiger mr32
verdi -p ${profile} computer setup --config $WDIR/imxgesrv1.epfl.ch/imxgesrv1/computer-setup.yaml -n
verdi -p ${profile} computer configure ssh imxgesrv1 --config $WDIR/imxgesrv1.epfl.ch/imxgesrv1/computer-configure.yaml
verdi -p ${profile} code setup --config $WDIR/imxgesrv1.epfl.ch/codes/pw-7.0.yaml --computer imxgesrv1
verdi -p ${profile} code setup --config $WDIR/imxgesrv1.epfl.ch/codes/ph-7.0.yaml --computer imxgesrv1


# 
verdi -p ${profile} config set caching.enabled_for -a aiida.calculations:quantumespresso.pw
verdi -p ${profile} config set caching.enabled_for -a aiida.calculations:quantumespresso.ph