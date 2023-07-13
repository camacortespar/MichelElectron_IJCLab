#!/bin/bash

lar -c ../fcl/protodunevd_gen.fcl                                                                           -o ../output/protodunevd_gen_electron.root
lar -c ../fcl/protodunevd_refactored_g4_stage1.fcl ../output/protodunevd_gen_electron.root                  -o ../output/protodunevd_refactored_g4_stage1_electron.root  
lar -c ../fcl/protodunevd_refactored_g4_stage2.fcl ../output/protodunevd_refactored_g4_stage1_electron.root -o ../output/protodunevd_refactored_g4_stage2_electron.root   
lar -c ../fcl/protodunevd_detsim.fcl               ../output/protodunevd_refactored_g4_stage2_electron.root -o ../output/protodunevd_detsim_electron.root                
lar -c ../fcl/protodunevd_reco.fcl                 ../output/protodunevd_detsim_electron.root               -o ../output/protodunevd_reco_electron.root                   
