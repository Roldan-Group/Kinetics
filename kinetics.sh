#!/bin/sh


        default="Matlab"

    if [ ! $2 ]; then
        software=$default
    else
        software=$2
    fi

    if [ $1 = "Maple" ]; then
        software=$1
    elif [ $1 = Mathematica ]; then
        software=$1
    elif [ $1 = Matlab ]; then
        software=$1
    else
        file=$1
      if [ ! $2 ]; then
          software=$default
      else 
          software=$2
      fi
   fi



#     UPDATED 07/2014
#  Increased efficientcy
# version Maple 3.X splits the maple file: 3.1 is for thermodynamics
# version Mathematica 0.X generates qtrans_3D in-situ as a function of the pressure



echo " $file $software"
   if [ $software = "Maple" ]; then
       ~/Software/KINETICS/kinetic2Maple_v3.1.pl $file
   elif  [ $software = "Mathematica" ]; then
       ~/Software/KINETICS/kinetic2Mathematica_v4.2.pl $file
   elif  [ $software = "Matlab" ]; then
       python3 ~/Software/KINETICS/Input2mk.py $file
       ~/Software/KINETICS/kinetic2Matlab_v4.2.pl $file.mk.in
       rm .m
   fi

