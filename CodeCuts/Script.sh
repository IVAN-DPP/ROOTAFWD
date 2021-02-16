#*******************************************#
# Author: Ivan Dario Piernagorda Peña     
# Date:   15/02/2020                  
# Title:  Photon Polarization Tables                              
#*******************************************#


#!/bin/sh

echo "---------- OPEN ----------"

rm ./TABLES/ListTables/*

# List of Tables X Energy of Photon and Electron Beam

# PARA

ls -d1 TABLES/Tables/*1.3_PARA* >  TABLES/ListTables/Beam_PARA_1.3_4199.dat
ls -d1 TABLES/Tables/*1.5_PARA* |grep -i "4072" > TABLES/ListTables/Beam_PARA_1.5_4072.dat 
ls -d1 TABLES/Tables/*1.5_PARA* |grep -i "4482" > TABLES/ListTables/Beam_PARA_1.5_4482.dat
ls -d1 TABLES/Tables/*1.7_PARA* |grep -i "4072" > TABLES/ListTables/Beam_PARA_1.7_4072.dat 
ls -d1 TABLES/Tables/*1.7_PARA* |grep -i "4726" > TABLES/ListTables/Beam_PARA_1.7_4726.dat
ls -d1 TABLES/Tables/*1.7_PARA* |grep -i "4756" > TABLES/ListTables/Beam_PARA_1.7_4756.dat
ls -d1 TABLES/Tables/*1.9_PARA* > TABLES/ListTables/Beam_PARA_1.9_5052.dat
ls -d1 TABLES/Tables/*2.1_PARA* |grep -i "5052" > TABLES/ListTables/Beam_PARA_2.1_5052.dat 
ls -d1 TABLES/Tables/*2.1_PARA* |grep -i "5166" > TABLES/ListTables/Beam_PARA_2.1_5166.dat 
ls -d1 TABLES/Tables/*2.3_PARA* |grep -i "5166" > TABLES/ListTables/Beam_PARA_2.3_5166.dat

# PERP 

ls -d1 TABLES/Tables/*1.3_PERP* >  TABLES/ListTables/Beam_PERP_1.3_4199.dat
ls -d1 TABLES/Tables/*1.5_PERP* |grep -i "4072" > TABLES/ListTables/Beam_PERP_1.5_4072.dat 
ls -d1 TABLES/Tables/*1.5_PERP* |grep -i "4482" > TABLES/ListTables/Beam_PERP_1.5_4482.dat
ls -d1 TABLES/Tables/*1.7_PERP* |grep -i "4072" > TABLES/ListTables/Beam_PERP_1.7_4072.dat 
ls -d1 TABLES/Tables/*1.7_PERP* |grep -i "4726" > TABLES/ListTables/Beam_PERP_1.7_4726.dat
ls -d1 TABLES/Tables/*1.7_PERP* |grep -i "4756" > TABLES/ListTables/Beam_PERP_1.7_4756.dat
ls -d1 TABLES/Tables/*1.9_PERP* > TABLES/ListTables/Beam_PERP_1.9_5052.dat
ls -d1 TABLES/Tables/*2.1_PERP* |grep -i "5052" > TABLES/ListTables/Beam_PERP_2.1_5052.dat 
ls -d1 TABLES/Tables/*2.1_PERP* |grep -i "5166" > TABLES/ListTables/Beam_PERP_2.1_5166.dat 
ls -d1 TABLES/Tables/*2.3_PERP* |grep -i "5166" > TABLES/ListTables/Beam_PERP_2.3_5166.dat


echo "TABLES/ListTables Files Save"

#ListFiles .root

rm -rf ./SKIMS/ListFiles.txt

ls -1d ./SKIMS/PARA/*.root > ./SKIMS/ListFiles.txt
ls -1d ./SKIMS/PERP/*.root >> ./SKIMS/ListFiles.txt 

echo "./SKIMS/ListFiles.txt Save"

chmod +x Script.sh 

echo "---------- EXIT ----------"
