# Continuous Frenkel Pair Creation in MoNbTaW 
#
################################################################################
#---------------------------- Simulation Setup --------------------------------#
################################################################################
#--- Initialization ---#
units           metal
dimension       3
boundary        p p p
atom_style      atomic

#--- Variables 1 ---#
variable        latparam equal 3.23
variable        comp     string MoNbTaW

#--- Create Lattice ---#
lattice        bcc ${latparam} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region         whole block -20 20 -20 20 -20 20 units lattice
create_box     5 whole
create_atoms   1 box

set type 1 type/ratio 2 0.75 73932
set type 2 type/ratio 3 0.6667 478140
set type 3 type/ratio 5 0.5 937201

mass 1 95.95    #Mo
mass 2 92.90637 #Nb
mass 3 180.94788 #Ta
mass 4 50.9415 #V
mass 5 183.84   #W

reset_timestep 0
timestep .0005

#--- Interactions ---#
pair_style      hybrid/overlay tabgap eam/fs
pair_coeff      * * tabgap /home/abarne65/scratch16-mfalk1/annie/defects/Mo-Nb-Ta-V-W.tabgap Mo Nb Ta V W no yes
pair_coeff      * * eam/fs /home/abarne65/scratch16-mfalk1/annie/defects/Mo-Nb-Ta-V-W.eam.fs Mo Nb Ta V W

#--- Neighbor Cuttoff ---#
neighbor 2.0 bin
neigh_modify delay 10 every 1 check yes

#--- Variables 2 ---#
variable        Temp    equal 573 #300C
variable        N_init  equal count(all)
variable        time_tot equal "(dt*step)/1000"

#--- Computes ---#1 2 
compute csym all centro/atom bcc

compute potEnergy all pe/atom
compute eatoms all reduce sum c_potEnergy
compute kinEnergy all ke/atom

compute sigma all stress/atom NULL virial       #per atom stress tensor, 6 components below
compute aveSigma1 all reduce sum c_sigma[1]     #xx
compute aveSigma2 all reduce sum c_sigma[2]     #yy
compute aveSigma3 all reduce sum c_sigma[3]     #zz
compute aveSigma4 all reduce sum c_sigma[4]     #xy
compute aveSigma5 all reduce sum c_sigma[5]     #xz
compute aveSigma6 all reduce sum c_sigma[6]     #yz

compute flux all heat/flux kinEnergy potEnergy sigma
variable Jx equal c_flux[1]/vol
variable Jy equal c_flux[2]/vol
variable Jz equal c_flux[3]/vol

compute disp all displace/atom

group type1_atoms type 1
group type2_atoms type 2
group type3_atoms type 3
group type5_atoms type 5

compute msd_t1 type1_atoms msd com no average no
compute msd_t2 type2_atoms msd com no average no
compute msd_t3 type3_atoms msd com no average no
compute msd_t5 type5_atoms msd com no average no

################################################################################
#---------------------------- Create Frenkel Pairs ----------------------------#
################################################################################
# 0 K relaxation
min_style cg
minimize 1e-12 1e-8 10000 100000

variable counter_interstitials equal 0
velocity all create ${Temp} 7392611 mom yes rot yes dist gaussian

write_dump all custom initial_${comp}_0005.dump id type element mass x y z c_csym c_potEnergy c_kinEnergy c_sigma[1] c_sigma[2] c_sigma[3] c_sigma[4] c_sigma[5] c_sigma[6] c_disp[1] c_disp[2] c_disp[3] c_disp[4]

#--- Loop to generate Frenkel pairs ---#
variable        n loop 100
label           loop

variable fppl equal 1000
variable deposteps equal 100000
variable rate  equal ${deposteps}/${fppl}
variable dumprate equal ${deposteps}/50
variable quenchsteps equal 1000

#--- NPT relaxation ---#
fix npt all npt temp ${Temp} ${Temp} $(100*dt) aniso 0.0 0.0 $(1000*dt)
fix JJ all ave/correlate 100 100 10000 c_flux[1] c_flux[2] c_flux[3] type auto file 4-HT_5dpa_J0J.dat ave running

dump d1 all custom ${dumprate} ${comp}_npt${n}_${Temp}C.*.dump id type element mass x y z c_csym c_potEnergy c_kinEnergy c_sigma[1] c_sigma[2] c_sigma[3] c_sigma[4] c_sigma[5] c_sigma[6] c_disp[1] c_disp[2] c_disp[3] c_disp[4]
dump_modify d1 element Mo Nb Ta V W
dump_modify d1 sort id

#--- Box Settings ---#
region lattice_region_${n} block EDGE EDGE EDGE EDGE EDGE EDGE units box

#-- Assign Interstitials ---#   
variable atype_unfiltered equal ceil(random(1,4,19299272))  # Randomly select atom type between 1 and 4, replacing 4(V) with 5(W) 
print "selected atom type = ${atype_unfiltered}"
variable atype equal ${atype_unfiltered}
if "${atype} == 4" then "variable atype equal 5"
print "deposited atom type = ${atype}"

#--- Create Vacancies ---#
variable num_delete_atoms equal ${fppl}
group vacancy_atoms type ${atype}
delete_atoms random count ${num_delete_atoms} yes vacancy_atoms NULL 263721

#--- Create Interstitials ---#
fix depo all deposit ${fppl} ${atype} ${rate} 261599 region lattice_region_${n} near 0.5 attempt 100 id max

thermo 50
thermo_style custom step temp atoms c_eatoms v_time_tot ke etotal press vol lx ly lz c_msd_t1[1] c_msd_t1[2] c_msd_t1[3] c_msd_t1[4] c_msd_t2[1] c_msd_t2[2] c_msd_t2[3] c_msd_t2[4] c_msd_t3[1] c_msd_t3[2] c_msd_t3[3] c_msd_t3[4] c_msd_t5[1] c_msd_t5[2] c_msd_t5[3] c_msd_t5[4]
thermo_modify lost ignore flush yes

run   ${deposteps}

#write_dump all xyz 4_${Temp}C.${n}.xyz modify sort id element Mo Nb Ta V W
#write_data 4_${Temp}C.${n}.data
write_restart ${comp}_${Temp}C.${n}.restart

unfix depo
undump d1
unfix npt
unfix JJ

next    n
jump    SELF loop
#------------------------------------------------------------------------------#

jump SELF END
label END
print       "All Done"


