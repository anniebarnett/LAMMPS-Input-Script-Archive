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
variable        latparam equal 3.52
variable        comp     string Ni

#--- Create Lattice ---#
lattice        fcc ${latparam} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region         whole block -10 10 -10 10 -10 10 units lattice
create_box     2 whole
create_atoms   1 box

#set type 1 type/ratio 2 0.050 73932

mass 1 58.69   #Ni
mass 2 26.98 #Al

timestep .0005 

#--- Interactions ---#
pair_style eam/alloy
pair_coeff * * ./Ni-Al-2009.eam.alloy Ni Al

#--- Neighbor Cuttoff ---#
neighbor 2.0 bin
neigh_modify delay 10 every 1 check yes

#--- Variables 2 ---#
variable        Temp    equal 573 #300C
variable        N_init  equal count(all)
variable        time_tot equal "(dt*step)/1000"

#--- Computes ---#1 2 
compute csym all centro/atom fcc

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

compute disp all displace/atom

group type1_atoms type 1
group type2_atoms type 2

compute msd_t1 type1_atoms msd com no average no
compute msd_t2 type2_atoms msd com no average no

################################################################################
#---------------------------- Create Frenkel Pairs ----------------------------#
################################################################################
# 0 K relaxation
reset_timestep 0

#fix f1 all box/relax iso 0.0 vmax 0.001

thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms

dump d1 all custom 500 ${comp}-CG.*.dump id type element mass x y z c_csym c_potEnergy c_kinEnergy c_sigma[1] c_sigma[2] c_sigma[3] c_sigma[4] c_sigma[5] c_sigma[6]
dump_modify d1 element Ni Al
dump_modify d1 sort id

min_style cg
minimize 1e-25 1e-25 5000 10000

variable counter_interstitials equal 0
#velocity all create ${Temp} 7392611 mom yes rot yes dist gaussian

#unfix f1
undump d1
#--- Loop to generate Frenkel pairs ---#
variable        n loop 1
label           loop

variable fppl equal 5	## FPPL = FRENKEL PAIRS PER LOOP -> CHANGE THIS VALUE TO CHANGE THE 'DOSE RATE' OF YOUR SIMULATION ###
variable deposteps equal 10
variable rate  equal ${deposteps}/${fppl}
variable dumprate equal ${deposteps}/2
variable quenchsteps equal 100

#--- NPT relaxation ---#
fix npt all npt temp ${Temp} ${Temp} $(100*dt) aniso 0.0 0.0 $(1000*dt)

dump d2 all custom ${dumprate} ${comp}_npt${n}_${Temp}C.*.dump id type element mass x y z c_csym c_potEnergy c_kinEnergy c_sigma[1] c_sigma[2] c_sigma[3] c_sigma[4] c_sigma[5] c_sigma[6] c_disp[1] c_disp[2] c_disp[3] c_disp[4]
dump_modify d2 element Ni Al
dump_modify d2 sort id

#--- Box Settings ---#
region lattice_region_${n} block EDGE EDGE EDGE EDGE EDGE EDGE units box

#-- Assign Interstitials ---#   
#variable atype_unfiltered equal ceil(random(1,2,19299272))  # Randomly select atom type between 1 and 2 
variable atype_unfiltered equal 1
print "selected atom type = ${atype_unfiltered}"
variable atype equal ${atype_unfiltered}

#--- Create Vacancies ---#
variable num_delete_atoms equal ${fppl}
group vacancy_atoms type ${atype}
delete_atoms random count ${num_delete_atoms} yes vacancy_atoms NULL 263721

#--- Create Interstitials ---#
fix depo all deposit ${fppl} ${atype} ${rate} 261599 region lattice_region_${n} near 0.5 attempt 100 id max

thermo 50
thermo_style custom step temp atoms c_eatoms v_time_tot ke etotal press vol lx ly lz c_msd_t1[1] c_msd_t1[2] c_msd_t1[3] c_msd_t1[4] c_msd_t2[1] c_msd_t2[2] c_msd_t2[3] c_msd_t2[4]
thermo_modify lost ignore flush yes

run   ${deposteps}

dump d3 all custom 500 ${comp}-quench.*.dump id type element mass x y z c_csym c_potEnergy c_kinEnergy c_sigma[1] c_sigma[2] c_sigma[3] c_sigma[4] c_sigma[5] c_sigma[6]
dump_modify d3 element Ni Al
dump_modify d3 sort id

min_style cg
minimize 1e-25 1e-25 5000 10000

write_restart ${comp}_${Temp}C.${n}.restart

unfix depo
undump d2
unfix npt

next    n
jump    SELF loop
#------------------------------------------------------------------------------#

jump SELF END
label END
print       "All Done"


