for {set i 1} {$i <= $NUM_INCREMENT} {incr i} {# loop of mass increment
puts "Increment $i"
# FIND MAX. ABS. DISPLACEMENT
set filename02 "DATA/NON_DYN_CONFINED/DFree_$i.txt"
set result02 [find_max_abs_value $filename02]
puts "Maximum absolute displacement: $result02"

# write each max. abs. displacement in each iterations
write_variable_to_file DATA/NON_DYN_CONFINED/MAD_$i.txt "$i" $result02

# #################################################################
# FIND MAX. ABS. VELOCITY
set filename02 "DATA/NON_DYN_CONFINED/VFree_$i.txt"
set result02 [find_max_abs_value $filename02]
puts "Maximum absolute velocity: $result02"

# write each max. abs. velocity in each iterations
write_variable_to_file DATA/NON_DYN_CONFINED/MAV_$i.txt "$i" $result02

# #################################################################
# FIND MAX. ABS. ACCELERATION 
set filename02 "DATA/NON_DYN_CONFINED/AFree_$i.txt"
set result02 [find_max_abs_value $filename02]
puts "Maximum absolute acceleration: $result02"

# write each max. abs. acceleration in each iterations
write_variable_to_file DATA/NON_DYN_CONFINED/MAA_$i.txt "$i" $result02

# #################################################################
# FIND MAX. ABS. BASE SHEAR
set filename02 "DATA/NON_DYN_CONFINED/RBase_$i.txt"
set result02 [find_max_abs_value $filename02]
puts "Maximum absolute base shear: $result02"

# write each max. abs. acceleration in each iterations
write_variable_to_file DATA/NON_DYN_CONFINED/MAB_$i.txt "$i" $result02


}