# ############################################
# Function to generate a random number from a uniform distribution
proc randUniform {min max} {
    return [expr {$min + ($max - $min) * rand()}]
}
# ############################################
# Function to generate a random number from a normal distribution
proc myRandNormal {mean stdev} {
    set u1 [expr {rand()}]
    set u2 [expr {rand()}]
    set z [expr {sqrt(-2 * log($u1)) * cos(2 * $u2 * 3.14159265358979323846)}]
    return [expr {$mean + $stdev * $z}]
}
# ############################################
# Function to generate a random number from a triangular distribution
proc randTriangular {a b c} {
    set u [expr {rand()}]
    if {$u <= ($c - $a) / ($b - $a)} {
        return [expr {$a + sqrt($u * ($b - $a) * ($c - $a))}]
    } else {
        return [expr {$b - sqrt((1 - $u) * ($b - $a) * ($b - $c))}]
    }
}
# ############################################
# Function to find Maximum Data
proc findMax {arr} {
    set maxVal [lindex $arr 0] ;# Initialize with the first element
    foreach num $arr {
        if {$num > $maxVal} {
            set maxVal $num
        }
    }
    return $maxVal
}
# ############################################
# Function to find Minimum Data
proc findMin {arr} {
    set maxVal [lindex $arr 0] ;# Initialize with the first element
    foreach num $arr {
        if {$num < $maxVal} {
            set minVal $num
        }
    }
    return $minVal
}