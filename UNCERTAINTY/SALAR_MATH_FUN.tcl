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
    set pi 3.141592653589793
    set z [expr {sqrt(-2 * log($u1)) * cos(2 * $u2 * pi)}]
    return [expr {$mean + $stdev * $z}]
}
##############################################
# Function to calculate the log-normal distribution
proc logNormalDist {x mu sigma} {
    # Ensure sigma is positive
    if {$sigma <= 0} {
        error "Sigma must be positive for a valid log-normal distribution"
    }

    # Constants
    set pi 3.141592653589793
    set sqrt_2pi [expr sqrt(2.0 * $pi)]

    # Calculate the log-normal PDF
    set z [expr (log($x) - $mu) / $sigma]
    set pdf [expr (1.0 / ($x * $sqrt_2pi * $sigma)) * exp(-0.5 * $z * $z)]
    return $pdf
}

# Function to generate a random sample from the log-normal distribution
proc logNormalRandom {mu sigma} {
    # Generate a random number from a standard normal distribution
    set z [expr rand() * 2.0 - 1.0] ;# Adjusting rand() output
    # Transform to log-normal
    set sample [expr exp($mu + $sigma * $z)]
    return $sample
}

# Example usage
#set mu 1.0   ;# Mean of the underlying normal distribution
#set sigma 0.5 ;# Standard deviation of the underlying normal distribution

# Calculate PDF at a given point
#set x 2.0
#set pdf [logNormalDist $x $mu $sigma]
#puts "PDF at x=$x: $pdf"

# Generate a random sample
#set random_sample [logNormalRandom $mu $sigma]
#puts "Random sample from Log-Normal Distribution: $random_sample"

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
