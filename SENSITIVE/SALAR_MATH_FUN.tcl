# ############################################

# Function to generate a random number from a uniform distribution
proc randUniform {min max} {
    return [expr {$min + ($max - $min) * rand()}]
}
# Example usage:
#set minValue 10
#set maxValue 50
#set randomValue [randUniform $minValue $maxValue]
#puts "Random value from uniform distribution: $randomValue"

# ############################################

# Function to generate a random number from a normal distribution
proc myRandNormal {mean stdev} {
    set u1 [expr {rand()}]
    set u2 [expr {rand()}]
    set z [expr {sqrt(-2 * log($u1)) * cos(2 * $u2 * 3.14159265358979323846)}]
    return [expr {$mean + $stdev * $z}]
}

# Example usage:
#set mean 10
#set stdev 2
#set numSamples 100
#for {set i 0} {$i < $numSamples} {incr i} {
#    set randomValue [myRandNormal $mean $stdev]
#    puts "Sample $i: $randomValue"
#}

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
# Example usage:
#set minVal 10
#set maxVal 30
#set modeVal 25

#set randomValue [randTriangular $minVal $maxVal $modeVal]
#puts "Random value from triangular distribution: $randomValue"

# ############################################

# Function to generate a random number from a Gamma and Beta distribution
proc randGamma {shape scale} {
    set alpha $shape
    set beta $scale
    set d [expr {$alpha - 1.0/3.0}]
    set c [expr {1.0 / sqrt(9.0 * $d)}]
    while {1} {
        set u [expr {rand()}]
        set v [expr {pow([expr {rand()}], 3)}]
        set x [expr {$c * ($u - 0.5) / $v + 1}]
        set v [expr {$v * $x * $x * $x}]
        set u [expr {$u - 0.5}]
        if {[expr {$u < $c * $v}] || [expr {log($u + 1.0) + 0.5 * $v <= $d * (1 - $v + log($v))}]} {
            return [expr {$d * $x * $beta}]
        }
    }
}

# Function to generate a random number from a Beta distribution
proc randBeta {alpha beta} {
    set x [randGamma $alpha 1.0]
    set y [randGamma $beta 1.0]
    return [expr {$x / ($x + $y)}]
}

# Example usage:
#set alpha 2.0
#set beta 5.0
#set randomNumber [randBeta $alpha $beta]
#puts "Random Beta($alpha, $beta): $randomNumber"

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

# Example usage:
#set myArray {12 2 5 4 2 6 7 55 8 9 6 4}
#set minValue [findMin $myArray]
#puts "Minimum value in the array: $minValue"

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
# Example usage:
#set myArray {12 2 5 4 2 6 7 55 8 9 6 4}
#set maxValue [findMax $myArray]
#puts "Maximum value in the array: $maxValue"

# ############################################

# Function to generate an array of absolute values from an input array
proc absoluteArray {arr} {
    set absoluteArray {}
    foreach num $arr {
        lappend absoluteArray [expr {abs($num)}]
    }
    return $absoluteArray
}

# ############################################

# Function to find Quantile Data
proc quantile {cdf p} {
    # Input:
    #   cdf: Cumulative distribution function (a list of probabilities)
    #   p: Probability value (between 0 and 1)
    #
    # Output:
    #   Returns the quantile value corresponding to the given probability

    set n [llength $cdf]
    set idx [expr {int($p * ($n - 1))}]
    set lowerProb [lindex $cdf $idx]
    set upperProb [lindex $cdf [expr {$idx + 1}]]

    # Linear interpolation
    set quantileValue [expr {$lowerProb + ($p - $idx / ($n - 1)) * ($upperProb - $lowerProb)}]
    return $quantileValue
}

# Example usage:
#set myCDF {10.1 9.3 8.6 6.8 9.5} ;# Example CDF (sorted probabilities)
#set myProbability 0.5 ;# Desired probability (e.g., median)

#set quantileResult [quantile $myCDF $myProbability]
#puts "Quantile at $myProbability: $quantileResult"

# ############################################

# Function to calculate the mean (average) of an array
proc calculateMean {numbers} {
    set sum 0.0
    set count 0

    foreach num $numbers {
        incr sum $num
        incr count
    }

    if {$count > 0} {
        set mean [expr {$sum / $count}]
        return $mean
    } else {
        return "Array is empty"
    }
}

# Example usage:
#set myNumbers {10 20 30 40 50}
#set result [calculateMean $myNumbers]
#puts "Mean of the array: $result"

# ############################################

# Function to calculate the standard deviation of an array
proc standardDeviation {data} {
    set n [llength $data]
    if {$n <= 1} {
        return -1 ;# Not enough data to calculate standard deviation
    }
    
    set mean [expr {[tcl::mathop::+ {*}$data] / double($n)}]
    set sumSq 0.0
    foreach x $data {
        set sumSq [expr {$sumSq + ($x - $mean) ** 2}]
    }
    set variance [expr {$sumSq / ($n - 1)}]
    set stdDev [expr {sqrt($variance)}]
    return $stdDev
}

# Example usage:
#set numbers {4 8 6 5 3 2 8 9 7 5}
#set result [standardDeviation $numbers]
#puts "Standard deviation: $result"

# ############################################

# Function to calculate the standard deviation of an array
proc standardDeviation {data} {
    set n [llength $data]
    if {$n <= 1} {
        return -1 ;# Not enough data to calculate standard deviation
    }
    
    set mean [expr {[tcl::mathop::+ {*}$data] / double($n)}]
    set sumSq 0.0
    foreach x $data {
        set sumSq [expr {$sumSq + ($x - $mean) ** 2}]
    }
    set variance [expr {$sumSq / ($n - 1)}]
    set stdDev [expr {sqrt($variance)}]
    return $stdDev
}

# Example usage:
#set numbers {4 8 6 5 3 2 8 9 7 5}
#set result [standardDeviation $numbers]
#puts "Standard deviation: $result"

# ############################################

# Function to perform simple linear regression
proc linearRegression {xValues yValues} {
    set n [llength $xValues]
    set sumX 0
    set sumY 0
    set sumXY 0
    set sumX2 0

    # Calculate the necessary sums
    for {set i 0} {$i < $n} {incr i} {
        set x [lindex $xValues $i]
        set y [lindex $yValues $i]
        set sumX [expr {$sumX + $x}]
        set sumY [expr {$sumY + $y}]
        set sumXY [expr {$sumXY + $x * $y}]
        set sumX2 [expr {$sumX2 + $x * $x}]
    }

    # Calculate the slope (a) and y-intercept (b)
    set a [expr {($n * $sumXY - $sumX * $sumY) / ($n * $sumX2 - $sumX * $sumX)}]
    set b [expr {($sumY - $a * $sumX) / $n}]

    # Return the coefficients as a list
    return [list $a $b]
}

# Example usage:
#set xValues {1 2 3 4 5}
#set yValues {2 4 6 8 10}
#set coefficients [linearRegression $xValues $yValues]
#puts "Slope (a): [lindex $coefficients 0]"
#puts "Intercept (b): [lindex $coefficients 1]"

# ############################################

# Newton-Raphson method
proc newtonRaphson {x0 tolerance maxIter} {
    set x $x0
    for {set i 0} {$i < $maxIter} {incr i} {
        set f [expr {pow($x, 5) - 2 * pow($x, 3) + $x - 1000}]
        set df [expr {5 * pow($x, 4) - 6 * pow($x, 2) + 1}]
        if {$df == 0} {
            puts "Derivative is zero. No solution found."
            return
        }
        set dx [expr -$f / $df]
        set x [expr $x + $dx]
        if {abs($dx) < $tolerance} {
            puts "Root found: $x"
            return
        }
    }
    puts "Maximum iterations reached without convergence."
}

# Initial guess and parameters
#set x0 2.0
#set tolerance 1.0e-6
#set maxIter 100

# Call the Newton-Raphson method
#newtonRaphson $x0 $tolerance $maxIter

# ############################################

proc write_variable_to_file {filename variableName value} {
    set outfile [open $filename a+]

    # Write the variable name and value to the file
    #puts $outfile "$variableName = $value"
	puts $outfile "$variableName $value"

    # Close the file
    close $outfile
}

# Example usage:
#set filename "report1.out"
#set k 0
#set y 3
#write_variable_to_file $filename "k" $k

# ############################################
