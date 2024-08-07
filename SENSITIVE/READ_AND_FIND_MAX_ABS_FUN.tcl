proc find_max_abs_value {filename} {
    # Open the input file for reading
    set inputFile [open $filename r]

    # Initialize the maximum absolute value
    set maxAbsValue 0.0000000

    # Read each line from the file
    while {[gets $inputFile line] != -1} {
        # Split the line by space or semicolon (adjust as needed)
        set columns [split $line " "]

        # Extract the value from the second column (index 1)
        set value [lindex $columns 1]

        # Calculate the absolute value
        set absValue [expr {abs($value)}]
		puts "absValue: $absValue"

        # Update the maximum absolute value if needed
        if {$absValue > $maxAbsValue} {
            set maxAbsValue $absValue
        }
    }

    # Close the input file
    close $inputFile

    # Return the result
    return $maxAbsValue
}


