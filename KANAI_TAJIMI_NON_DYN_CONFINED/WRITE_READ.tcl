
# Write CSV file
# Step 1: Generate an array of 10 random numbers
set randomNumbers {}
for {set i 0} {$i < 10} {incr i} {
    lappend randomNumbers [expr {int(rand() * 100)}]
}

# Step 2: Write the array to a CSV file
set filename "random_numbers.csv"
set file [open $filename "w"]
puts $file [join $randomNumbers ","]
close $file
puts "Generated random numbers and wrote to $filename"

# Step 3: Open the CSV file and read the contents back into an array
set file [open $filename "r"]
set data [read $file]
close $file

# Read CSV file

# Convert CSV data back into an array
set readNumbers [split $data ","]
puts "Read numbers from $filename: $readNumbers"
