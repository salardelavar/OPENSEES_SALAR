wipe all
file mkdir Truss; 				# create data directory
model BasicBuilder -ndm 2 -ndf 2
# create material
set Fy 4000.0;
set Fu 6000.0;
set E 200000.0;
set Sy [expr $Fy/$E];
set Su 0.1;
set Esh [expr ($Fu-$Fy)/($Su-$Sy)];
set b [expr ($Esh/$E)];
set NI 1000;#Analysis Number
set DISini 0.5;#Initial Displacement
set PushNode 3;#Incremental Displacement Node
set PushNodeDir 1;#Direction Incremental Displacement Node 1:x - 2:y
set x {0 3000.0 0 3000.0};
set y {0 0 3000.0 3000.0};
set a {400.0 500.0 600.0 700.0 800.0 400.0}; # Area of All Truss Element
for {set i 0} {$i<4} {incr i} {
node [expr $i+1] [lindex $x $i] [lindex $y $i]
}
uniaxialMaterial Steel01 1 $Fy $E $b
set conn {{2 1} {3 1} {4 1} {3 2} {4 2} {3 4}}; # Node Connection
for {set i 0} {$i < 6} {incr i} {
	set I [lindex $conn $i 0];
	set J [lindex $conn $i 1];
	set A [lindex $a $i];
	element truss [expr $i+1] $I $J $A 1
}
fix 1 1 1;
fix 2 0 1;
pattern Plain 1 "Linear" {
	load 3 1 0;
	load 3 1 -1;
	load 4 0 -1;
}
# create an analysis
constraints Plain
numberer RCM
test NormDispIncr 1.0e-12 6 0
algorithm Newton
system BandGen
integrator DisplacementControl $PushNode $PushNodeDir $DISini
analysis Static
# perform the analysis & print the results
# open output file
set OUT [open Truss/Truss_DamageIndex_OUTPUT.csv "w"];# write the data
	puts $OUT "Increment,DisX01,DisY01,DisX02,DisY02,DisX03,DisY03,DisX04,DisY04,STRAIN01,STRESS01,STRAIN02,STRESS02,STRAIN03,STRESS03,STRAIN04,STRESS04,STRAIN05,STRESS05,STRAIN06,STRESS06,ELEF01,ELEF02,ELEF03,ELEF04,ELEF05,ELEF06,DI01,DI02,DI03,DI04,DI05,DI06";
	puts "	==========================";
	puts "	          DIS X           ";
	puts "	==========================";
for {set i 0} {$i < $NI} {incr i 1} {
	analyze 1
	for {set k 1} {$k <= 4} {incr k 1} {
	set dx0$k [nodeDisp $k 1];
	set dy0$k [nodeDisp $k 2];
	}
	for {set j 1} {$j <= 6} {incr j 1} {
	set STRAIN0$j [eleResponse $j material strain];# Element Strain
	set STRESS0$j [eleResponse $j material stress];# Element Stress
	set Sd [eleResponse $j material strain];
	set X [expr (abs($Sd)-$Sy)/($Su-$Sy)];# Element Damage Index
	if { $X < 0 } { set X 0};# Lower Bound of Damage index
	if { $X > 1 } { set X 1};# Upper bound of Damage index
	set DI0$j $X;
	if { $X == 1 } {puts "Element $j Damage Index Reached to Ultimate Capacity"}
	#set ELEF01 [expr $STRESS$j*[lindex $a [expr $j-1]]];# Element Axial Force01
	}
	set ELEF01 [eleResponse 1 forces];# or set ELEF01 [expr $STRESS01*[lindex $a 0]];# Element Axial Force01
	set ELEF02 [eleResponse 2 forces];# or set ELEF02 [expr $STRESS02*[lindex $a 1]];# Element Axial Force02
	set ELEF03 [eleResponse 3 forces];# or set ELEF03 [expr $STRESS03*[lindex $a 2]];# Element Axial Force03
	set ELEF04 [eleResponse 4 forces];# or set ELEF04 [expr $STRESS04*[lindex $a 3]];# Element Axial Force04
	set ELEF05 [eleResponse 5 forces];# or set ELEF05 [expr $STRESS05*[lindex $a 4]];# Element Axial Force05
	set ELEF06 [eleResponse 6 forces];# or set ELEF06 [expr $STRESS06*[lindex $a 5]];# Element Axial Force06
	puts $OUT "[expr $i+1],$dx01,$dy01,$dx02,$dy02,$dx03,$dy03,$dx04,$dy04,$STRAIN01,$STRESS01,$STRAIN02,$STRESS02,$STRAIN03,$STRESS03,$STRAIN04,$STRESS04,$STRAIN05,$STRESS05,$STRAIN06,$STRESS06,$ELEF01,$ELEF02,$ELEF03,$ELEF04,$ELEF05,$ELEF06,$DI01,$DI02,$DI03,$DI04,$DI05,$DI06";
	puts "$dx03";
}
close $OUT;#close the file 