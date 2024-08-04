
#puts "       >>   IN THE NAME OF ALLAH  <<     "
# puts "    ELASTIC RESPONSE SPECTRUM ANALYSIS      "
set tStart [clock clicks -milliseconds]
for {set i 1} {$i <= 1000} {incr i} {# loop of mass increment
# UNPUT DATA =================================================================
# units: N, mm, sec
wipe;						# clear opensees model
model basic -ndm 2 -ndf 3;				# 2 dimensions, 3 dof per node
file mkdir SPECTRUM/E_SPEC; 					# create data directory

set B 200; #Width of Column Section
set Ela 23500;# Modulus of Elasticity of Section
set Le 3000;# Length of Column
set Dpr 0.02;#Damping ratio
set Ma [expr 0.2*$i];
set A [expr $B*$B];# area of Section
set I [expr ($B*$B*$B*$B)/12];# Moment Intria of Section
set Ke [expr (3*$Ela*$I)/pow($Le,3)];# Sitffness of Structure
set PERIOD [expr  2*3.1415*sqrt($Ma/$Ke)];# Sitffness of Structure
 puts "INCREMENT: $i -> PERIOD: $PERIOD"

# define GEOMETRY ===========================================================
# nodal coordinates:
node 1 0 0;					# node#, X Y
node 2 0 $Le

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			# node DX DY RZ

# nodal masses:
mass 2 $Ma 1.e-9 0.;					# node#, Mx My Mz, Mass=Weight/g.

# Define ELEMENTS ===========================================================
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
geomTransf Linear 1;  		# associate a tag to transformation

# connectivity: (make A very large, 10e6 times its actual value)
element elasticBeamColumn 1 1 2  $A $Ela $I 1;	# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag

# Define RECORDERS ==========================================================
recorder EnvelopeNode -file SPECTRUM/E_SPEC/DFree_$i.txt -time -node 2 -dof 1 disp;			# displacements of free nodes
recorder EnvelopeNode -file SPECTRUM/E_SPEC/VFree_$i.txt -time -node 2 -dof 1 vel;			# vel of free nodes
recorder EnvelopeNode -file SPECTRUM/E_SPEC/AFree_$i.txt -time -node 2 -dof 1 accel;			# accel of free nodes
recorder EnvelopeNode -file SPECTRUM/E_SPEC/RBase_$i.txt -time -node 1 -dof 1 reaction;			# support reaction

# define GRAVITY ============================================================
pattern Plain 1 Linear {
   load 2 0. 0. 0.;			# node#, FX FY MZ --  superstructure-weight
}
constraints Plain;     				# how it handles boundary conditions
numberer Plain;					# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;				# how to store and solve the system of equations in the analysis
test NormDispIncr 1.0e-8 6 ; 				# determine if convergence has been achieved at the end of an iteration step
algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
integrator LoadControl 0.1;				# determine the next time step for an analysis, # apply gravity in 10 steps
analysis Static					# define type of analysis static or transient
analyze 10;					# perform gravity analysis
loadConst -time 0.0;				# hold gravity constant and restart time

# DYNAMIC ground-motion analysis ===========================================
# create load pattern
#set accelSeries "Series -dt 0.01 -filePath BM68elc.acc -factor 1";	# define acceleration vector from file (dt=0.01 is associated with the input file gm)
set accelSeries "Series -dt 0.01 -filePath Northridge_EQ.acc -factor 1";	# define acceleration vector from file (dt=0.01 is associated with the input file gm)
pattern UniformExcitation 2 1 -accel $accelSeries;		# define where and how (pattern tag, dof) acceleration is applied
rayleigh 0. 0. 0. [expr 2*$Dpr/pow([eigen 1],0.5)];		# set damping based on first eigen mode

# create the analysis
wipeAnalysis;					# clear previously-define analysis parameters
constraints Plain;     				# how it handles boundary conditions
numberer Plain;					# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;					# how to store and solve the system of equations in the analysis
test NormDispIncr 1.0e-8 10;				# determine if convergence has been achieved at the end of an iteration step
algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
integrator Newmark 0.5 0.25 ;			# determine the next time step for an analysis
analysis Transient;					# define type of analysis: time-dependent
analyze 5736 0.01;					# apply 1000 0.02-sec time steps in analysis
}
set tEnd [clock clicks -milliseconds]
set duration [expr $tEnd-$tStart]
puts "Anaysis Duration: $duration"