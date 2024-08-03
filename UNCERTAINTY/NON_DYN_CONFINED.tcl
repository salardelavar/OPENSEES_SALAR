#  #######################################################################
#  #                           IN THE NAME OF ALLAH                      #
#  # NONLINEAR DYNAMIC ANALYSIS ON CONCRETE CONFINED SECTION COLUMN      #
#  # WITH UNCERTAINTY CONDITIONS USING PROBABILITY DISTRIBUTION FUNCTION #
#  #                        MONTE CARLO METHOD                           #
#  #---------------------------------------------------------------------#
#  #            THIS PROGRAM WRITTEN BY SALAR DELAVAR QASHQAI            #
#  #                 EMAIL: salar.d.ghashghaei@gmail.com                 #
#  #######################################################################

# Math Functions
source SALAR_MATH_FUN.tcl

set tStart [clock clicks -milliseconds]
set NUM_INCREMENT 1000; #Number of Analysis Increment 
for {set i 1} {$i <= $NUM_INCREMENT} {incr i} {# loop of mass increment
# SET UP ----------------------------------------------------------------------------
# units: N, mm, sec
wipe;					# clear memory of all past model definitions
file mkdir DATA/NON_DYN_CONFINED; 				# create data directory
model BasicBuilder -ndm 2 -ndf 3;		# Define the model builder, ndm=#dimension, ndf=#dofs


# define GEOMETRY ========================================================
# define section geometry
set B [randUniform 199 201]; #Width of Column Section
set Ela 23500;# Modulus of Elasticity of Section
set Le [randUniform 2985 3005];# Length of Column
set Dpr [randTriangular 0.018 0.02 0.022];#Damping ratio
set RZD [randUniform 15.9 16.1];# Column Rebar Size Diameter
set RANDOM_MASS [randTriangular 0.15 0.18 0.22]; #Random Mass
#set Ma [expr 0.2];
set A [expr $B*$B];# area of Section
set I [expr ($B*$B*$B*$B)/12];# Moment Intria of Section
set Ke [expr (3*$Ela*$I)/pow($Le,3)];# Sitffness of Structure
set PERIOD [expr  2*3.1415*sqrt($RANDOM_MASS/$Ke)];# Sitffness of Structure
 puts "INCREMENT: $i -> PERIOD: $PERIOD"


# nodal coordinates:
node 1 0 0;			# node#, X, Y
node 2 0 $Le 		

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			# node DX DY RZ

# nodal masses:
mass 2 $RANDOM_MASS  1e-9 0.;		# node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes

# Define ELEMENTS & SECTIONS ========================================================
set ColSecTag 1;			# assign a tag number to the column section	
# define section geometry
set coverCol [randUniform 23 27];			# Column cover to reinforcing steel NA.
set numBarsCol 8;			# number of longitudinal-reinforcement bars in column. (symmetric top & bot)
set barAreaCol [expr (3.1415*$RZD*$RZD)/4];	# area of longitudinal-reinforcement bars
# MATERIAL parameters -------------------------------------------------------------------
set IDconcCore 1; 				# material ID tag -- confined core concrete
set IDconcCover 2; 				# material ID tag -- unconfined cover concrete
set IDreinf 3; 				# material ID tag -- reinforcement
# nominal concrete compressive strength
set RANDOM_FC [randTriangular 22 25 27]
set fc 		[expr -$RANDOM_FC];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set Ec 		[expr 4700*sqrt(-$fc)];	# Concrete Elastic Modulus
# confined concrete
set Kfc 		[randUniform 1.28 1.32];			# ratio of confined to unconfined concrete strength
set fc1C 		[expr $Kfc*$fc];		# CONFINED concrete (mander model), maximum stress
set eps1C	[expr 2.*$fc1C/$Ec];	# strain at maximum stress 
set fc2C 		[expr 0.2*$fc1C];		# ultimate stress
set eps2C 	[expr 5*$eps1C];		# strain at ultimate stress 
# unconfined concrete
set fc1U 		$fc;			# UNCONFINED concrete (todeschini parabolic model), maximum stress
set eps1U	-[randUniform 0.0025 0.003];			# strain at maximum strength of unconfined concrete
set fc2U 		[expr 0.2*$fc1U];		# ultimate stress
set eps2U	-[randUniform 0.008 0.012];			# strain at ultimate stress
set lambda 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
set ftC [expr -0.55*$fc1C];		# tensile strength +tension
set ftU [expr -0.55*$fc1U];		# tensile strength +tension
set Ets [expr $ftU/0.002];		# tension softening stiffness
# -----------
set RANDOM_FY [randTriangular 3500 3800 4000]
set RANDOM_CY [randTriangular 0.015 0.018 0.02]
set Fy $RANDOM_FY;		# STEEL yield stress
set Cy $RANDOM_CY;		# STEEL yield strain
set Es [expr $Fy/$Cy];		# modulus of steel
set Bs		0.01;			# strain-hardening ratio 
set R0 18;				# control the transition from elastic to plastic branches
set cR1 0.925;				# control the transition from elastic to plastic branches
set cR2 0.15;				# control the transition from elastic to plastic branches
uniaxialMaterial Concrete02 $IDconcCore $fc1C $eps1C $fc2C $eps2C $lambda $ftC $Ets;	# build core concrete (confined)
uniaxialMaterial Concrete02 $IDconcCover $fc1U $eps1U $fc2U $eps2U $lambda $ftU $Ets;	# build cover concrete (unconfined)
uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2;				# build reinforcement material

# FIBER SECTION properties -------------------------------------------------------------
# RC section: 
   set coverY [expr $B/2.0];	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coverZ [expr $B/2.0];	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
   set coreY [expr $coverY-$coverCol];	# The distance from the section z-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
   set coreZ [expr $coverZ-$coverCol ];	# The distance from the section y-axis to the edge of the core concrete --  edge of the core concrete/inner edge of cover concrete
   set nfCoreY 16;			# number of fibers for concrete in y-direction -- core concrete
   set nfCoreZ 4;			# number of fibers for concrete in z-direction
   set nfCoverY 16;			# number of fibers for concrete in y-direction -- cover concrete
   set nfCoverZ 4;			# number of fibers for concrete in z-direction
   section fiberSec $ColSecTag   {;	# Define the fiber section
	# Define the core patch
	patch quadr $IDconcCore $nfCoreZ $nfCoreY -$coreY $coreZ -$coreY -$coreZ $coreY -$coreZ $coreY $coreZ
      
	# Define the four cover patches
	patch quadr $IDconcCover 1 $nfCoverY -$coverY $coverZ -$coreY $coreZ $coreY $coreZ $coverY $coverZ
	patch quadr $IDconcCover 1 $nfCoverY -$coreY -$coreZ -$coverY -$coverZ $coverY -$coverZ $coreY -$coreZ
	patch quadr $IDconcCover $nfCoverZ 1 -$coverY $coverZ -$coverY -$coverZ -$coreY -$coreZ -$coreY $coreZ
	patch quadr $IDconcCover $nfCoverZ 1 $coreY $coreZ $coreY -$coreZ $coverY -$coverZ $coverY $coverZ

	# Define reinfocement layers
	layer straight $IDreinf $numBarsCol $barAreaCol  $coreY $coreZ  $coreY -$coreZ;	# top layer reinforcement
	layer straight $IDreinf $numBarsCol $barAreaCol -$coreY $coreZ -$coreY -$coreZ;	# bottom layer reinfocement
    };	# end of fibersection definition

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
set ColTransfTag 1; 			# associate a tag to column transformation
geomTransf Linear $ColTransfTag  ; 	

# element connectivity:
set numIntgrPts 5;								# number of integration points for force-based element
element nonlinearBeamColumn 1 1 2 $numIntgrPts $ColSecTag $ColTransfTag;	# self-explanatory when using variables

# Define RECORDERS ========================================================
recorder EnvelopeNode -file DATA/NON_DYN_CONFINED/DFree_$i.txt -time -node 2 -dof 1 disp;			# displacements of free nodes
recorder EnvelopeNode -file DATA/NON_DYN_CONFINED/VFree_$i.txt -time -node 2 -dof 1 vel;			# vel of free nodes
recorder EnvelopeNode -file DATA/NON_DYN_CONFINED/AFree_$i.txt -time -node 2 -dof 1 accel;			# accel of free nodes
recorder EnvelopeNode -file DATA/NON_DYN_CONFINED/RBase_$i.txt -time -node 1 -dof 1 reaction;			# support reaction
# define GRAVITY ========================================================
pattern Plain 1 Linear {
   load 2 0 0 0
}

# Gravity-analysis parameters -- load-controlled static analysis
set Tol 1.0e-8;			# convergence tolerance for test
constraints Plain;     		# how it handles boundary conditions
numberer Plain;			# renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;		# how to store and solve the system of equations in the analysis
test NormDispIncr $Tol 6 ; 		# determine if convergence has been achieved at the end of an iteration step
algorithm Newton;			# use Newton's solution algorithm: updates tangent stiffness at every iteration
set NstepGravity 10;  		# apply gravity in 10 steps
set DGravity [expr 1./$NstepGravity]; 	# first load increment;
integrator LoadControl $DGravity;	# determine the next time step for an analysis
analysis Static;			# define type of analysis static or transient
analyze $NstepGravity;		# apply gravity
# ------------------------------------------------- maintain constant gravity loads and reset time to zero
loadConst -time 0.0

puts "Model Built"

# DYNAMIC EQ ANALYSIS ========================================================
# Uniform Earthquake ground motion (uniform acceleration input at all support nodes)
set GMdirection 1;				# ground-motion direction
set GMfile "Northridge_EQ.acc" ;			# ground-motion filenames
set GMfact 1;				# ground-motion scaling factor

# set up ground-motion-analysis parameters
set DtAnalysis	0.01;	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr 10.];	# maximum duration of ground-motion analysis -- should be 50*$sec

# DYNAMIC ANALYSIS PARAMETERS
constraints Transformation ; 
numberer Plain
system SparseGeneral -piv
set Tol 1.e-8;              # Convergence Test: tolerance
set maxNumIter 10;          # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
set printFlag 0;            # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
set TestType EnergyIncr;	# Convergence-test type
test $TestType $Tol $maxNumIter $printFlag;
set algorithmType ModifiedNewton 
algorithm $algorithmType;        
set NewmarkGamma 0.5;	# Newmark-integrator gamma parameter (also HHT)
set NewmarkBeta 0.25;	# Newmark-integrator beta parameter
integrator Newmark $NewmarkGamma $NewmarkBeta 


analysis Transient

# define DAMPING========================================================
# apply Rayleigh DAMPING from $xDamp
# D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial
set xDamp $Dpr;				# 2% damping ratio
set lambda [eigen 1]; 			# eigenvalue mode 1
set omega [expr pow($lambda,0.5)];
set alphaM 0.;				# M-prop. damping; D = alphaM*M
set betaKcurr 0.;         			# K-proportional damping;      +beatKcurr*KCurrent
set betaKcomm [expr 2.*$xDamp/($omega)];   	# K-prop. damping parameter;   +betaKcomm*KlastCommitt
set betaKinit 0.;         			# initial-stiffness proportional damping      +beatKinit*Kini
# define damping
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm; 				# RAYLEIGH damping

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
# Uniform EXCITATION: acceleration input
set IDloadTag 400;			# load tag
set dt 0.01;			# time step for input ground motion
set GMfatt 1.0;			# data in input file is in g Unifts -- ACCELERATION TH
set AccelSeries "Series -dt $dt -filePath $GMfile -factor  $GMfatt";			# time series information
pattern UniformExcitation  $IDloadTag  $GMdirection -accel  $AccelSeries  ;		# create Unifform excitation

set Nsteps 5736;#[expr int($TmaxAnalysis/$DtAnalysis)];
set ok [analyze $Nsteps $DtAnalysis];			# actually perform analysis; returns ok=0 if analysis was successful

if {$ok != 0} {      ;					# if analysis was not successful.
	# change some analysis parameters to achieve convergence
	# performance is slower inside this loop
	#    Time-controlled analysis
	set ok 0;
	set controlTime [getTime];
	while {$controlTime < $TmaxAnalysis && $ok == 0} {
		set ok [analyze 1 $DtAnalysis]
		set controlTime [getTime]
		set ok [analyze 1 $DtAnalysis]
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test NormDispIncr   $Tol 1000  0
			algorithm Newton -initial
			set ok [analyze 1 $DtAnalysis]
			test $TestType $Tol $maxNumIter  0
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmType
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch .8
			set ok [analyze 1 $DtAnalysis]
			algorithm $algorithmType
		}
	}
};      # end if ok !0


puts "Ground Motion Done. End Time: [getTime]"
}
set tEnd [clock clicks -milliseconds]
set duration [expr $tEnd-$tStart]
puts "Anaysis Duration: $duration (s)"