# ==============================================================================
# COMPOSITE_FRP_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA_3D
# Creates a circular FRP composite section with confined concrete, steel pipe,
# and steel rebars in OpenSees (Tcl version).
#
# Usage:
#   COMPOSITE_FRP_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA secTag RI RO COVER \
#       fc Kfc THICKNESS STEEL_TYPE STEEL_DENSITY CONCRETE_DENSITY plotFlag \
#       SECTION_HEIGHT_VAR ELE_MASS_VAR
#
# Input arguments:
#   secTag         – integer, tag for the fiber section
#   RI             – inner radius of concrete core (mm)
#   RO             – outer radius of concrete core (mm) = radius to steel pipe inner face
#   COVER          – concrete cover thickness (mm) – used only to compute rc = RO - COVER
#   fc             – unconfined concrete compressive strength (positive value, N/mm²)
#   Kfc            – ratio confined/unconfined concrete strength
#   THICKNESS      – thickness of the external steel pipe (mm)
#   STEEL_TYPE     – "ELASTIC" or "INELASTIC" (for rebars and I‑section, though I‑section not used)
#   STEEL_DENSITY  – density of steel (kg/mm³) – actually used as kg/mm³? typical: 7.85e-6 kg/mm³
#   CONCRETE_DENSITY – density of concrete (kg/mm³) – typical: 2.4e-6 kg/mm³
#
# Output (via upvar):
#   SECTION_HEIGHT_VAR – variable name that receives the outer diameter (mm)
#   ELE_MASS_VAR       – variable name that receives the mass per unit length (kg/mm)
#
# Notes:
#   - All lengths are in mm, forces in N, stresses in N/mm².
#   - The steel I‑section (crossed shapes) from the Python plot are NOT added to the
#     structural section because the original code omitted them. Only the concrete core,
#     concrete cover, FRP pipe, and rebars are included.
#   - Rebar positions and diameters are hard‑coded as in the Python function.
#   - The FRP tube is modeled as an elastic material (ElasticUniaxialMaterial).
#   - Inelastic steel uses the Hysteretic material with the same parameters.
#
# Originally written for Python by Salar Delavar Ghashghaei (Qashqai)
# Tcl translation by request.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
# ==============================================================================

proc COMPOSITE_FRP_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA_3D {secTag RI RO COVER fc Kfc THICKNESS STEEL_TYPE STEEL_DENSITY CONCRETE_DENSITY} {
    # --------------------------------------------------------------------------
    # 1. Material parameters (units: N, mm)
    # --------------------------------------------------------------------------
    # Unconfined concrete
    set fcU [expr {-1.0 * $fc}]                     # negative compressive strength
    set EcU [expr {4700.0 * sqrt(-$fcU)}]           # N/mm²
    set ec0U [expr {2.0 * $fcU / $EcU}]
    set fcUU [expr {0.2 * $fcU}]
    set ecuU [expr {5.0 * $ec0U}]
    set LambdaU 0.1
    set ftU [expr {0.7 * sqrt(-$fcU)}]
    set EtsU [expr {$ftU / abs($ec0U)}]

    # Confined concrete
    set fcC [expr {$Kfc * $fcU}]                    # negative
    set EcC [expr {4700.0 * sqrt(-$fcC)}]
    set ec0C [expr {2.0 * $fcC / $EcC}]
    set fcUC [expr {0.65 * $fcC}]
    set ecuC [expr {15.0 * $ec0C}]
    set LambdaC 0.1
    set ftC [expr {0.7 * sqrt(-$fcC)}]
    set EtsC [expr {$ftC / abs($ec0C)}]

    # Steel rebar properties (same for all rebars)
    set fy 400.0
    set Es 200000.0
    set ey [expr {$fy / $Es}]
    set fu [expr {1.1818 * $fy}]
    set esu 0.09
    set Esh [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs [expr {$Esh / $Es}]

    # Steel I‑section properties (not used in section, but defined for material)
    set fyI 240.0
    set EsI 200000.0
    set eyI [expr {$fyI / $EsI}]
    set fuI [expr {1.1818 * $fyI}]
    set esuI 0.25
    set EshI [expr {($fuI - $fyI) / ($esuI - $eyI)}]
    set BsI [expr {$EshI / $EsI}]

    # FRP tube properties
    set fyF 4800.0
    set eyF 0.02
    set EsF [expr {$fyF / $eyF}]

    # Material tags (unique per section)
    set coreTag   [expr {$secTag + 1000}]
    set coverTag  [expr {$secTag + 2000}]
    set steelTag  [expr {$secTag + 3000}]
    set steelITag [expr {$secTag + 4000}]
    set steelFTag [expr {$secTag + 5000}]

    # --------------------------------------------------------------------------
    # 2. Define uniaxial materials
    # --------------------------------------------------------------------------
    # Steel rebar material (depending on STEEL_TYPE)
    if {[string equal -nocase $STEEL_TYPE "ELASTIC"]} {
        uniaxialMaterial Elastic $steelTag $Es
        uniaxialMaterial Elastic $steelITag $EsI
    } else {
        # INELASTIC: Hysteretic material with symmetric tension/compression
        set pinchX 0.8
        set pinchY 0.5
        set damage1 0.0
        set damage2 0.0
        set beta 0.1
        uniaxialMaterial Hysteretic $steelTag \
            $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
            [expr {-$fy}] [expr {-$ey}] [expr {-$fu}] [expr {-$esu}] [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
            $pinchX $pinchY $damage1 $damage2 $beta
        uniaxialMaterial Hysteretic $steelITag \
            $fyI $eyI $fuI $esuI [expr {0.2*$fuI}] [expr {1.1*$esuI}] \
            [expr {-$fyI}] [expr {-$eyI}] [expr {-$fuI}] [expr {-$esuI}] [expr {-0.2*$fuI}] [expr {-1.1*$esuI}] \
            $pinchX $pinchY $damage1 $damage2 $beta
    }

    # Concrete materials
    uniaxialMaterial Concrete02 $coreTag $fcC $ec0C $fcUC $ecuC $LambdaC $ftC $EtsC
    uniaxialMaterial Concrete02 $coverTag $fcU $ec0U $fcUU $ecuU $LambdaU $ftU $EtsU

    # FRP tube (elastic with small damping)
    uniaxialMaterial Elastic $steelFTag $EsF

    # --------------------------------------------------------------------------
    # 3. Create fiber section
    # --------------------------------------------------------------------------
    section Fiber $secTag -GJ 1.0e7

    # Radii for patches
    set rc [expr {$RO - $COVER}]       # radius of concrete core (confined part)
    set rPipeOuter [expr {$RO + $THICKNESS}]

    # Numbers of subdivisions (as in Python)
    set nfCoreR 8
    set nfCoreT 8
    set nfCoverR 4
    set nfCoverT 8

    # Core concrete (confined)
    patch circ $coreTag $nfCoreT $nfCoreR 0.0 0.0 $RI $rc 0.0 360.0

    # Cover concrete (unconfined)
    patch circ $coverTag $nfCoverT $nfCoverR 0.0 0.0 $rc $RO 0.0 360.0

    # FRP pipe (outer tube)
    patch circ $steelFTag $nfCoverT $nfCoverR 0.0 0.0 $RO $rPipeOuter 0.0 360.0

    # --------------------------------------------------------------------------
    # 4. Add steel rebars (same coordinates and diameters as Python code)
    #    Each: fiber x y area matTag
    # --------------------------------------------------------------------------
    # Hard‑coded list: (diameter, x, y)
    set rebarList {
        {25.0  200.0000   0.0000}
        {16.0  184.7759  76.5367}
        {25.0    0.0000 200.0000}
        {16.0 -184.7759  76.5367}
        {25.0 -200.0000   0.0000}
        {16.0 -184.7759 -76.5367}
        {25.0    0.0000 -200.0000}
        {16.0  184.7759 -76.5367}
        {16.0  141.4214 -141.4214}
        {16.0   76.5367 -184.7759}
        {16.0  -76.5367 -184.7759}
        {16.0 -141.4214 -141.4214}
        {16.0  -76.5367  184.7759}
        {16.0  141.4214  141.4214}
        {16.0 -141.4214  141.4214}
        {16.0   76.5367  184.7759}
    }

    foreach rebar $rebarList {
        lassign $rebar dia x y
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # --------------------------------------------------------------------------
    # 5. Compute section height (outer diameter) and mass per unit length
    # --------------------------------------------------------------------------
    set SECTION_HEIGHT [expr {2.0 * $rPipeOuter}]

    set CONCRETE_AREA   [expr {3.141592653589793 * ($RO*$RO - $RI*$RI)}]
    set STEEL_PIPE_AREA [expr {3.141592653589793 * ($rPipeOuter*$rPipeOuter - $RO*$RO)}]

    # Mass per unit length (kg/mm)
    set ELE_MASS [expr {$CONCRETE_AREA * $CONCRETE_DENSITY + $STEEL_PIPE_AREA * $STEEL_DENSITY}]

    puts "Section height = $SECTION_HEIGHT mm"

    # Return both values as a Tcl list
    return [list $SECTION_HEIGHT $ELE_MASS]
}

# Example: create a section with tag 1
#set secTag 1
#set RI 0
#set RO 200
#set COVER 40
#set fc 30          ;# unconfined strength (positive)
#set Kfc 1.3
#set THICKNESS 10
#set STEEL_TYPE "INELASTIC"
#set STEEL_DENSITY 7.85e-6     ;# kg/mm³ (7850 kg/m³)
#set CONCRETE_DENSITY 2.4e-6   ;# kg/mm³ (2400 kg/m³)

#COMPOSITE_FRP_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA_3D $secTag $RI $RO $COVER $fc $Kfc $THICKNESS $STEEL_TYPE $STEEL_DENSITY $CONCRETE_DENSITY

#puts "Section outer diameter = $secHeight mm"
#puts "Mass per unit length   = $massPerUnit kg/mm"