set proteinsel "4"
set files [glob -directory "wnt${proteinsel}" -- *]
set nFiles [llength $files]
mol new wnt${proteinsel}/total.psf
for {set i 1} {$i < $nFiles} {incr i} {
if {$i < 10} {
set selstr "00${i}"
} elseif {$i >= 10 && $i < 100} {
set selstr "0${i}"
} else {
set selstr ${i}
}
mol addfile wnt${proteinsel}/Wnt${proteinsel}WlsPc_copy_01_run_${selstr}.dcd waitfor all
}
set nf [molinfo top get numframes]
for {set i 0} {$i < 5} {incr i} { 
    # Load the frame
    animate goto $i
    set sel [atomselect top all frame $i]
    # Calculate the center of mass
    set com [measure center $sel]
    set transfo [transoffset [vecinvert $com]]
    puts $transfo
    # Move the molecule to the center of mass
    $sel moveby $transfo
    set newcom [measure center $sel]
    puts $newcom
} 
