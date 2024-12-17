set proteinsel "8a"
set copynum "01"
set files [glob -directory "/Users/masauer2/Library/CloudStorage/Box-Box/Summerinternship_2024/WNT-WLS-project/Data_Backup/wnt${proteinsel}/copy${copynum}/dcd_files" -- *]
set nFiles [llength $files]
for {set i 1} {$i <= $nFiles} {incr i} {
if {$i < 10} {
set selstr "00${i}"
} elseif {$i >= 10 && $i < 100} {
set selstr "0${i}"
} else {
set selstr ${i}
}
mol new /Users/masauer2/Library/CloudStorage/Box-Box/Summerinternship_2024/WNT-WLS-project/Data_Backup/wnt${proteinsel}/copy${copynum}/Wnt${proteinsel}WlsPc_copy_${copynum}.psf
mol addfile /Users/masauer2/Library/CloudStorage/Box-Box/Summerinternship_2024/WNT-WLS-project/Data_Backup/wnt${proteinsel}/copy${copynum}/dcd_files/Wnt${proteinsel}WlsPc_copy_${copynum}_run_${selstr}.dcd waitfor all
set proteinonly [atomselect top "segname PROA or segname PROB"]
animate write dcd "output/wnt${proteinsel}/Wnt${proteinsel}WlsPc_copy_${copynum}_run_${selstr}.dcd" sel $proteinonly
mol delete top
}
mol new /Users/masauer2/Library/CloudStorage/Box-Box/Summerinternship_2024/WNT-WLS-project/Data_Backup/wnt${proteinsel}/copy${copynum}/Wnt${proteinsel}WlsPc_copy_${copynum}.pdb
set proteinonly [atomselect top "segname PROA or segname PROB"]
$proteinonly writepdb output/wnt${proteinsel}/Wnt${proteinsel}WlsPc_copy_${copynum}.pdb
animate write psf output/wnt${proteinsel}/total_${copynum}.psf sel $proteinonly
