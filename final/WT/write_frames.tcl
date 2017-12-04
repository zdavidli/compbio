mol load psf ionized.psf dcd prod.dcd
set nf [molinfo 0 get numframes]
for {set i 0 } {$i < $nf} {incr i} {
[atomselect top protein frame $i] writepdb frames/frame_$i.pdb 
}
