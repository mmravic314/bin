## Find best fit of helices in helical pair cluster 4 (Zhang et al 2015, Structure) to DF1 helix 2 (Chain A)
## This runs 1 backbone atom superposition at a time

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, os, pymol, time

pymol.finish_launching()

if len(sys.argv)<3 or sys.argv[1][-3:]!='pdb': 
	print "ERROR!! arguments must be pdb files and user must supply two \ncmd.py staticPath mobilePath"
	sys.exit()

staticStructurePath = sys.argv[1]
staticStructureName = os.path.basename( staticStructurePath ).split('.')[0]
mobileStructurePath = sys.argv[2]
mobileStructureName = os.path.basename( mobileStructurePath ).split('.')[0]
 
# Load Structures
 
pymol.cmd.load(mobileStructurePath, mobileStructureName)
pymol.cmd.load(staticStructurePath, staticStructureName)


# CEAlign STATIC, MOBILE
# CEAlign produces same alignment as the complex SUPER below
#pymol.cmd.do("run /home/joao/Software/cealign-0.9/cealign.py") # Import Module
#output=pymol.cmd.do("align %s, %s" %(staticStructureName, mobileStructureName))
time.sleep(1) # Dunno why, but if I don't wait, structures do not align properly..

print ">PyMOL align %s, %s" %(staticStructureName, mobileStructureName)

output=pymol.cmd.align(staticStructureName, mobileStructureName, cutoff=2.0,
          cycles=5,  #gap=-10.0, extend=-0.5,
          #max_gap=50, object=None, 
          mobile_state=0, target_state=0,  quiet=1,
           max_skip=0,  transform=1, reset=0 )
time.sleep(10)
print output[0], "over", output[1], 'atoms'

# Save Superimposition
#print output, 'HIIIIIIIIIIIIIIIII'
# save(file, selection, state (0 default), format)
#pymol.cmd.save("%s_%s.pdb" %(mobileStructureName, staticStructureName), mobileStructureName, 0, 'pdb')
 
## SUPER - old
#pymol.cmd.super((staticStructureName and (resn ZN around 5 and (resn CYS or resn HIS))), (mobileStructureName and (resn ZN around 5 and (resn CYS or resn HIS))))
#pymol.cmd.save("%s_%s_SUPER.pdb" %(mobileStructureName, staticStructureName), mobileStructureName, 0, 'pdb')
#print sys.argv[3], 'max % sequence ID of helical chains'
#print  
# Get out!
pymol.cmd.quit()

