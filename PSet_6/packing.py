import rosetta
rosetta.init()
from toolbox import generate_resfile_from_pose 

static_pose=pose_from_pdb('1KJG.clean.pdb')
#get full atom score function
scorefxn=get_fa_scorefxn()
move_map=MoveMap()
#disallow backbone movements
move_map.set_bb(False)
#allow chi angle movements
move_map.set_chi(True)
min_mover=MinMover()
min_mover.movemap(move_map)
#use Davidson-Fletcher-Powell minimization
min_mover.min_type('dfpmin')

def make_mover(pose,scorefxn,pack,mut_num):
	'''Creates repacking mover based on a given pose and score function. 'pack' variable 
	specifies whether or not packing mover will allow mutations or not. If pack='mutant' the mut_num
	is the residue number of the P chain that will be mutated.'''
	
	#pdb numbering for active site residues
	resi_lst=[8,23,25,29,30,32,45,47,50,53,82,84]
	
	#list of residue rosetta numbers for active site
	refined_r_lst=[]
	
	#for each chain
	for i in range(2):
		if i==0:
			chain='A'
		else:
			chain='B'
		#for every active site residue
		for j in range(len(resi_lst)):
		
			#update list of rosetta numbers for each active site residue in chain A and B
			refined_r_lst.append(pose.pdb_info().pdb2pose(chain,resi_lst[j]))
			
	#p chain residues		
	for k in range(8):
		#update list of rosetta numbers for residues 2-9 on chain P
		refined_r_lst.append(pose.pdb_info().pdb2pose('P',k+2))

	#create standard packer object
	task_pack = standard_packer_task(pose)

	if pack=='native':
		#only repacking (turn off design/mutations)
		task_pack.restrict_to_repacking()
		#fix all residues
		task_pack.temporarily_fix_everything()
	
	for h in range(len(refined_r_lst)):
		#unfix residue
		task_pack.temporarily_set_pack_residue(refined_r_lst[h],True)
		
	#if allowing residue mutation on p chain
	if pack=='mutant':
		#read in residue allowed to be mutated from resfile
		#task_pack.design_residue(mut_num)
		parse_resfile(pose, task_pack, '/Users/Drew/Desktop/Prot_Predict_and_Structure/PSet_6/1KJG-P'+str(mut_num)+'.resfile')
		
	pack_mover=PackRotamersMover(scorefxn,task_pack)
	
	return pack_mover
	
			


def pack_chain(pack):
	'''Repacks side chain. Takes in 'pack' variable that specifies whether or not to pack a native
	structure or one with mutations allowed'''
	all=[]
	if pack=='native':
		#list of decoy scores
		scores=[]
		#make 10 decoys
		for i in range(1):
			##reset pose made from PDB file
			iter_pose=pose_from_pdb('1KJG.clean.pdb')
			#make pack_mover
			pack_mover=make_mover(iter_pose,scorefxn,pack,0)
			#repack pose
			pack_mover.apply(iter_pose)
			#use minimize mover
			min_mover.apply(iter_pose)
			#add pose score to list of scores
			scores.append(scorefxn(iter_pose))
			scorefxn.show(iter_pose)
		#return lowest score
		return min(scores)
	
	
	if pack=='mutant':
		scores=[]
		#for all 8 residues to be mutated on peptide
		for k in range(8):
			low_score=0
			low_mut=str(k+2)
			#scores of each decoy
			#make 10 decoys
			for i in range(1):
				#reset pose from PDB file
				#iter_pose=pose_from_pdb('1KJG.clean.pdb')
				iter_pose=Pose()
				iter_pose.assign(static_pose)
				#make pack_mover
				pack_mover=make_mover(iter_pose,scorefxn,pack,k+2)
				#repack/design pose
				pack_mover.apply(iter_pose)
				#use minimize mover
				min_mover.apply(iter_pose)
				
				if scorefxn(iter_pose)<low_score:
					#residue after design 
					mut_resi=iter_pose.residue(static_pose.pdb_info().pdb2pose('P',k+2)).name()
					#residue before design
					orig_resi=static_pose.residue(static_pose.pdb_info().pdb2pose('P',k+2)).name()
					#low score
					low_score=scorefxn(iter_pose)
					scorefxn.show(iter_pose)
					#mutation that causes lowest score
					low_mut=orig_resi+str(k+2)+mut_resi
					print low_mut
					
			scores.append((low_score,low_mut))
				
			#add lowest score from residue to list of all residue low scores as tuple
			#tuple has format: (low score,residue number) 
		return scores
	
def main():
	'''Repacks side chains on a native structure and mutant structures of 1KJG. Prints 
	low energies of native and mutants'''
	native=pack_chain('native')
	designs=pack_chain('mutant')
	print 'Native Low Score: '+str(native)
	
	for j in range(len(designs)):
		if float(designs[j][0])==float(native):
			print 'Residue '+str(designs[j][1][3])+' Low Score: '+str(designs[j][0])+ ' No Score Improvement'
		elif float(designs[j][0])<float(native):
			print 'Residue '+str(designs[j][1][3])+' Low Score: '+str(designs[j][0])+' Top Mutation '+designs[j][1]
			
main()