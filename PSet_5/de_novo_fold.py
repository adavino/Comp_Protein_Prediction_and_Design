from rosetta import *
rosetta.init()
import math
import random

#domain sequence to be folded
seq='SYKGEKIGQGKANATAWLKDNPETAKEIEKKVRELLLSNPNSTPDFSVDDSEGVAETNEDF'

pose=pose_from_sequence(seq)
kT=1.0
n_moves=5
frag_runs=10000
ref_runs=10000

scorefxn=get_fa_scorefxn()

movemap=MoveMap()
movemap.set_bb(True)
	
#make small and shear mover
small_mover=SmallMover(movemap,kT,n_moves)
shear_mover=ShearMover(movemap,kT,n_moves)
	
#make minimization mover
min_mover=MinMover()
mm_40_60=MoveMap()
mm_40_60.set_bb_true_range(40,60)
min_mover.movemap(mm_40_60)
min_mover.score_function(scorefxn)

#combine small and min mover into one
small_and_min=SequenceMover()
small_and_min.add_mover(small_mover)
small_and_min.add_mover(min_mover)

#combine small and min mover into one
shear_and_min=SequenceMover()
shear_and_min.add_mover(shear_mover)
shear_and_min.add_mover(min_mover)

#make monte carlo mover
mc=MonteCarlo(pose,scorefxn,kT)

#make 9 aa frag mover
frag9_move=MoveMap()
frag9_move.set_bb(True)
fragset9=ConstantLengthFragSet(9)
fragset9.read_fragment_file("aat000_09_05.200_v1_3.txt")
move_9mer=ClassicFragmentMover(fragset9,frag9_move)

#combine frag mover and MC mover to accept/reject each move
trial_9mer=TrialMover(move_9mer,mc)
#execute n=frag_runs trials when run
rep_9mer=RepeatMover(trial_9mer,frag_runs)

#make 3 aa frag mover
frag3_move=MoveMap()
frag3_move.set_bb(True)
fragset3=ConstantLengthFragSet(3)
fragset3.read_fragment_file("aat000_03_05.200_v1_3.txt")
move_3mer=ClassicFragmentMover(fragset3,frag3_move)
#combine frag mover and MC mover to accept/reject each move
trial_3mer=TrialMover(move_3mer,mc)
#execute n=frag_runs trials when run
rep_3mer=RepeatMover(trial_3mer,frag_runs)

#combine fragment movers into one sequence
seq_mover=SequenceMover()
seq_mover.add_mover(rep_9mer)
seq_mover.add_mover(rep_3mer)

#create mover for refinement (small and shear moves)
refine_mover=SequenceMover()
refine_mover.add_mover(small_mover)
refine_mover.add_mover(shear_mover)

#combine refinement mover and MC mover to accept/reject each move
refine_trial=TrialMover(refine_mover,mc)

#execute 1000 trials when run
rep_refine=RepeatMover(refine_trial,ref_runs)
	
pmm=PyMolMover()
for i in range(3):
	#initialize pose from given sequence
	low_res=Pose()
	high_res=Pose()
	pose=Pose()
	reset_pose=pose_from_sequence(seq)
	pose.assign(reset_pose)
	
	#make fragment insertions
	seq_mover.apply(pose)
	
	#low resolution copy
	low_res.assign(pose)
	print '======'
	print 'trial: ',i
	print 'low res score = ', scorefxn(low_res)
	
	
	for i in range(5):
		#decrease angle for small move by 5 degrees each time through
		small_mover.angle_max('H',(5-i)*5)
		small_mover.angle_max('E',(5-i)*5)
		small_mover.angle_max('L',(5-i)*5)
		
		#decrease angle for shear move by 5 degrees each time through
		shear_mover.angle_max('H',(5-i)*5)
		shear_mover.angle_max('E',(5-i)*5)
		shear_mover.angle_max('L',(5-i)*5)
		
		#make refinements
		rep_refine.apply(pose)
	#high resolution copy
	high_res.assign(pose)
	print 'high res score = ', scorefxn(high_res)

	
	
