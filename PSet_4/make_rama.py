import random
import math
import matplotlib.pyplot as plt
from rosetta import *
rosetta.init()

#get score function from Rosetta
scorefxn=get_fa_scorefxn()


def make_move(resi,which_angle):
	'''Takes in the pose and seq. Uses random module to determine which residue will be altered, which
	angle of that residue and by what angle. Returns these three values as int and float values. '''

	if which_angle==0:
		new_pose.set_phi(resi,(new_pose.phi(resi)+random.gauss(0,25))%360)
		#print new_phi
		#print '======'

	else:
		new_pose.set_psi(resi,(new_pose.psi(resi)+random.gauss(0,25))%360)
		#print new_psi
		#print '====='



def accept_reject(d_score):
	'''Determines whether to accept or reject score. Takes in the difference in pose scores from the
	previous pose to the current pose. Returns boolean statement, False rejects the move, True accepts it '''
	#initialize Accept at False

	if d_score<0:
		#changes Accept to True
		return True
	else:
		#P is probability that move is accepted
		P=math.exp(-1*d_score)
		#R is random float between 0 and 1
		R=random.random()
		#accepts move
		if R<P:
			return True
			
		#rejects move
		else:
			return False


#lists of angles in low_poses, used to make Ramachandran
phi_angles=[]
psi_angles=[]
#number of decoys generated
for j in range(100):
	
	seq='A'*10
	
	#created pose
	low_pose=Pose()
	new_pose=Pose()
	
	pose=pose_from_sequence(seq)
	
	new_pose.assign(pose)
	
	low_energy=scorefxn(pose)
	
	#pmm=PyMOL_Mover()
	#pmm.apply(pose)
	
	#number of angle perturbations on each decoy
	for i in range(50000):
		
		which_angle=random.randint(0,1)
		
		resi=random.randint(1,pose.total_residue())

		make_move(resi,which_angle)
		
		#calculate score before next move
		old_score=scorefxn(pose)

		#calculate score after move
		new_score=scorefxn(new_pose)
		
		#find difference in energy made by move
		d_score=new_score-old_score
		
		#decides whether to accepts or reject move
		accept=accept_reject(d_score)
		
		#if it should not be accepted
		
		if accept==True:
			#current pose updated
			pose.assign(new_pose)
		else:
			#new_pose reset to last accepted pose
			new_pose.assign(pose)
		
		if new_score<low_energy:
			#changes lowest energy recorded
			low_energy=new_score
			#designates current pose as 'low_pose'
			low_pose.assign(pose)
	
	#take in phi and psi angles of lowest energy poses for each trajectory
	for k in range(low_pose.total_residue()):
		phi=low_pose.phi(k+1)
		psi=low_pose.psi(k+1)
		if phi>180:
			phi-=360
		if psi>180:
			psi-=360
		phi_angles.append(phi)
		psi_angles.append(psi)
		
	print low_energy
#use matplotlib to make Ramachandran Plot

fig,ax=plt.subplots()
set ticks on x and y and axis
ax.set_xticks([-180,-90,0,90,180])
ax.set_yticks([-180,-90,0,90,180])
ax.xaxis.grid(True)
ax.yaxis.grid(True)
plt.axis([-180,180,-180,180])
label x and y axis
plt.xlabel('Phi Angles')
plt.ylabel('Psi Angles')
plot phi and psi angles
plt.scatter(phi_angles,psi_angles)
show plot
plt.show()
		
		
		