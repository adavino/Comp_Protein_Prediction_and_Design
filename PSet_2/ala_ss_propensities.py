import sys
from rosetta import *
rosetta.init()

def determine_sec_struct(resi,resi_num,next_res,next_res_num,pose):
	'''Determines secondary structure that residue is likely involved in based on phi psi angles of residue and surrounding residues. Takes in residue number and type, residue number and type for next residue and the pose. Uses angles commonly found in different secondary structures,derived from Ramachandran plots. Different secondary structure criteria used for G and all other amino acids. If two consecutive amino acids ai and ai+1 have sec. structure, ai reported to have that structure. If they differ, ai reported to be in loop region. Returns a number that codes for secondary structure that the residue is likely involved in. Coded as follows:0=non_ala in helix, 1= non_ala in beta sheet, 2=ala in loop, 3=ala in beta sheet, 4= ala in in loop, 5= non_ala in loop'''
	
	#Retrieve torsion angles of amino acid i
	phi=pose.phi(resi_num)
	psi=pose.psi(resi_num)
	
	#retrieve torsion angles of next amino acid, i+1
	n_phi=pose.phi(next_res_num)
	n_psi=pose.psi(next_res_num)

	#angles for rama specific for glycine
	if resi=='G':
		if psi<60 and psi>-90 :
			if phi<-30 or phi>30:
				#G i in alpha helical region
				return 0
		

		elif psi <-90 or psi>90:
			if phi<-30 or phi>30:
				#G i in beta region
				return 1
		else:
			return 5
			



	#angles for rama generic for all other residues
	else:
		#resi i in alpha helical region
		if phi<0 and psi<60 and psi>-60:
			#resi in a helix
			if resi=='A':
				return 2
			else:
				return 0
		#resi i in beta sheet region
		elif phi<0 and psi>90:
			if resi=='A':
				return 3
				#resi in a beta sheet
			else:
				return 1
		else:
			if resi=='A':
				return 4
			else:
				return 5

def check_valid(pdb):
	'''Checks that input is possibly a valid PDB ID number. Takes in the user's input and checks that all of the characters are letter or numbers and that it is the correct length. '''
	#initializes valid as true, must have something wrong to return false
	valid=True
	for i in range(4):
		#checks that every character is letter or number
		if pdb[i].isalpha()==False and pdb[i].isdigit()==False:
			return False
	#if phrase (XXXX.pdb) isn't correct length, rejects it
	if len(pdb)!=8:
		return False
	#checks that file type is correct
	if pdb[4:]!='.pdb':
		return False
		
	return valid


	
def main():
	#initializes input as False, must be found to be valid to proceed
	valid=False
	while valid==False:
		#prompts user for input
		#accepts XXXX or XXXX.pdb
		pdb=raw_input('Please input PDB ID: ')
		#gives user the ability to leave the loop
		if pdb=='stop':
			print 'You have chosen to leave the program. Goodbye!'
			sys.exit()
		if len(pdb)==4:
			pdb+='.pdb'

		valid=check_valid(pdb)
		if valid==False:
			print "PDB ID was not valid. Please input only the 4 character PDB code"
			print "If you would like to stop inputting ID's, type in 'stop' when prompted for ID"
			
		
		
	
	from toolbox import cleanATOM, pose_from_rcsb
	
	cleanATOM(pdb)
	
	#create pose
	pose=pose_from_pdb(pdb[0:4]+'.clean'+pdb[4:])

	seq=pose.sequence()
	

	#initialize all variables that will count sec. struct. types for ala and non_ala resi

	non_ala_h=0
	non_ala_s=0
	non_ala_l=0
	ala_h=0
	ala_s=0
	ala_l=0

	for i in range(len(seq)-1):
		#find sec_struct of particular residue
		sec_struct=determine_sec_struct(seq[i],i+1,seq[i+1],i+2,pose)
		#classify residue structure, update appropriate variable
		if sec_struct==0:
			non_ala_h+=1
		elif sec_struct==1:
			non_ala_s+=1
		elif sec_struct==2:
			ala_h+=1
		elif sec_struct==3:
			ala_s+=1
		elif sec_struct==4:
			ala_l+=1
		elif sec_struct==5:
			non_ala_l+=1
	
	#divide ala totals by totals to get propensity
	p_ala_h=ala_h/float(non_ala_h+ala_h)
	p_ala_s=ala_s/float(non_ala_s+ala_s)
	p_ala_l=ala_l/float(non_ala_l+ala_l)

	print "Helix Propensity=",p_ala_h*100
	print "Sheet Propensity=",p_ala_s*100
	print "Loop Propensity=",p_ala_l*100

main()