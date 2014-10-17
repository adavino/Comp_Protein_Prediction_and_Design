from rosetta import *
rosetta.init()
from toolbox import cleanATOM



def calc_VDW(structure):
	'''Calculates full atom VDW using PyRosetta scoring function'''
	#attractive VDW
	VDW_atr=scorefxn.score_by_scoretype(structure,fa_atr)
	#repulsive VDW
	VDW_rep=scorefxn.score_by_scoretype(structure,fa_rep)
	VDW=VDW_atr+VDW_rep
	
	return VDW

def calc_Hbond(structure):
	''''Calculates hydrogen using PyRosetta scoring function'''
	E=0
	#H_bond types
	H_bond_lst=[hbond_sr_bb,hbond_lr_bb,hbond_bb_sc,hbond_sc]
	for i in range(len(H_bond_lst)):
		#gets Hbond, adds to total energy
		E+=scorefxn.score_by_scoretype(structure,H_bond_lst[i])
	return E

scorefxn=get_fa_scorefxn()

cleanATOM('1YY9.pdb')
cleanATOM('EGFR.pdb')
cleanATOM('centuximab.pdb')

#Gets complexes, makes them into pose
complex=pose_from_pdb('1YY9.clean.pdb')
EGFR=pose_from_pdb('EGFR.clean.pdb')
cent=pose_from_pdb('centuximab.clean.pdb')

structures=[complex,EGFR,cent]
#for labeling energies as they are printed
struct_names=['1YY9','EGFR','Centuximab']
def main():
	for j in range(len(structures)):
		#total pose energy FA score
		Energy=scorefxn(structures[j])
		#FA pose VDW
		E_VDW=calc_VDW(structures[j])
		#FA pose solvation energy
		sol=scorefxn.score_by_scoretype(structures[j],fa_sol)
		#FA pose Hbond energy
		E_Hbond=calc_Hbond(structures[j])
		
		print '=========='
		print struct_names[j]
		print "Van der Waals Energy=",E_VDW
		print "Solavation Energy=",sol
		print "Hydogen Bonding Energy=",E_Hbond
		print "Total Energy=",Energy
		print '+++++++'
main()