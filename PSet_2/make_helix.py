from rosetta import *
rosetta.init()

def check_valid(num_resi):
'''Function takes in user's input as a string. Checks that input is an integer value or a 
float value with only zeroes after the decimal (ie '8' and '8.0' will be accepted, '8.4' will not).
Returns valid=True if input value is an integer. Returns valid=False if input value is not an integer.
Descriptions of error are given to user for invalid inputs. '''
	#initializes variable counting number of decimals to False
	decimal=False
	#iterate through every character in input
	for i in range(len(num_resi)):
		#char is the HTML number for the ASCII character
		char=ord(num_resi[i])
		
		#checks whether or not character is a decimal 
		if char==46 and decimal==False:
			#if character is decimal, marks that a decimal has occurred by changing decimal variable to 'True'
			decimal=True
		#if character is decimal, but another decimal has already occurred in the input sequence, cannot be valid
		elif char==46 and decimal==True:
			#alerts user to their error
			print ('ERROR: Invalid input. Please submit an integer value. Your input contained two decimal points')
			return False
		
		#checks whether character falls within HTML numbers for the 10 integers
		if char<48 or char>57:
			#only valid input that falls outside of this range is '46', the decimal point
			#ensures that this character isn't a decimal
			if decimal==False:
				print ('ERROR: Invalid input. One or more of the characters that you submitted were not integers. Please input only integers or whole number numbers (ex "8" and "8.0" ar both valid, but not "8.4")')
				return False
				
	#changes input string to type(float) so it can be divided my modulo division
	num_resi=float(num_resi)
	#checks that decimal number is a whole number
	if num_resi%1==0:
		return True
	else:
		print ('ERROR: Your value was not a whole number. Please submit only integer values')	
		return False


#initialize 'valid' variable at False
#User will be prompted continuously until they input a valid number or elect to stop inputting
valid=False
while valid==False:
	#prompt user for input
	num_resi=raw_input('What is the length of your helix? Please input an integer or whole number. ')
	#check that input is valid
	valid=check_valid(num_resi)

#Creates String of A's that is of user's specified length
A_str=''
for i in range(int(num_resi)):
	A_str+='A'

#creates helix of correct length called 'pose' 
pose=pose_from_sequence(A_str,'fa_standard')
#saves helix
pose.dump_pdb('helix.pdb')