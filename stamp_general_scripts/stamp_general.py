import os
import shutil
import re
from PIL import Image, ImageFilter

cwd = os.getcwd()
print(cwd)
###### If Possible Add a schematic of all then ask which one, like png vwesion of just schematic configuration.
model = input('Enter Cotransporter Type:\n symporter = sym \n antiporter = anti \n pumps = pump \n uniporter = uni \n')
if(str(model) != ('anti'and 'sym')): # Add uniporter and Pump
	print('Invalid Answer. Exiting the program.')
if (str(model) == 'sym'):
	im = Image.open( 'Symporter.jpg' )
	im.show()
	outmodel1 = input('Choose symporter model type:\n slippage = Slp \n TwoStatesAFirst = AFirst \n TwoStatesBFirst = BFirst \n ABAC_ordered = ABAC \n ABCB_ordered = ABCB \n TwoSolutesAndCompetitor = Compete \n SimplifiedModel = Simple \n' )
	if(str(outmodel1) == 'Slp' or str(outmodel1) == 'AFirst' or str(outmodel1) == 'BFirst' or str(outmodel1) == 'ABAC' or str(outmodel1) == 'ABCB' or str(outmodel1) == 'Compete' or str(outmodel1) == 'Simple'):
		print('The Symporter ' + str(outmodel1) + ' model is chosen')
		os.system('python3 ' + 'Sym' + str(outmodel1) + '.py')
	else:
		print('Invalid Answer. Exiting the program.')
if (str(model) == 'anti'):
	im = Image.open( 'Antiporter.jpg' )
	im.show()
	outmodel2 = input('Choose antiporter model type:\n EightState = EightState \n SixState = SixState \n TwoSolutesAndCompetitor= Compete \n SimplifiedModel = Simple \n' )
	if(str(outmodel2) == 'Slp' or str(outmodel2) == 'EightState' or str(outmodel2) == 'SixState' or str(outmodel2) == 'Compete' or str(outmodel2) == 'Simple'):
		print('The Antiporter ' + str(outmodel2) + ' model is chosen')
		os.system('python3 ' + 'Anti' + str(outmodel2) + '.py')
	else:		
		print('Invalid Answer. Exiting the program.')
if (str(model) == 'pump'):
	outmodel3 = input('Choose pump model:\n One_Substrate_Detailed = one \n OneSubstratesHills = hill \n One_Substrate_Detailed = one \n TwoSubstrates = two \n ' )
	if(str(outmodel3) == 'one' or str(outmodel3) == 'hill' or str(outmodel3) == 'two'):
		print('The Pump ' + str(outmodel3) + ' model is chosen')
		os.system('python3 ' + 'pump' + str(outmodel3) + '.py')
	else:		
		print('Invalid Answer. Exiting the program.')

if (str(model) == 'uni'):
	outmodel4 = input('Choose uniporter model:\n Detailed = det \n MichaelisMenten = MM \n ' )
	if(str(outmodel4) == 'det' or str(outmodel4) == 'MM'):
		print('The Uniporter ' + str(outmodel4) + ' model is chosen')
		os.system('python3 ' + 'uni' + str(outmodel4) + '.py')
	else:		
		print('Invalid Answer. Exiting the program.')

