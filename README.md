# stamp-cell
<p align="center">
<img src="https://user-images.githubusercontent.com/33631502/159146148-1bbebdd1-e70e-4fd9-860f-9d075cfb5ae2.png" width=30% height=30%>
  
 
STAMP can be used for two different types of studies:

  1. **stamp_general**.py version provides different types of models for general membrane transporters (e.g., symporters slippage, symporters simplified model, antiporters simplified, etc.) discussed in [A comprehensive approach to the mathematical modeling of mass transport in biological systems: Fundamental concepts and models](https://www.sciencedirect.com/science/article/abs/pii/S0017931019365603) and
  2. **stamp_specific**.py version includes the models for the specific type of each membrane transporter and ion channels specified in [A mathematical modeling toolbox for ion channels and transporters across cell membranes](https://www.sciencedirect.com/science/article/abs/pii/S0017931021005263) (e.g., sodium chloride cotransporters, sodium-potassium ATPase pumps, etc.). Both general and specific stamps follow the same steps.
  
  
  ## Workflow
  To run these scripts, Python 3 must be installed. To begin, the user will be asked to choose a family of transporters to explore (e.g., symporters, antiporters, etc). Next, the user will be asked to select a subcategory of the chosen transporter (e.g., symporter, slippage, simplified, etc.). The user will then be prompted to enter the required data values for the parameters. TransporterDB provides users with a number of data files to choose from or they can enter their own data. Using the given data, Stamp will then run the simulation. Results from STAMP will be saved both in csv format and plots in the same directory as STAMP itself. The user will then be asked if any additional stimulation is desired. The process will be terminated if no response is received, and the process will be restarted if yes response is received. At each step, STAMP provides short keywords for the user to select from. STAMP recommends that users use the keywords in the exact way that STAMP suggests. The process will be stopped if an incorrect entry is made, and the message 'Invalid Answer. Exiting the program.' will display.
  
  While STAMP has not been fully completed, it contains moere than 40 models of human ion channels and membrane transporters. 

STAMP is developed by Shadi Zaheri and Fatemeh Hassanipour.
This work was supported by the National Science Foundation under grant number 1454334
  
  ## Copyright:
  If using the stamp_general pipeling,  please cite [A comprehensive approach to the mathematical modeling of mass transport in biological systems: Fundamental concepts and models](https://www.sciencedirect.com/science/article/abs/pii/S0017931019365603)
  
If using the stamp_specifc pipeling,  please cite [A mathematical modeling toolbox for ion channels and transporters across cell membranes](https://www.sciencedirect.com/science/article/abs/pii/S0017931021005263)
