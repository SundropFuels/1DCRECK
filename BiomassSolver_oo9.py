'''
Created on May 31, 2011

@author: Chris Perkins & Jordan Jalving
'''
#!/usr/bin/python2.7


from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline
from scipy.optimize import newton
import numpy as np

from Cantera import *
from Cantera.constants import StefanBoltz
from Cantera.constants import GasConstant as R

import sys
import csv
import time



class Tube:
    """
        Creates a tube object that can be used in class LSODA simulator.  Class tube contains data members that define tube diameter,
        length,emissivity,wall temperature,wall flux,etc..., and calculates cross sectional area.
        Additionally, the Overall Heat Transfer Coefficient may be set, which only applies to Gas Phase Models.
    """
    def __init__(self, diameter = 3.5, length = 7.0, emissivity = 1.0, temperature = 1400):
        self.diameter = diameter*2.54/100
        self.length = length*0.3048
        self.emissivity = emissivity
        self.temperature = temperature + 273.15
        self.XSA = 3.14159*(self.diameter/2.0)**2.0
        self.U = 100   #Overall Heat Transfer Coefficient W/m^2 - K

    def setTube(self, diameter = 3.5, length = 7.0, emissivity = 1.0, temperature = 1400):
        self.diameter = diameter*2.54/100.0
        self.length = length*0.3048
        self.emissivity = emissivity
        self.temperature = temperature + 273.15
        self.XSA = 3.14159*(self.diameter/2.0)**2.0
        
    def setDiameter(self,diameter = 3.5):
        self.diameter = diameter*2.54/100.0
        self.XSA = 3.14159*(self.diameter/2.0)**2.0
    def setLength(self,length = 7.0):
        self.length = length*0.3048
    def setEmissivity(self,emissivity = 1.0):
        self.emissivity = emissivity
    def setWallTemperature(self,temperature = 1450):
        self.temperature = temperature + 273.15
    def setHeatTransferCoefficient(self,U = 150):
        self.U = U
    def setWallFlux(self,heatflux = 60):
        self.heatflux = heatflux

#MODULE LEVEL FUNCTION.  
def clearPhase():
    """
        This is a module-level function that clears class attributes in class Phase.
        This could be implemented as a class method, but the module level function makes enough sense.
        This method is called in the Runner module to clear this list so simulations can be run multiple times from the 
        GUI
    """
    Phase.phaseCounter = 0
    Phase.phaseList = []
    Phase.phaseDict = {}

class Phase:
    """
        Creates a Phase object that can be used in subclasses SolidGas and Gas.  The Phase object acts
        as a wrapper for a cantera phase such as a gas object.  Additionally, it stores the inlet data of
        a phase such as flow rates of each species and temperature.  It also holds the filename for the 
        reaction mechanism and a few other helpful attributes.  These values are set with the setInletConditions
        method in Class LSODASimulator
    """
    phaseCounter = 0
    phaseList = []
    phaseDictionary = {}
    def __init__(self,**kwargs):
        """
            kwargs: name
        """
        Phase.phaseCounter +=1
        self.name = kwargs.get('name','Phase %s '%(Phase.phaseCounter))  #Either set name to the input, or set to 'Phase' with a counter
        Phase.phaseDictionary[self.name]=(self)
        Phase.phaseList.append(self)
        self.filename = None                #Filename for reaction mechanism of phase
        self.cantera_phase = None           #The cantera phase object
        self.m_0 = None                     #Dictionary of inlet flow rates for the phase
        self.phase_string = None            #A string of the inlet flow rates that can be used to set the inlet composition
        self.TempIn = None                  #The inlet temperature of the phase
        self.pressure = None                #The inlet pressure of the phase (Note: Only 1 pressure can be entered for all phases)
        self.total_flow = 0.0               #[kg/s] Inlet mass flow rate of phase
        self.total_molar_flow = 0.0         #[kmol/s] Inlet molar flow rate of phase 
        self.total_molar_flow_out = 0.0     #[kmol/s] Outlet molar flow rate of a phase.  This is used in the heat duty calculation.
        self.species_indices = None         #Index of the species in the cantera phase 
        self.species_names = None           #Names of the species in the cantera phase
        self.moleFractions = np.zeros(1)    #Mole fractions of the species in the cantera phase (Where do I use this?)
        self.temperatureIndex = 0           #Very Important.  This is the temperature index of the phase.
                                            #The solver determines which phase temperature its referencing based on this index.
        #self.temperatureProfiles = None    #No longer using this
        
    def printInfo(self):
        """
            Print phase information.
        """
        print "name: %s "%self.name
        print "m_0: %s "%self.m_0
        print "TempIn %s "%self.TempIn
        print "pressure %s "%self.pressure
        

        
#class ExperimentalHandler:
#    """
#    This class will act as a handler for reading experimental data.  Before coding this, a reasonable 
#    input format needs to be decided.  (csv columns, names of fields, etc....)  It is possible that this 
#    class may mesh better with the plotting toolbox.
#    """
#    def __init__(self):
#        #CSV File Data
#        self.steamTempData = []
#        self.pressureData = []      
#        self.biomassFeedRateData = []
#        self.steamFeedRateData = []
#        self.exitTemperatureData = []
#        self.exitConversionData = []
#        self.wallTempData = []
#        self.exitH2Data = []
#        self.exitCH4Data = []
#        self.exitCOData = []
#        self.exitCO2Data = []
#        
#    def readExperimentalData(self,inputfile):
#        reader = csv.reader(open(inputfile,'r'),delimiter = ';',quoting = csv.QUOTE_NONNUMERIC)
#        #reader = csv.reader(open(inputfile,'r'),delimiter = ',',quoting = csv.QUOTE_NONNUMERIC)
#        reader.next()
#        for row in reader:
#            
#            self.biomassFeedRateData.append(row[9])
#            self.steamFeedRateData.append(row[33])
#            self.exitTemperatureData.append(row[4])  #Celsius
#            self.exitConversionData.append((float(row[7])*100.0/self.maxConversion))
#            self.exitH2Data.append(row[18])
#            self.exitCH4Data.append(row[34])
#            self.exitCOData.append(row[32])
#            self.exitCO2Data.append(row[12])
#            self.steamTempData.append(row[24])
#            self.pressureData.append(row[14])
#            self.wallTempData.append((((row[36]+273.15)**4 + (row[37]+273.15)**4 + (row[38]+273.15)**4 + (row[39]+273.15)**4 + (row[40]+273.15)**4)/5)**0.25-273.15) #Celsius
#            #self.wallFluxData = []
#            #self.effectivenessData = []
#    
#            #Data from Sunulator 1
#            #self.biomassFeedRateData.append(row[1])
#            #self.steamFeedRateData.append(row[2])
#            #self.exitConversionData.append((float(row[3])*100.0/self.maxConversion))
#            #self.exitCH4Data.append(row[4])
#        inletData = {'biomassFeedRateData':self.biomassFeedRateData,'steamTempData':self.steamTempData,'pressureData':self.pressureData,'wallTempData':self.wallTempData,'steamFeedRateData':self.steamFeedRateData}
#        #inletData = {'biomassFeedRateData':self.biomassFeedRateData,'steamFeedRateData':self.steamFeedRateData}
#        return inletData
            
class LSODASimulator:
    """
        Creates a simulator object that contains methods to setup and solve a mass and energy balance for a tube gasifier.
        The initial inputs for the class include a tube object instantiated from class Tube.
    """

    def __init__(self, tube):

        #Create an initial tube 
        self.tube = tube                      #The simulator must take an instance of Class Tube
        self.normalizedSpecies = []           #A list of species to normalize concentration data against
        self.rtol = 1E-12                     #Integration relative tolerance
        self.atol = 1E-12                     #Integration absolute tolerance
        self.steps = 1000                     #Number of steps the solver will return values for
        
        #Initiate the simulation count.  The coupler class will increment these values.
        self.simulationCount = 0              #Biomass Simulation Counter            
        self.MCCount = 0                      #Monte-Carlo Simulation Counter

             
        self.maxConversion = 87.0               #PROBABLY NEED TO ADD THIS TO THE GUI
        self.gasConversion = []                 #THIS DEPENDS ON SPECIES OF INTEREST.  I WILL NEED TO FIND A WAY TO IMPLEMENT THIS BETTER.

    def setInletConditions(self, phaseDict):
        """
           Sets up system initial conditions.  phaseDict is a dictionary passed by a subclass containing
           a Phase instance and parameters to set the phase.  The method loops through every phase in the
           dictionary and sets their parameters.  It also generates several data members for use throughout
           the simulation.
        """
        #Setup lists to grab relevant phase data.  The solver inlet condition is set up using these lists
        phaseMoleFractions = []
        phaseTemperatures = []
        phaseTotalMolarFlow = []
        phaseSetCount = 0
        for phase,parameters in phaseDict.items():                   #Loop through each phase in phaseDict
            if phase.filename != parameters['phaseMechanism']:       #Determine if the reaction mechanism changed
                phase.cantera_phase = importPhase(parameters['phaseMechanism']) #Create the cantera phase if its a new mechanism
                phase.filename = parameters['phaseMechanism']
            phase.m_0 = parameters['m_0']                            #Dictionary of mass flow rates for each phase species
            #-----------------------------------------------------------------------------------
            #-----------------------------------------------------------------------------------
            #BUILD STRING OF SPECIES FROM DICTIONARY IN THE FORM 'SPECIES:AMOUNT'
            phase.phase_string = ""
            for k,v in phase.m_0.items():
                phase.phase_string += "%s:%s," % (k,v)
            phase.phase_string = phase.phase_string[0:-1]            #Remove last comma from string
            #-----------------------------------------------------------------------------------
            #-----------------------------------------------------------------------------------
            #SET GAS INLET CONDITIONS (TEMPERATURES, PRESSURE, COMPOSITION)
            phase.TempIn = parameters['TempIn'] + 273.15                 #self.TempIn is in K
            phase.pressure = (parameters['pressure']+14.7)*101325.0/14.7 #self.pressure is in Pa
            phase.cantera_phase.set(T = phase.TempIn, P = phase.pressure, Y = phase.phase_string)
            #-----------------------------------------------------------------------------------
            #-----------------------------------------------------------------------------------  
            #CALCULATE TOTAL INLET FLOW RATES FOR THE PHASE
            phase.total_flow = 0.0
            for value in phase.m_0.values():
                phase.total_flow += value/2.204/3600.0                         #[kg/s] Converts lb/hr of each species to kg/s and sums
            phase.total_molar_flow = phase.total_flow/phase.cantera_phase.meanMolarMass() #[kmol/s]convert to kmol/s
            #CALCULATE AND GRAB DATA NEEDED FOR SETTING INLET CONDITION
            phase.moleFractions = phase.cantera_phase.moleFractions()          
            phaseMoleFractions.append(phase.moleFractions)
            phaseTemperatures.append(phase.TempIn)                             
            phaseTotalMolarFlow.append(phase.total_molar_flow)
            #MAP SPECIES INDICES INTO A DICTIONARY     
            phase.species_indices = {}
            phase.species_names = phase.cantera_phase.speciesNames()
            for name in phase.species_names:
                phase.species_indices[name] = phase.species_names.index(name)
            #SET THE PHASE'S TEMPERATURE INDEX.  THIS IS NEEDED SINCE I CAN'T ENFORCE AN ORDER IN WHICH PHASES ARE SET
            phase.temperatureIndex = phaseSetCount
            phaseSetCount +=1
 
        #DATE MEMBERS USED THROUGHOUT THE PROGRAM
        self.species_indices = phase.species_indices  #This should give all species and their index since every phase should have the same species
        self.cantera_species = phase.cantera_phase    #This should give a cantera phase object corresponding to the last phase entered.  It should be used for referencing species data.
        #create initial value matrix for flow
        self.n_0 = np.sum([phaseMoleFractions[phase]*phaseTotalMolarFlow[phase] for phase in range(Phase.phaseCounter)],axis = 0)
        #create complete matrix with phase temperatures appended to the end
        self.y0 = np.append(self.n_0,phaseTemperatures)
        self.tempIndex = len(self.y0)-len(self.n_0)   #Index in y0 where temperature data begins
    
    #THIS METHOD NEEDS TO 
    def calc_carbon_conversion(self,productSpecies):
        """
            Calculate the carbon conversion based on carbon turned to a product species.
            This method will only work once method solveLSODA has been called.
        """
        try:
            n_C_inlet = 0.0
            n_C_outlet = 0.0
            for species, index in self.species_indices.items():
                n_C_inlet += self.n_0[index]*self.cantera_species.nAtoms(species = species, element = 'C')
                if species in productSpecies:     #productSpecies is defaulting to 'CO,CO2, and CH4'
                    n_C_outlet += self.n[:,index]*self.cantera_species.nAtoms(species = species, element = 'C')
            return (n_C_outlet/n_C_inlet)*100.0*(100/self.maxConversion)   #max conversion is defaulting to 87.0
        except AttributeError:
            print 'No outlet data available.  Run solve first.'


    def systemEnthalpy(self,phaseFlows):
        """
            This method calculates the system enthalpy of each Phase instance and sums them to get a 
            total.  It takes the dictionary argument phaseFlows, which is a total molar flow rate of each
            phase.
        """
        
        phaseEnthalpies = []
        for phase, molarFlow in phaseFlows.items():
            phaseEnthalpies.append(sum(phase.cantera_phase.moleFractions()*molarFlow*
                                       phase.cantera_phase.enthalpies_RT()*R*phase.cantera_phase.temperature()))
        return sum(phaseEnthalpies)    

    def calc_net_energy(self):
        """
            Calculate heat duty of gasifier.  This method sets the Phase objects to the inlet and outlet
            conditions and then calculates the enthalpy at each condition.
        """
        try:
            phaseFlows = {}
            for phase in Phase.phaseList:
            #need to set the gas at each temp, pressure here to be accurate
                phase.cantera_phase.set(T = phase.TempIn,P = phase.pressure, Y = phase.phase_string)
                phaseFlows[phase] = phase.total_molar_flow                
            H_inlet = self.systemEnthalpy(phaseFlows)

            phaseFlows = {}
            for phase in Phase.phaseList:
                phase.cantera_phase.set(T = self.phaseTemps[-1,phase.temperatureIndex],P = phase.pressure, Y = phase.total_molar_flow_out)
                phaseFlows[phase] = sum(phase.total_molar_flow_out)
            H_outlet = self.systemEnthalpy(phaseFlows)

            return (H_outlet - H_inlet)
        
        except AttributeError:
            print 'No outlet data available.  Run solve first.'

    def calcMixingEnergy(self, Tmix,gas,flows,Temps):
        """
            Calculate the Mixing enthalpy of two gases.  Tmix is a guess of the mixture temperature, gas is a cantera phase object,
            flows is the mass flow rate of each species, and Temps is the initial temperature of each species.
            A root solver can be applied to this method to solve for the actual mixing temperature
        """
        inletEnergy = []
        mixingEnergy = []
        for key,value in flows.items():
            gas.set(T = Temps[key],Y = '%s:%s'%(key,value))
            m = value/2.204/3600.0
            inletEnergy.append(m*gas.enthalpy_mass())
            gas.set(T = Tmix)
            mixingEnergy.append(m*gas.enthalpy_mass())
        return sum(inletEnergy) - sum(mixingEnergy)
        

    def setTolerance(self, rtol = 1E-12, atol = 1E-12):
        """
            Set the absolute and relative tolerance for integration
        """
        self.rtol = rtol
        self.atol = atol
    
    def setSteps(self, steps = 10000):
        """
        Set the number of steps to return integration solutions at.  This is not the number of steps used to solve mass and energy balance.
        A smaller number of steps however, can increase run speed
        """
        self.steps = steps

    def solveLSODA(self,*args):
        """
        Solve the differential equation using LSODA from ODEPack.  This is implemented through scipy odeint.
        Any args given to this function will be passed to the odeint call, and subsequently, its mass and energy balance.
        """
        #DEFINE TUBE DISTANCES TO RETURN DATA AT AND RUN THE SIMULATION
        self.z = np.linspace(0,self.tube.length, self.steps)
        self.wallflux = [] #List for storing wall flux values from the mass and energy balance
        self.wallTemp = [] #List for storing wall temperature values for the mass and energy balance
        self.length = []   #This list stores all of the length values used in the mass and energy balance, not the stripped list
        
        print "Simulation Count %s"% self.simulationCount
        print "Solving the differential equation..."
        y = odeint(self.massAndEnergyBalance, self.y0, self.z, args, rtol = self.rtol, atol = self.atol,mxstep = 2000)
        print "ME Balance solution complete"
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #CUT UP THE SOLUTION MATRIX INTO RELEVANT DATA
        #SOLID TEMPERATURE, GAS TEMPERATURE, SOLID MOLAR FLOW, GAS MOLAR FLOW
        self.phaseTemps = y[:,-self.tempIndex:]
        self.n = y[:,0:-self.tempIndex]            
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #CALCULATE CARBON CONVERSION AND HEAT DUTY
        self.X = self.calc_carbon_conversion(['CH4','CO','CO2'])
        self.Q = self.calc_net_energy()

        #REDUCE THE SIZE OF THE WALLFLUX AND WALLTEMP LISTS TO THE SAME SIZE RETURNED FROM ODEINT
        getSampleSet = lambda x,y,x2: list(UnivariateSpline(x,y)(x2))
        self.wallflux = getSampleSet(self.length,self.wallflux,self.z)  
        self.wallTemp = getSampleSet(self.length,self.wallTemp,self.z)
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #CALCULATE THE COMPOSITION AND MOLAR FLOW OF EACH SPECIES FOR EVERY DISTANCE POINT 
        #THIS ALSO CALCULATES THE COMPOSITION NORMALIZED
        self.mol_frac = {}
        self.mol_fracNormalized = {}
        self.molarFlow = {}
        for species, index in self.species_indices.items():
            self.mol_frac[species] = (self.n[:,index]/self.n.sum(axis=1))    #column of species molar flow rates divided by column of sums of each row
        for species, composition in self.mol_frac.items():
            self.mol_fracNormalized[species] = composition/(1 - sum(self.mol_frac[species] for species in self.normalizedSpecies))
        for species, index in self.species_indices.items():
            self.molarFlow[species] = self.n[:,index]

    
    def write_data(self,filename = "/home/jordan/code/python/simulationResults.csv"):
        """
            Write simulation data to a csv file.  The optional input is a filepath to a .csv file
            to write the data into.  This method should be run after all simulations are finished.
            This write method will only write run conditions and exit values
        """
        #WRITE THE CSV FILE HEADER DURING THE FIRST SIMULATION
        if self.simulationCount == 1: 
            inputfile = open(filename , 'w')
            headerwriter = csv.writer(inputfile,delimiter = ',',quoting=csv.QUOTE_MINIMAL)
            phaseTemperatures = ['Exit ' + phase.name + ' Temperature [C]' for phase in Phase.phaseList]
            phaseTotalFlows = [phase.name + 'Total Flow [lb/hr]' for phase in Phase.phaseList]
            speciesConcentrations = [name + 'Concentration' for name in self.cantera_species.speciesNames()]
            speciesConcentrations2 = []
            for name in speciesConcentrations:
                if name.strip('Concentration') in self.normalizedSpecies:
                    name = 'N- ' + name
                speciesConcentrations2.append(name)
            speciesFlows = [name + 'mol/min' for name in self.cantera_species.speciesNames()]
            headerwriter.writerow((np.hstack(('RunID','Length [ft]','Pressure [psig]','Exit Carbon Conversion','Heat Duty [W]',phaseTotalFlows,phaseTemperatures,speciesConcentrations2,speciesFlows))))
        else: 
            inputfile = open(filename,'a')

        datawriter = csv.writer(inputfile,delimiter = ',',quoting = csv.QUOTE_NONNUMERIC)
        exitTemps = [self.phaseTemps[-1,phase.temperatureIndex] - 273.15 for phase in Phase.phaseList]
        exitCompositions = [self.mol_fracNormalized[species][-1] for species in self.cantera_species.speciesNames()]
        exitFlowRates = [self.molarFlow[species][-1]*1000.0*60.0 for species in self.cantera_species.speciesNames()] #mol/min
        totalMassFlowRates = [phase.total_flow*2.204*3600.0 for phase in Phase.phaseList]
        datawriter.writerow((np.hstack((self.simulationCount,self.z[-1],Phase.phaseList[0].pressure*14.7/101325-14.7,self.X[-1],self.Q,totalMassFlowRates,exitTemps,exitCompositions,exitFlowRates))))
        inputfile.close()
    
    def writeProfileData(self,filename):
        """
            Write the profile data for each simulation to a csv file.  This will work like the write method 
            I used for Tim.  The PlottingToolBox will chop up the profiles depending on simulation count.
        """
        #WRITE THE CSV FILE HEADER DURING THE FIRST SIMULATION
        if self.simulationCount == 1: 
            inputfile = open(filename , 'w')
            headerwriter = csv.writer(inputfile,delimiter = ',',quoting=csv.QUOTE_MINIMAL)
            phaseTemperatures = [phase.name + ' Temperature [C]' for phase in Phase.phaseList]
            speciesConcentrations = [name + 'Concentration' for name in self.cantera_species.speciesNames()]
            speciesConcentrations2 = []
            #I could use a lambda function to make this into a list comprehension
            #PUT A 'N- ' IN FRONT OF NORMALIZED SPECIES IN THE CSV FILE 
            for name in speciesConcentrations: 
                if name.strip('Concentration') in self.normalizedSpecies:
                    name = 'N- ' + name
                speciesConcentrations2.append(name)
            speciesFlows = [name + 'mol/min' for name in self.cantera_species.speciesNames()]
            headerwriter.writerow((np.hstack(('RunID','Distance','Carbon Conversion','Wall Flux [kW/m^2]',
                                              'Wall Temperature [C]',phaseTemperatures,speciesConcentrations2,
                                              speciesFlows))))
        else:
            inputfile = open(filename, 'a')
        datawriter = csv.writer(inputfile,delimiter = ',',quoting=csv.QUOTE_NONNUMERIC)  
        #NOTE: Converted Array will be inverted.
        phaseTempProfiles = np.asarray([self.phaseTemps[:,phase.temperatureIndex] - 273.15 for phase in Phase.phaseList])
        compositionProfiles = np.asarray([self.mol_fracNormalized[species] for species in self.cantera_species.speciesNames()])
        molarRateProfiles = np.asarray([self.molarFlow[species] for species in self.cantera_species.speciesNames()])
        for count in range(len(self.z)):
            datawriter.writerow((np.hstack((self.simulationCount,self.z[count],self.X[count],
                                self.wallflux[count],self.wallTemp[count],phaseTempProfiles[:,count],
                                compositionProfiles[:,count],molarRateProfiles[:,count]))))
        inputfile.close()
         
    #THIS METHOD SHOULD BE CHANGED TO DETERMINE GAS CONVERSION BASED ON INLET AND OUTLET
    def calc_gas_conversion(self,selectSpecies):
        """
            Calculate the conversion of the reactant species.  
        """
        #Data for Gas Cracking plots
        gasMolsIn = self.n_0[self.species_indices[selectSpecies]]#*1000.0*60.0
        otherMolsIn = 0.0
        for species, index in self.species_indices.items():
            if selectSpecies != species:
                otherMolsIn += self.n_0[index]
        maxPercentMolsOut = gasMolsIn/(gasMolsIn + otherMolsIn)
        gasMolsOut = self.exitCompositions[selectSpecies][-1]/maxPercentMolsOut*gasMolsIn
        self.gasConversion.append(1 - (gasMolsOut/gasMolsIn))
        print "Gas Molar Conversion = %s"%(1 - (gasMolsOut/gasMolsIn))

   
 


#    def plot_gasModel(self,showplots = 'yes',selectSpecies = 'CH4'):
#        figure(101)
#        #data = [0.58,0.49,0.52,0.51,0.39,0.32,0.31,0.27,0.20,0.09,0.04,0.03,0.02,0.01,0.005]
#        #scatter(self.wallTemp,data, marker = 'o', s = 30, color = 'red')
#        #plot(self.wallTemp,self.exitCompositions[selectSpecies],color = 'blue')
#        NUM_COLORS = 15
#        cm = get_cmap('Set1')
#        colorIndex = 0
#        for species, compositions in self.exitCompositions.items():
#            if self.exitCompositions[species][-1] > 0.001:
#                color = cm(1.*colorIndex/NUM_COLORS)  # color will now be an RGBA tuple
#                plot(self.wallTemp,self.exitCompositions[species],color = color, label = species,linewidth = 2.0)
#                colorIndex += 1
#        suptitle('Gas Fraction versus Wall Temperature',fontsize = 24)
#        xlabel('Wall Temperature [C]',fontsize = 20)
#        ylabel('Normalized vol %s'% (selectSpecies,),fontsize = 20)
#        legend()
#
#        figure(102)
#        #data = [0.31,0.41,0.38,0.39,0.52,0.61,0.62,0.68,0.77,0.89,0.93,0.96,0.98,0.99,0.98]
#        plot(self.wallTemp,self.gasConversion,color = 'blue')
#        #scatter(self.wallTemp,data,marker = 'o', s= 30,color = 'red')
#        suptitle('Gas Cracking Conversion versus Temperature',fontsize = 24)
#        xlabel('Wall Temperature [C]',fontsize = 20)
#        ylabel('Cracking Conversion %s'% (selectSpecies,),fontsize = 20)
#
#        if showplots == 'yes':
#            show()
#
#    def plot3D(self,X,Y,figureDirectory = None, showplots = 'yes'):
#        X,Y = np.meshgrid(X,Y)
#        hydrogenslpm = round(self.n_0[self.species_indices['H2']]*1000*60*24,3)
#        TempIn = round(self.y0[-1]-273.15) #C
#        diameter = round(self.tube.diameter/2.54*100.0,3)
#        
#        CgrOut = np.asarray([matrix[:,self.species_indices['C(gr)']][-1] for matrix in self.molarFlowsOut])
#        CH4In = np.asarray([array[self.species_indices['CH4']] for array in self.molarFlowsIn])
#        CH4Out = np.asarray([matrix[:,self.species_indices['CH4']][-1] for matrix in self.molarFlowsOut]) 
#        H2In = np.asarray([array[self.species_indices['H2']] for array in self.molarFlowsIn])
#        H2Out = np.asarray([matrix[:,self.species_indices['H2']][-1] for matrix in self.molarFlowsOut])
#        TOut = np.asarray([array[-1] for array in self.T_phases])-273.15 #C
#        CO2Out = np.asarray([matrix[:,self.species_indices['CO2']][-1] for matrix in self.molarFlowsOut])
#        COOut = np.asarray([matrix[:,self.species_indices['CO']][-1] for matrix in self.molarFlowsOut])
#
#        #HYDROGEN RATIO PLOT
#        fig = figure(100 + self.simulationCount,figsize = (20,10))
#        ax = Axes3D(fig)
#        Z = (H2Out/H2In).reshape(len(X),len(Y))
#        surf = ax.plot_surface(X, Y, Z, cstride = 0.25,rstride = 0.25, cmap=get_cmap('jet'),linewidth=0, antialiased=False)
#        xLabel = ax.set_xlabel('Methane [slpm]', linespacing=3.2,fontsize = 16)
#        yLabel = ax.set_ylabel('Oxygen [slpm]', linespacing=3.1,fontsize = 16)
#        zLabel = ax.set_zlabel('H2 Out/H2 In', linespacing=3.4,fontsize = 16)
#        title = fig.suptitle('Hydrogen Production [H2 In = %s slpm, T In = %s C, Diameter = %s in]'%(hydrogenslpm,TempIn,diameter),fontsize = 16)
#        ax.set_zlim3d(0.01, amax(Z)*1.05)
#        ax.w_zaxis.set_major_locator(LinearLocator(10))
#        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        fig.savefig(figureDirectory + 'HydrogenProductionPlot_H2_%s_slpm.png'%hydrogenslpm)
#        
#        #CARBON OUT TO METHANE IN PLOT
#        fig = figure(200 + self.simulationCount,figsize = (20,10))
#        ax = Axes3D(fig)
#        Z = (CgrOut/CH4In).reshape(len(X),len(Y))*100
#        surf = ax.plot_surface(X, Y, Z, cstride = 0.25,rstride = 0.25, cmap=get_cmap('jet'),linewidth=0, antialiased=False)
#        xLabel = ax.set_xlabel('Methane [slpm]', linespacing=3.2,fontsize = 16)
#        yLabel = ax.set_ylabel('Oxygen [slpm]', linespacing=3.1,fontsize = 16)
#        zLabel = ax.set_zlabel('C Out/CH4 In', linespacing=3.4,fontsize = 16)
#        title = fig.suptitle('Ratio of Carbon Out to Methane In [H2 In = %s slpm, T In = %s C, Diameter = %s in]'%(hydrogenslpm,TempIn,diameter),fontsize = 16)
#        ax.set_zlim3d(0.00, 100.0)
#        ax.w_zaxis.set_major_locator(LinearLocator(10))
#        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        fig.savefig(figureDirectory + 'CarbonToMethane_H2_%s_slpm.png'%hydrogenslpm)
#        
#        #METHANE CONVERSION PLOT
#        fig = figure(300 + self.simulationCount,figsize = (20,10))
#        ax = Axes3D(fig)
#        Z = (1.0 - (CH4Out/CH4In)).reshape(len(X),len(Y))*100.0
#        surf = ax.plot_surface(X, Y, Z, cstride = 0.25,rstride = 0.25, cmap=get_cmap('jet'),linewidth=0, antialiased=False)
#        xLabel = ax.set_xlabel('Methane [slpm]', linespacing=3.2,fontsize = 16)
#        yLabel = ax.set_ylabel('Oxygen [slpm]', linespacing=3.1,fontsize = 16)
#        zLabel = ax.set_zlabel('Methane Conversion %', linespacing=3.4,fontsize = 16)
#        title = fig.suptitle('Methane Conversion [H2 In = %s slpm, T In = %s C, Diameter = %s in]'%(hydrogenslpm,TempIn,diameter),fontsize = 16)
#        ax.set_zlim3d(0.00, 100.00)
#        ax.w_zaxis.set_major_locator(LinearLocator(10))
#        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        fig.savefig(figureDirectory + 'MethaneConversion_H2_%s_slpm.png'%hydrogenslpm)
#        
#        #EXIT TEMPERATURE PLOT
#        fig = figure(400 + self.simulationCount,figsize = (20,10))
#        ax = Axes3D(fig)
#        Z = TOut.reshape(len(X),len(Y))
#        surf = ax.plot_surface(X, Y, Z, cstride = 0.25,rstride = 0.25, cmap=get_cmap('jet'),linewidth=0, antialiased=False)
#        xLabel = ax.set_xlabel('Methane [slpm]', linespacing=3.2,fontsize = 16)
#        yLabel = ax.set_ylabel('Oxygen [slpm]', linespacing=3.1,fontsize = 16)
#        zLabel = ax.set_zlabel('Exit Temperature [C]', linespacing=3.4,fontsize = 16)
#        title = fig.suptitle('Exit Temperature [H2 In = %s slpm, T In = %s C, Diameter = %s in]'%(hydrogenslpm,TempIn,diameter),fontsize = 16)
#        ax.set_zlim3d(amin(Z)*0.95, amax(Z)*1.05)
#        ax.w_zaxis.set_major_locator(LinearLocator(10))
#        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
#        fig.colorbar(surf, shrink=0.5, aspect=5.0)
#        fig.savefig(figureDirectory + 'ExitTemperature_H2_%s_slpm.png'%hydrogenslpm)
#        
#        #CO2 OUT TO METHANE IN PLOT
#        fig = figure(500 + self.simulationCount,figsize = (20,10))
#        ax = Axes3D(fig)
#        Z = (CO2Out/CH4In).reshape(len(X),len(Y))*100
#        surf = ax.plot_surface(X, Y, Z, cstride = 0.25,rstride = 0.25, cmap=get_cmap('jet'),linewidth=0, antialiased=False)
#        xLabel = ax.set_xlabel('Methane [slpm]', linespacing=3.2,fontsize = 16)
#        yLabel = ax.set_ylabel('Oxygen [slpm]', linespacing=3.1,fontsize = 16)
#        zLabel = ax.set_zlabel('CO2 Out/CH4 In', linespacing=3.4,fontsize = 16)
#        title = fig.suptitle('Ratio of CO2 Out to Methane In [H2 In = %s slpm, T In = %s C, Diameter = %s in]'%(hydrogenslpm,TempIn,diameter),fontsize = 16)
#        ax.set_zlim3d(0.00, 100.0)
#        ax.w_zaxis.set_major_locator(LinearLocator(10))
#        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        fig.savefig(figureDirectory + 'CO2ToMethane_H2_%s_slpm.png'%hydrogenslpm)
#        
#        #CO OUT TO METHANE IN PLOT
#        fig = figure(600 + self.simulationCount,figsize = (20,10))
#        ax = Axes3D(fig)
#        Z = (COOut/CH4In).reshape(len(X),len(Y))*100
#        surf = ax.plot_surface(X, Y, Z, cstride = 0.25,rstride = 0.25, cmap=get_cmap('jet'),linewidth=0, antialiased=False)
#        xLabel = ax.set_xlabel('Methane [slpm]', linespacing=3.2,fontsize = 16)
#        yLabel = ax.set_ylabel('Oxygen [slpm]', linespacing=3.1,fontsize = 16)
#        zLabel = ax.set_zlabel('CO Out/CH4 In', linespacing=3.4,fontsize = 16)
#        title = fig.suptitle('Ratio of CO Out to Methane In [H2 In = %s slpm, T In = %s C, Diameter = %s in]'%(hydrogenslpm,TempIn,diameter),fontsize = 16)
#        ax.set_zlim3d(0.00, 100.0)
#        ax.w_zaxis.set_major_locator(LinearLocator(10))
#        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.00f'))
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        fig.savefig(figureDirectory + 'COToMethane_H2_%s_slpm.png'%hydrogenslpm)
#        if showplots == 'yes':
#            show()
#        
#    def plot_data(self,figureDirectory = '/home/jordan/code/python/Figures',showplots = 'yes', temperature = 'yes',concentration = 'yes',conversion = 'yes'):
#        """
#            Plot Temperature, Concentration, and Carbon Conversion profiles depending on the 
#            optional arguments.  figureDirectory optional argument is a filepath to put generated figures.  
#        """
#        solids = Phase.phaseDictionary['Solids']
#        gas = Phase.phaseDictionary['Gas']
##        
#        figure(000 + self.simulationCount,figsize = (20,10))
#        suptitle('WallFlux Profile_RunID=_%s_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)),fontsize = 16) 
#        z_long = np.linspace(0,self.tube.length,len(self.wallflux))
#        plot(z_long,self.wallflux)
#        #plot(self.z,self.wallflux)
#        xlabel('Length [m]')
#        ylabel('Wall Flux [kW/m^2]')
#        savefig(figureDirectory + ('RunID_%s_WallFlux_Profile_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)))) 
#        
#        figure(050 + self.simulationCount,figsize = (20,10))
#        suptitle('Wall Temp Profile_RunID=_%s_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)),fontsize = 16) 
#        z_long = np.linspace(0,self.tube.length,len(self.wallTemp))
#        plot(z_long,self.wallTemp)
#        #plot(self.z,self.wallflux)
#        xlabel('Length [m]')
#        ylabel('Wall Temperature [C]')
#        savefig(figureDirectory + ('RunID_%s_WallFlux_Profile_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)))) 
#        
#        if temperature == 'yes':
#            #labels = []
#            colors = ['red','blue']
##            for phase in Phase.phaseList:
##                labels.append(phase.name)
#            #Temperature plot
#            figure(100 + self.simulationCount,figsize = (20,10))
#            i = 0
#            for phase in Phase.phaseList:
#                plot(self.z,self.phaseTemps[:,phase.temperatureIndex], label = phase.name, color = colors[i])
#                i += 1
#            
#            suptitle('Temperature Profile_RunID=_%s_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)),fontsize = 16) 
##            suptitle('TWall=%sK_TGasIn=%sC_P=%spsig_CH4In=%s_O2In=%s_H2In=%s.jpeg'
##                    %(round(self.tube.temperature),round(self.y0[-1]),round(self.gas.pressure*14.7/101327 - 14.7),
##                      round(self.n_0[self.species_indices['CH4']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['O2']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['H2']]*1000*60*24,3)),fontsize = 16) 
#            xlabel('Distance [m]',fontsize = 16)
#            ylabel('Temperature [K]',fontsize = 16)
#            legend(loc = 'upper right')
#            savefig(figureDirectory + ('RunID_%s_Temperature_Profile_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)))) 
##            savefig(figureDirectory +('TemperatureProfile_TWall=%sK_TGasIn=%sC_P=%spsig_CH4In=%s_O2In=%s_H2In=%s.jpeg'
##                    %(round(self.tube.temperature),round(self.y0[-1]),round(self.gas.pressure*14.7/101327 - 14.7),
##                      round(self.n_0[self.species_indices['CH4']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['O2']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['H2']]*1000*60*24,3)))) 
#        if concentration == 'yes':
#            figure(200 + self.simulationCount,figsize = (20,10))
#            NUM_COLORS = 15
#            cm = get_cmap('Set1')
#            colorIndex = 0
#            for species, compositions in self.mol_frac.items():
#                if compositions[-1] > 0.001 or compositions[0] > 0.01: 
#                    if species != 'N2':
#                        color = cm(1.*colorIndex/NUM_COLORS)  # color will now be an RGBA tuple
#                        plot(self.z,self.mol_fracNormalized[species],color = color, label = species,linewidth = 2.0)
#                        colorIndex += 1
#            suptitle('Concentration Profile_RunID=_%s_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)),fontsize = 16) 
#            xlabel('Distance [m]',fontsize = 16)
#            ylabel('Mole Fraction',fontsize = 16)
#            legend(loc='upper right')
#            savefig(figureDirectory + ('RunID_%s_Concentration_Profile_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7))))
#            #plot(self.z,self.mol_fracNormalized['H2O'], label = 'H2O',linewidth=2.0,color = 'magenta')
#            #plot(self.z,self.mol_frac['NO'], label = 'NO',linewidth=2.0,color = 'orange')
#            #plot(self.z,self.mol_frac['NO2'], label = 'NO2',linewidth=2.0,color = 'yellow')
#            #plot(self.z,self.mol_frac['N2O'], label = 'N2O',linewidth=2.0,color = 'magenta')
#            #plot(self.z,self.mol_frac['N2'], label = 'N2',linewidth=2.0,color = 'indigo')
#            
##            plot(self.z,self.mol_fracNormalized['H2'], label = 'H2',linewidth = 2.0,color = 'cyan')
##            plot(self.z,self.mol_fracNormalized['O2'], label = 'O2',linewidth = 2.0,color = 'red')
##            plot(self.z,self.mol_fracNormalized['CO'], label = 'CO',linewidth=2.0,color = 'blue')
##            plot(self.z,self.mol_fracNormalized['CO2'], label = 'CO2',linewidth = 2.0,color = 'purple')
##            plot(self.z,self.mol_fracNormalized['CH4'], label = 'CH4',linewidth = 2.0,color = 'green')
##            plot(self.z,self.mol_fracNormalized['C(gr)'],label = 'C(gr)',linewidth=2.0,color='orange')
##            plot(self.z,self.mol_fracNormalized['C2H2'], label = 'C2H2', linewidth=2.0,color='black')
##            plot(self.z,self.mol_fracNormalized['C2H4'], label = 'C2H4', linewidth = 2.0, color = 'brown')
#            
##            suptitle('TWall=%sK_TGasIn=%sC_P=%spsig_CH4In=%s_O2In=%s_H2In=%s.jpeg'
##                    %(round(self.tube.temperature),round(self.y0[-1]),round(self.gas.pressure*14.7/101327 - 14.7),
##                      round(self.n_0[self.species_indices['CH4']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['O2']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['H2']]*1000*60*24,3)),fontsize = 16) 
#             
##            savefig(figureDirectory +('ConcentrationProfile_TWall=%sK_TGasIn=%sC_P=%spsig_CH4In=%s_O2In=%s_H2In=%s.jpeg'
##                    %(round(self.tube.temperature),round(self.y0[-1]),round(self.gas.pressure*14.7/101327 - 14.7),
##                      round(self.n_0[self.species_indices['CH4']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['O2']]*1000*60*24,3),
##                      round(self.n_0[self.species_indices['H2']]*1000*60*24,3))))
#        if conversion == 'yes':
#            #Conversion plot
#            figure(300,figsize = (20,10))
#            plot(self.z,self.X)
#            suptitle('Conversion Profile_RunID=%s_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)),fontsize = 16) 
#            xlabel('Distance [m]',fontsize = 16)
#            ylabel('Conversion (%)',fontsize = 16)
#            savefig(figureDirectory + ('RunID_%s_Conversion_Profile_mBio=%slbhr_mSteam=%slbhr_mN2=%slbhr_TWall=%sK_TSteam=%sC_P=%spsig.jpeg'
#                    %(self.simulationCount,round(sum(solids.m_0.values())),round(gas.m_0['H2O']),round(gas.m_0['N2']),round(self.tube.temperature),round(gas.TempIn),round(gas.pressure*14.7/101327 - 14.7)))) 
#        if showplots != 'yes':
#            clf()
#            clf()
#            clf()
#
#        if showplots == 'yes':
#            show()
#            close()
#            close()
#            close()

class LSODAGas(LSODASimulator):
    """
        Class Instance which subclasses LSODASimulator.  This class instantiates a Gas Phase instance
        and runs the Base Class setInletConditions function and solveLSODA function.  It also includes
        its own energy balance function for a gas phase.
    """
    def __init__(self,tube,**kwargs):
        """
            kwargs: gas_filename,energyBalance
        """
        LSODASimulator.__init__(self,tube = tube)
        self.gas_filename = kwargs.get('gas_filename')
        self.energyBalance = kwargs.get('energyBalance')
        self.gas = Phase(name = 'Gas')
        
    def setInletConditions(self, **kwargs):
        #Create a cantera phase object here to do the mixing temperature calculation
        if self.simulationCount ==0:
            self.gas.cantera_phase = importPhase(self.gas_filename)
            self.gas.fileName = self.gas_filename
        Tmix = newton(self.calcMixingEnergy,x0 = np.mean(kwargs.get('T_gas').values()),args=(self.gas.cantera_phase,kwargs.get('m_gas'),kwargs.get('T_gas')))
        phaseDict = {self.gas:{'m_0':kwargs.get('m_gas'),'TempIn':Tmix,'pressure':kwargs.get('pressure'),'phaseMechanism':self.gas_filename}}
        LSODASimulator.setInletConditions(self,phaseDict)
    
    def massAndEnergyBalance(self,y,z,*args):
        sys.stdout.write("Simulation progress: %d%%   \r" % ((z/self.tube.length)*100) )
        A = self.tube.XSA                    #Tube Cross Sectional Area
        D = self.tube.diameter               #Tube Diameter
        T_wall = self.tube.temperature       #Tube Wall Temperature
        n_gas = y[0:-1]                      #Gas Molar Rates
        T_g = y[-1]                          #Gas Temperature
        pressure = self.gas.pressure         #Gas Pressure
        self.gas.cantera_phase.set(T=T_g, P=pressure, X=n_gas)
        self.gas.total_molar_flow_out = n_gas                           #Array of molar rate for each species
        H1 =self.gas.cantera_phase.enthalpies_RT()*R*T_g                #Gas Enthalpy
        Cp1 = self.gas.cantera_phase.cp_R()*R                           #Gas Heat Capacity
        r1 = self.gas.cantera_phase.netProductionRates()                #Gas Net Production Rates
        U = self.tube.U                                                 #Tube Heat Transfer Coefficient
        
        sourceTerm = D*U*(T_wall - T_g)                                 #Convection Heat source Term for energy balance
        if self.energyBalance == 'ConstantWallTemp':
            dT_dz = (D*U*(T_wall - T_g) - sum(H1*r1*A))/ sum(n_gas*Cp1)             
            dn_dz = A*r1
            dy = np.append(dn_dz, dT_dz)
        elif self.energyBalance == 'Adiabatic':
            dT_dz = -sum(H1*r1*A)/ sum(n_gas*Cp1) 
            dn_dz = A*r1
            dy = np.append(dn_dz, dT_dz) 
        elif self.energyBalance == 'Isothermal':
            dT_dz = 0
            dn_dz = A*r1
            dy = np.append(dn_dz, dT_dz)
        #Get wall flux and wall temp data
        self.length.append(z)
        self.wallTemp.append(T_wall-273.15)             #[C]
        self.wallflux.append(sourceTerm*(D/4.0)/1000.0) #[kW/m^2]
        return dy
    
class LSODASolidGas(LSODASimulator):
    """
        Class Instance which subclasses LSODASimulator for simulating solid-gas kinetics.  The subclass instantiates
        a solid and gas phase and uses its own energy balance for solid gas reactions and heat transfer  
    """
    def __init__(self,tube,**kwargs):
        """
            kwargs: gas_filename,solids_filename,energyBalance
        """
        LSODASimulator.__init__(self,tube = tube)
        self.solids = Phase(name = 'Solids')
        self.gas = Phase(name = 'Gas')
        
        self.solidSpecies = kwargs.get('SolidSpecies',['ASH','CELL','HCE','LIGH','LIGO','LIGC'])   #Selected Through GUI
        
        self.phi = 0.5                                  #Changed at Coupler level
        self.particleDiameter = 100E-6        #[meters] #Changed at Coupler level
        self.gasThermalConductivity  = 0.050  #[W/m-K]  #Changed at Coupler level
        self.solidDensity = 1000.0            #[kg/m^3] #Changed at Coupler level
        
        self.gas_filename = kwargs.get('gas_filename')
        self.solids_filename = kwargs.get('solids_filename')
        self.energyBalance = kwargs.get('energyBalance',"MultiPhase")
        
    def setInletConditions(self, **kwargs):
        """
            kwargs: m_solids,T_solids,pressure,m_gas
        """
        if self.simulationCount ==0:
            self.gas.cantera_phase = importPhase(self.gas_filename)
            self.gas.fileName = self.gas_filename
        Tmix = newton(self.calcMixingEnergy,x0 = np.mean(kwargs.get('T_gas').values()),args=(self.gas.cantera_phase,kwargs.get('m_gas'),kwargs.get('T_gas')))
        phaseDict = {self.solids:{'m_0':kwargs.get('m_solids'),'TempIn':kwargs.get('T_solids'),'pressure':kwargs.get('pressure'),
                                  'phaseMechanism':self.solids_filename},
                     self.gas:{'m_0':kwargs.get('m_gas'),'TempIn':Tmix,'pressure':kwargs.get('pressure'),'phaseMechanism':self.gas_filename}}
        LSODASimulator.setInletConditions(self,phaseDict)
        
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------     
        #PICK OUT SOLID AND GAS INDICES SPECIFIC TO ACTUAL GIVEN SPECIES.  (CELL,HCE,ect....)
        #solidSpeciesIndices AND gasSpeciesIndices ARE DICTIONARIES WITH THE SPECIES NAME MAPPED TO THE 
        #SPECIES INDEX.  UNLIKE THE SECTION ABOVE, THESE DICTIONARIES ONLY CONTAIN EITHER SOLID OR GAS 
        self.solidSpeciesIndices = {}
        self.gasSpeciesIndices = {}
        for species in self.solidSpecies:    #solidSpecies is a list of solids given by the user
            if species in self.solids.species_names:
                self.solidSpeciesIndices[species] = self.solids.species_names.index(species)
        for species in self.gas.species_indices.keys():
            if species not in self.solidSpeciesIndices.keys():
                self.gasSpeciesIndices[species] = self.gas.species_names.index(species)
    

    def solveLSODA(self,qInterpolator,zleft,zright):
        LSODASimulator.solveLSODA(self,qInterpolator,zleft,zright)
        #RETURN THE SOLID TEMPERATURE ARRAY FOR THE MONTE CARLO SIMULATION
        return self.phaseTemps[:,self.solids.temperatureIndex]    #return solid temperture for monte carlo calcultions
             
    def massAndEnergyBalance(self,y,z,qInterpolator,z0,zf):
        sys.stdout.write("Simulation progress: %d%%   \r" % ((z/self.tube.length)*100))
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #CUT THE INPUT ARRAY INTO RELEVANT VALUES  (SAME ROUTINE AS IN solveLSODA)
        n_solids = np.zeros(len(self.solids.species_indices.values()))
        n_gas = np.zeros(len(self.gas.species_indices.values()))      #array of zeros equal to the length of species array
        #GET APPROPRIATE SOLID AND GAS SPECIES FROM INPUT ARRAY
        n_solids[np.ix_(self.solidSpeciesIndices.values())] = y[np.ix_(self.solidSpeciesIndices.values())]
        n_gas[np.ix_(self.gasSpeciesIndices.values())] = y[np.ix_(self.gasSpeciesIndices.values())]
        n = y[0:-self.tempIndex]
        phaseTemps = y[-self.tempIndex:]
        T_s = phaseTemps[self.solids.temperatureIndex]                #Solids Temperature
        T_g = phaseTemps[self.gas.temperatureIndex]                   #Gas Temperature 

        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #SETUP NECESSARY DATA FOR THE ENERGY BALANCE TO SOLVE
        #Need solids and gas fractions to set correct pressures.  Overall pressure should stay constant.
        solids_fraction = sum(n_solids)/sum(n)
        gas_fraction = sum(n_gas)/sum(n)
        A = self.tube.XSA                    #Tube Cross Sectional Area
        D = self.tube.diameter               #Tube Diameter
        sigma = 5.67E-8                      #SI units for Stefan-Boltzmann constant
        epsilon = self.tube.emissivity       #Tube Emissivity
        
        pressure = self.gas.pressure             #System Pressure 
        m_solid = self.solids.total_flow           #[kg/s]
        #m_gas = self.gas.total_gas_flow                #[kg/s]
        rho_solid = self.solidDensity              #[kg/m^3]
        particleDiameter = self.particleDiameter   #[meters]
        k_g = self.gasThermalConductivity          #[W/m-K]Thermal Conductivity of gas
        phi = self.phi                             #Effectiveness Factor
        T_wall = self.tube.temperature             #Tube Wall Temperature
        heatflux = self.tube.heatflux/((D/4.0))
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #SET SOLID AND GAS OBJECT TO THE SYSTEM CONDITIONS IN ORDER TO GRAB REACTION AND THERMO PROPERTIES
        self.solids.cantera_phase.set(T=T_s, P=pressure*solids_fraction, X=n_solids)
        self.gas.cantera_phase.set(T=T_g, P=pressure*gas_fraction, X=n_gas)
        #The following is used for calculating heat duty
        self.solids.total_molar_flow_out = n_solids
        self.gas.total_molar_flow_out = n_gas
        
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #CALCULATE REACTION RATES, ENTHALPIES, AND HEAT CAPACITIES OF SOLID AND GAS PHASES
        r1 = self.solids.cantera_phase.netProductionRates()    #[kmol/m^3-s] reaction rate for each species in the solids
        r2 = self.gas.cantera_phase.netProductionRates()       #[kmol/m^3-s] reaction rate for each species in the gas
        r1_gas_only = np.zeros(len(self.gas.species_indices.values()))
        r1_gas_only[np.ix_(self.gasSpeciesIndices.values())] = r1[np.ix_(self.gasSpeciesIndices.values())]
        H1 = self.solids.cantera_phase.enthalpies_RT()*R*T_s     #[J/kmol] enthalpy of each gas component
        H2 = self.gas.cantera_phase.enthalpies_RT()*R*T_g        #[J/kmol] enthalpy of each solid component
        Cp1 = self.solids.cantera_phase.cp_R()*R                 #[J/kmol-K] molar heat capacities of each species
        Cp2 = self.gas.cantera_phase.cp_R()*R                    #[J/kmol-K] molar heat capacities of each species
        #-----------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------
        #CALCULATE ENERGY BALANCE SOURCE AND COUPLING TERMS
        #SOLVE THE ENERGY BALANCE CORRESPONDING TO THE USER SET PHASE
        sourceTerm = D/A*np.pi*sigma*epsilon*phi*(T_wall**4 - T_s**4)  #Radiation Source Term
        couplingTerm = (T_s-T_g)*12.0*k_g*pressure/R/T_g*m_solid/rho_solid/sum(n_gas)/(particleDiameter**2) #Energy Transfer from solids to gas Term
        
        if self.energyBalance == 'ConstantWallTemp':
            if self.MCCount == 0:
                dT1_dz = (sourceTerm - couplingTerm - sum(H1*r1))/sum(n_solids*Cp1/A)
                dT2_dz = (couplingTerm + sum(H1*r1_gas_only) - sum(H2*(r2+r1_gas_only)))/sum(n_gas*Cp2/A)  
            else:
                dT1_dz = (qInterpolator(z)[0] - couplingTerm - sum(H1*r1))/sum(n_solids*Cp1/A)   #qInterpolator(z)[0] because any value returned from a spline is a list
                dT2_dz = (couplingTerm + sum(H1*r1_gas_only) - sum(H2*(r2+r1_gas_only)))/sum(n_gas*Cp2/A)   
            dn_dz = A*(r1+r2)

        elif self.energyBalance =='Isothermal':
            dT1_dz = 0
            dT2_dz = 0
            dn_dz = A*(r1+r2)

        elif self.energyBalance == 'ConstantHeatFlux':
            sourceTerm = heatflux
            dT1_dz = (sourceTerm - couplingTerm - sum(H1*r1))/sum(n_solids*Cp1/A)
            dT2_dz = (couplingTerm + sum(H1*r1_gas_only) - sum(H2*(r2+r1_gas_only)))/sum(n_gas*Cp2/A)
            dn_dz = A*(r1+r2)
            T_wall = (sourceTerm/(D/A*np.pi*sigma*epsilon*phi) + T_s**4)**(1.0/4.0)
            
        if self.solids.temperatureIndex < self.gas.temperatureIndex:
            dy = np.append(dn_dz, [dT1_dz, dT2_dz])
        else:
            dy = np.append(dn_dz, [dT2_dz, dT1_dz])
        self.length.append(z)   
        self.wallTemp.append(T_wall-273.15)                #[C]
        self.wallflux.append(sourceTerm*(D/4.0)/1000.0)    #[kW/m^2]
        return dy
    
            
if __name__ == '__main__':
    tube1 = Tube(diameter = 1.5/39.37, length = 24.0/39.37, emissivity = 0.9, temperature = 1723)
    sim1 = LSODASolidGas(tube = tube1, gas_filename = "completeMechanism_gas.xml", solids_filename = "completeMechanism_solid.xml")
    sim1.setInletConditions(m_solids =  {'ASH':20.0,'CELL':44.0,'HCE':16.0,'LIGH':20.0/3.0,'LIGO':20.0/3.0,'LIGC':20.0/3.0}, m_gas = {'H2O':50.0,'CO2':16.53}, T_solids = 25, T_gas = 300, pressure = 50)
    #sim1.setTolerance(rtol = 1E-12, atol = 1E-12)
    sim1.solve()
    sim1.plot_data()
