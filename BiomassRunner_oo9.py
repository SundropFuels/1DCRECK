

"""
    Biomass Running Module.  This Module's purpose is to grab data from the Biomass GUI Callback Function and 
    set the data to be run.  It does this by creating a nested loop of all possible configurations from the input.
"""

#!/usr/bin/python2.7
import BiomassCoupler_oo9 as Coupler
import BiomassPlottingToolbox as Plotter
from itertools import product
import csv
import numpy as np

class SimulationRunner():
    """
        Class for Setting and Running the simulation programs
    """
    #Pull in options from the Biomass GUI
    def __init__(self,runDictionary):
        """
            Get all of the parameters from the GUI and right a csv file with a row corresponding to each possible
            run condition
        """
        self.generalParameters = runDictionary['GeneralParameters']
        self.MCParameters = runDictionary['MCParameters']
        self.runParameters = runDictionary['RunParameters']
        self.solidFlowRates = runDictionary['SolidFlowRates']
        self.gasFlowRates = runDictionary['GasFlowRates']
        self.gasTemperatures = runDictionary['GasTemperatures']
        self.tubeParameters = runDictionary['TubeParameters']
        self.plottingProfiles = runDictionary['PlottingProfiles']
        self.comparisonProfiles = runDictionary['ComparisonProfiles']
        #WRITE RUN PARAMETERS TO CSV IF NOT USING A PRECREATED SPREADSHEET
        if self.generalParameters['UseParametersFromSpreadsheet'] == False:
        #BUILD A NESTED TABLE FOR EACH RUN THROUGH THE SIMULATION IF PARAMETERS WERE ENTERED USING THE GUI
            parameterFile = open(self.generalParameters['ParameterSpreadsheetName'],'wb')
            parameterwriter = csv.writer(parameterFile,delimiter = ',',quoting = csv.QUOTE_NONNUMERIC)
            headerwriter = csv.writer(parameterFile,delimiter = ',', quoting = csv.QUOTE_ALL)
            headerwriter.writerow(np.hstack(('RunID',self.runParameters.keys(),self.tubeParameters.keys(),
                                         self.solidFlowRates.keys(),self.gasFlowRates.keys(),self.gasTemperatures.keys())))

            masterlist = [self.runParameters[k] for k in self.runParameters.keys()]
            masterlist.extend([self.tubeParameters[k] for k in self.tubeParameters.keys()])
            masterlist.extend([self.solidFlowRates[k] for k in self.solidFlowRates.keys()])
            masterlist.extend([self.gasFlowRates[k] for k in self.gasFlowRates.keys()])
            masterlist.extend([self.gasTemperatures[k] for k in self.gasTemperatures.keys()]) 
#        for a,b,c,d,e,f,g,h in product(*np.concatenate(([self.runParameters[k] for k in self.runParameters.keys()],
#                                               [self.tubeParameters[k] for k in self.tubeParameters.keys()]))):
            RunID = 1      
            for values in product(*masterlist):     #product from itertools.  Runs a cartesian product and effectively hits all the required nested loops
                parameterwriter.writerow(np.hstack((RunID,values)))
                RunID += 1
            parameterFile.close()
        
#        methaneslpm = 10.0  #slpm
#        mMethane = methaneslpm/24.0*16.0*60.0/1000.0*2.204  #lb/hr
    def runSimulation(self):    
        """
            This method executes the simulation.  It grabs the parameters from the parameter spreadsheet 
            and runs through each set of run conditions, writing the results to .csv files
        """
        #READ THE PARAMETER FILE AND SET THE CONDITIONS
        parameterFile = open(self.generalParameters['ParameterSpreadsheetName'], 'rb')
        parameterReader = csv.reader(parameterFile,delimiter = ',',quoting = csv.QUOTE_NONNUMERIC)
        parameterHeader = parameterReader.next()
        parameterIndexDict = {}
        self.parameterDataDict = {}
        for i in range(len(parameterHeader)):
            parameterIndexDict[parameterHeader[i]] = i 
            self.parameterDataDict[parameterHeader[i]] = []
        for row in parameterReader:
            for key,value in parameterIndexDict.items():
                self.parameterDataDict[key].append(row[value])
        parameterFile.close()
        
        #Make a dictionary of solid species by selecting the solids pulled from the spreadsheet.
        #This dictionary could map solids to either mass flow rates or mass fractions
        solidsDict = dict((k.strip('solid '),v) for (k,v) in self.parameterDataDict.items() if 'solid ' in k)
        totalSolidsFlow = self.parameterDataDict['Total Solids Flow [lb/hr]']
        #Set the solid flowrates in the solidsDict based on a total solids flow if it isn't blank
        if '' not in self.parameterDataDict['Total Solids Flow [lb/hr]']:
            for key, value in solidsDict.items():
                solidsDict[key] = [value[i]*totalSolidsFlow[i] for i in range(len(totalSolidsFlow))]
        solidSpecies = solidsDict.keys()    #List of solid species the simulator uses to determine what is a solid
        
        #Make a dictionary mapping gas species to their respective inlet flow rates
        gasesDict = dict((k.strip('gas '),v) for (k,v) in self.parameterDataDict.items() if 'gas ' in k) 
        #Make a dictionary mapping gas species to their respective inlet temperatures
        gasTempDict = dict((k.strip(' GasTemperature'),v) for (k,v) in self.parameterDataDict.items() if ' GasTemperature' in k)          
        #Instantiate a a coupler and set up the simulation
        coupler = Coupler.simulationCoupler(MCSimulation = self.generalParameters['MCSimulation'],
                                            phaseType = self.generalParameters['PhaseType'],MCFilePath = self.MCParameters['MCFilePath'],
                                            MCInputFile = self.MCParameters['MCInputFileName'],MCOutputFile = self.MCParameters['MCOutputFileName'],
                                            qtol = self.MCParameters['MCqTolerance'],
                                            Ttol = self.MCParameters['MCTempTolerance'],energyBalance = self.generalParameters['EnergyBalance'],gas_filename = self.generalParameters['GasFileName'],
                                            solids_filename = self.generalParameters['SolidFileName'],SolidSpecies = solidSpecies)

        coupler.setIntegrationSteps(self.generalParameters['IntegrationSteps'])
        coupler.setIntegrationTolerances(rtol = self.generalParameters['RelativeTolerance'],atol = self.generalParameters['AbsoluteTolerance'])
        coupler.setNormalizedSpecies(normalizedSpecies = self.generalParameters['NormalizedSpecies'])
        #Run a simulation for each set of conditions
        for i in range(len(self.parameterDataDict['RunID'])):
            coupler.setInletConditions(m_gas = dict((k,v[i]) for (k,v) in gasesDict.items()), T_gas = dict((k,v[i]) for (k,v) in gasTempDict.items()),
                                       m_solids = dict((k,v[i]) for (k,v) in solidsDict.items()), T_solids = self.parameterDataDict['Solid Temperatures [C]'][i],
                                       pressure = self.parameterDataDict['Pressures [psig]'][i])
            coupler.setDiameter(diameter = self.parameterDataDict['Diameter [in]'][i])
            coupler.setLength(length = self.parameterDataDict['Length [ft]'][i])
            coupler.setEmissivity(emissivity = self.parameterDataDict['Emissivity'][i])
            coupler.setWallTemperature(temperature = self.parameterDataDict['Wall Temperatures [C]'][i])
            coupler.setWallFlux(heatflux = self.parameterDataDict['Wall Fluxes [W/m^2]'][i])
            coupler.setHeatTransferCoefficient(U = self.parameterDataDict['Overall Heat Transfer Coefficient [W/m^2 - K]'][i])
            
            coupler.setParticleDiameter(diameter = self.parameterDataDict['Particle Diameters [microns]'][i]*10**-6.0)
            coupler.setPhi(phi = self.parameterDataDict['Effectiveness Factor [phi]'][i])
            coupler.setGasThermalConductivity(k_g = self.parameterDataDict['Gas Thermal Conductivity [W/m - K]'][i])
            coupler.setSolidDensity(rho_solid = self.parameterDataDict['Solid Density [kg/m^3]'][i])

            coupler.runSimulations()  
            coupler.write_data(self.generalParameters['ExitDataName'])
            coupler.writeProfileData(self.generalParameters['ProfileDataFileName'])
        Coupler.Solver.clearPhase()   #Calls a Solver module function.  Allows multiple simulations in a GUI session.
            
    def plotData(self):
        """
            This method tells the plottingtoolbox was to plot.  It tells it where to get the profile data
            and exit data and exactly what to plot
        """
        plotter = Plotter.PlottingToolbox(exitFile = self.generalParameters['ExitDataName'],profileFile = self.generalParameters['ProfileDataFileName'],)
        for profile, flag in self.plottingProfiles.items():
            if flag == True:
                plotter.plotProfile(profile,figureDirectory = self.generalParameters['SaveFigureDirectory'])

        for datalist in self.comparisonProfiles.values():
            plotter.plot2D(datalist,figureDirectory = self.generalParameters['SaveFigureDirectory'])     
            
    



