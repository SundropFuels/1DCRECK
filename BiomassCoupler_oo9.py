'''
Created on Apr 28, 2011

@author: root
'''
#!/usr/bin/python2.7

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import csv
import BiomassSolver_oo9 as Solver
import time
import os
import sys
import subprocess


class simulationCoupler:
    """
        A class for coupling the monte carlo solver and LSODA mass and energy balance solver.  The initializer takes arguments for the two simulations.
        This class also contains a wrapper method for executing the monte carlo c++ code in the terminal.  The initializer defines initial tolerances for the 
        volumetric heat flux and temperature profiles between the two simulations.
    """

    def __init__(self,**kwargs):
        
        #self.energyBalance = energyBalance   #This doesn't actually go anywhere.  I have to run set simulator to update the simulator
        self.tube = Solver.Tube()   #Parameters set in the Runner
        simulatorClass = getattr(Solver,"LSODA%s"%kwargs.get('phaseType','Gas'))
        
        self.simulator = simulatorClass(self.tube,**kwargs)  #kwargs are arguments that can be passed to the subclasse's __init__ function
        self.MCSimulation = kwargs.get('MCSimulation',False)
        
        
        self.qtol = kwargs.get('qtol',0.1)
        self.Ttol = kwargs.get('Ttol',0.1)
        #print 'qtol is %s'%self.qtol
        #print 'Ttol is %s'%self.Ttol
        self.MCFilePath = kwargs.get('MCFilePath',"/home/jordan/code/c++/radiation/MCTest")
        self.MCInputFile = self.MCFilePath.rsplit('/',1)[0] + '/'+ kwargs.get('MCInputFile',"MCInput.csv")
        self.MCOutputFile = self.MCFilePath.rsplit('/',1)[0] + '/' + kwargs.get('MCOutputFile',"MCOutput.csv")

        self.zSlices = kwargs.get('zSlices',15)
        print self.zSlices
        bufferFraction = kwargs.get('BufferFraction',0.25)
        self.MCSlices = int(np.ceil(self.zSlices*bufferFraction)*2)
        self.PhiSlices = kwargs.get('PhiSlices',1.0)
        self.RadialSlices = kwargs.get('RadialSlices',1.0)

        #qold, Told, q_data_MC, and T_data_MC are initially set to be different so the tolerance check will clear
        self.qold = np.ones(self.zSlices + self.MCSlices)*10                  
        self.Told = np.ones(self.zSlices + self.MCSlices)*10
        self.q_data_MC = np.ones(self.zSlices + self.MCSlices)           
        self.T_data_MC = np.ones(self.zSlices + self.MCSlices)
        
        #Ends of the tube
        self.z0 = 0
        self.zf = self.tube.diameter
        
        #Define zleft and zright for the qinterpolator.  These are passed to the simulator to handle interpolating outside the range of montecarlo data returned.
        self.zleft = 0
        self.zright = 0
        #qInterpolator will be reset once the data for it are generated 
        self.qInterpolator = None  

        
    def runSimulations(self):
        """
            Setup each simulation with the given parameters and run them to return temperature data for the energy balance, and volumetric heat flux
            data from monte carlo to be input back to the energy balance
        """
        if self.MCSimulation == False:
            self.T_data_new = self.simulator.solveLSODA(self.qInterpolator,self.zleft,self.zright)  #Run LSODA solver and get array of temperature data (Tsolid)
            self.simulator.simulationCount += 1
        elif self.MCSimulation == True:
            convergence = False
            while convergence == False:
                self.T_data_new = self.simulator.solveLSODA(self.qInterpolator,self.zleft,self.zright)
                self.simulator.simulationCount += 1
                self.Tinterpz()
                self.writeMCInputFile()
                self.runMonteCarlo() 
                self.simulator.MCCount += 1   
                self.readMCOutputFile()
                self.qinterpz()
                self.removeMCFiles()
                convergence = self.toleranceCheck()
                
    def toleranceCheck(self):
        """
            Read the old and new temperature and heat flux data and test if their profiles fall within the specified tolerance.
        """
        qnew = self.q_data_MC
        Tnew = self.T_data_MC
        qold = self.qold
        Told = self.Told
        print "T old is %s" % Told
        print "T new is %s" % Tnew
        print "q old is %s" % qold
        print "q new is %s" % qnew

        #Strip out the beginning and ending slices that were added for MC
        qnew = qnew[(self.MCSlices/2):(self.MCSlices/2*-1)]
        Tnew = Tnew[(self.MCSlices/2):(self.MCSlices/2*-1)]
        qold = qold[(self.MCSlices/2):(self.MCSlices/2*-1)]
        Told = Told[(self.MCSlices/2):(self.MCSlices/2*-1)]
        if np.alltrue(abs((qold - qnew)/qold) < self.qtol) and np.alltrue(abs((Told - Tnew)/Told) < self.Ttol):
            return True
        else:
        
            #print 'q relative difference is' + str(abs((qold - qnew)/qold))
            #print 'T relative difference is' + str(abs((Told - Tnew)/Told))
            self.qold = self.q_data_MC
            self.Told = self.T_data_MC
            print 'RAWR!'
            return False


    def removeMCFiles(self):
        """
            Remove the csv files created for the MC simulation
        """
        try:
            print self.MCInputFile
            print self.MCOutputFile
            retcode1 = subprocess.call("rm -f %s" %self.MCInputFile,shell = True)
            retcode2 = subprocess.call("rm -f %s" %self.MCOutputFile,shell = True)
            if retcode1 < 0 or retcode2 < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode1, -retcode2
            else:
                print sys.stderr, "Child returned, File Removed" , retcode1, retcode2
        except OSError, e:
            print >>sys.stderr, "Couldn't Delete MC Files!:", e
            
    def runMonteCarlo(self):
        """
            Input the Monte Carlo file path (MCTest.cpp) and the csv file containing the temperature data and 
            run the simulation and return the necessary results from the csv
        """
        
        print file(self.MCInputFile).read()
        try:
            retcode = subprocess.call("%s -i %s -o %s" % (self.MCFilePath,self.MCInputFile,self.MCOutputFile),shell=True)
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode
            else:
                print >>sys.stderr, "Child returned", retcode
        except OSError, e:
            print >>sys.stderr, "Couldn't Run Monte Carlo Simulation!:", e

        
    def writeMCInputFile(self):
        """
            Write Temperature results to a csv file for the monte carlo simulation to read
        """
        self.Cell = np.arange(1,self.zSlices+self.MCSlices + 1)
        self.abs_coeff = np.ones(self.zSlices + self.MCSlices)*8.0
        self.scatt_coeff = np.ones(self.zSlices + self.MCSlices)*0.2

        inputfile = open(self.MCInputFile , 'w')
        writer = csv.writer(inputfile,escapechar = '"',quoting=csv.QUOTE_NONE,delimiter = ',')
        writer.writerow(["Length","Radius","Emissivity","TubeTemperature"])
        writer.writerow(["%s,%s,%s,%s"% (self.MC_Length,self.tube.diameter/2.0,self.tube.emissivity,self.tube.temperature)])
        writer.writerow(["RadialSlices","Phi Slices","Z Slices"])
        writer.writerow(["%s,%s,%s"% (self.RadialSlices,self.PhiSlices,(self.zSlices+self.MCSlices))])
        writer.writerow(["Cell#","abs_coeff","scatt_coeff","T"])
        
        for slice in range(self.zSlices+self.MCSlices):
            writer.writerow(["%s,%s,%s,%s"% (self.Cell[slice],self.abs_coeff[slice],self.scatt_coeff[slice],self.T_data_MC[slice])])                  
        inputfile.close()
        


    def readMCOutputFile(self):
        """
            Read the Temperature Data from the Monte Carlo output file and return an array of q data
        """

        outputFile = open(self.MCOutputFile, 'r')
        filereader = csv.reader(outputFile)
        print file(self.MCOutputFile).read()
        filereader.next()
        z = []
        q = []
    
        for row in filereader:
            z.append(float(row[2]))
            q.append(float(row[3]))
    
        outputFile.close()
        self.z_MC = np.asarray(z)
        self.q_data_MC = np.asarray(q)


    def qinterpz(self):
        """
            This method interpolates the q values returned from monte carlo into a z-dependent function for input in the energy balance
        """
        inside_z = self.z_MC[len(self.leadingslices):-len(self.trailingslices)]
        translated_z = inside_z - self.z0
        self.zleft = translated_z[0]
        self.zright = translated_z[-1]
        inside_q = self.q_data_MC[len(self.leadingslices):-len(self.trailingslices)]
        self.qInterpolator = InterpolatedUnivariateSpline(translated_z,inside_q,k=3)
        

#        figure(1)
#        scatter(translated_z,inside_q,color = 'green')
#        title('Monte Carlo Data Points and Spline')
#        xlabel('Length')
#        ylabel('Volumetric Heat Flux [W/m^3]')
#
#        znew = self.simulator.z
#        plot(znew,self.qInterpolator(znew))
#        show()

    def Tinterpz(self):
        """
            This method interpolates the Temperature values returned from the energy balance to the slices the MC program takes in
        """
        T_interp = InterpolatedUnivariateSpline(self.simulator.z,self.T_data_new,k=3)
        sliceLength = self.tube.length/self.zSlices
        #print sliceLength
        inside_z = np.linspace(0 + sliceLength/2,self.tube.length-sliceLength/2, self.zSlices)
        inside_T = T_interp(inside_z)
        #print np.diff(inside_z)[0]
        self.MC_Length = self.tube.length + self.MCSlices*sliceLength

        z0prime = 0
        self.z0 = 0 + self.MCSlices/2*sliceLength
        zfprime = self.MC_Length
        zf = self.tube.length + self.MCSlices/2*sliceLength

        #Translate inner slices to insert leading slices at zero
        inside_z = inside_z + self.z0
        self.leadingslices = np.arange(z0prime,self.z0,sliceLength)
        self.trailingslices = np.arange(zf + sliceLength,zfprime+sliceLength,sliceLength)
        leading_T = np.ones(len(self.leadingslices))*inside_T[0]
        trailing_T = np.ones(len(self.trailingslices))*inside_T[-1]

        
        self.z_MC = np.hstack((self.leadingslices,inside_z,self.trailingslices))
        self.T_data_MC = np.hstack((leading_T,inside_T,trailing_T))
        self.T_data_MC = np.floor(self.T_data_MC)
 
        #Test the temperature spline fit
        #figure(1)
        #plot(self.simulator.z,self.T_data_new,color = 'green')
        #title('LSODA Data Plot and Interpolated MC Points')
        #xlabel('length')
        #ylabel('Temperature')

        #z_test = self.z_MC
        #scatter(self.z_MC,T_interp(self.z_MC))
        #show()
        

    def setqTtolerance(self, qtol = 0.01, Ttol = 0.01):
        """
            Set the Temperature and Heat Flux Tolerance for the coupled simulations
        """
        self.qtol = qtol
        self.Ttol = Ttol
        
        
    def setIntegrationTolerances(self,rtol = 1E-12,atol = 1E-12):
#        toleranceLength = len(self.simulator.y0)
#        rtol = np.ones(toleranceLength)*10E-12
#        atol = np.ones(toleranceLength)*10E-12
#        
#        selectTolerance = 1E-11
# 
#        rtolSpecies = {}
#        atolSpecies = {}
#        tolSpecies = {}
#        for species, index in self.simulator.species_indices.items():
#            rtolSpecies[species] = massRelative
#            atolSpecies[species] = massAbsolute    
#            tolSpecies[species] = massAbsolute
#
#        #tol['C8H10O3'] = selectTolerance
#        #tol['C3H2'] = selectTolerance
#        tolSpecies['CO'] = selectTolerance
#        tolSpecies['CO2'] = selectTolerance
#        tolSpecies['CH4'] = selectTolerance
#        tolSpecies['H2'] = selectTolerance
#        tolSpecies['H2O'] = selectTolerance
#
#        rtol[-2:] = energyRelative
#        atol[-2:] = energyAbsolute
#
#        rtol[0:-2] = tolSpecies.values()
#        atol[0:-2] = tolSpecies.values()


        self.simulator.setTolerance(rtol = rtol, atol = atol)

    def setIntegrationSteps(self,steps = 10000):
        self.simulator.steps = steps

    def setTube(self,diameter = 3.5, length = 7, emissivity = 0.9, temperature = 1400):
        """
            Set tube parameters.  (diameter, length, emissivity, temperature)
        """
        self.simulator.tube.setTube(diameter = diameter, length = length, emissivity = emissivity, temperature = temperature)
    
    def setDiameter(self,diameter = 3.5):
        self.simulator.tube.setDiameter(diameter = diameter)
    
    def setLength(self,length = 7.0):
        self.simulator.tube.setLength(length = length)
    
    def setEmissivity(self,emissivity = 1.0):
        self.simulator.tube.setEmissivity(emissivity = emissivity)
    
    def setWallTemperature(self,temperature):
        self.simulator.tube.setWallTemperature(temperature = temperature)
    def setWallFlux(self,heatflux):
        self.simulator.tube.setWallFlux(heatflux = heatflux)
    
    def setHeatTransferCoefficient(self,U = 100):
        self.simulator.tube.setHeatTransferCoefficient(U = U) 
        
    def setParticleDiameter(self,diameter = 200E-6):
        self.simulator.particleDiameter = diameter
         
    def setPhi(self,phi):
        self.simulator.phi = phi
    def setGasThermalConductivity(self,k_g):
        self.simulator.gasThermalConductivity = k_g
    def setSolidDensity(self,rho_solid):
        self.simulator.solidDensity = rho_solid
    
        
#    def setSimulator(self,gas_filename = "completeMechanism_gas.cti", solids_filename = "completeMechanism_solid.cti",energyBalance = 'MultiPhase',**kwargs):
#        """
#           
#        """
#        self.simulator.gas_filename = kwargs.get("gas_filename",None)
#        self.simulator.solids_filename = kwargs.get("solids_filename",None)
#        self.simulator.energyBalance = kwargs.get("energyBalance",None)

        
    def setInletConditions(self,**kwargs):
        """
            Set system inlet conditions
            kwargs include: m_gas, m_solids,T_gas,T_solid, and pressure
        """
        self.simulator.setInletConditions(**kwargs)
        
    def setNormalizedSpecies(self,normalizedSpecies):
        self.simulator.normalizedSpecies = normalizedSpecies
    
    def plot_data(self,figureDirectory = '/home/jordan/code/python/Figures',showplots = 'yes', temperature = 'yes',concentration = 'yes',conversion = 'yes'):
        self.simulator.plot_data(figureDirectory = figureDirectory,showplots = showplots,temperature = temperature,concentration = concentration,conversion = conversion)

    def write_data(self,filename):
        self.simulator.write_data(filename = filename)
        
    def writeProfileData(self,filename):
        self.simulator.writeProfileData(filename = filename)

#    def store_data(self):
#        self.simulator.store_data()
#        
#    def clear_data(self):
#        self.simulator.clear_data()

#    def readExperimentalData(self,inputfile):
#        inletData = self.simulator.readExperimentalData(inputfile)
#        return inletData

#    def plotExperimentalComparisons(self):
#        pass
#
#    def plot_gasModel(self,showplots = 'yes',selectSpecies = 'CH4'):
#        self.simulator.plot_gasModel(showplots = showplots, selectSpecies = selectSpecies)
    
    def plot3D(self,X,Y,**kwargs):
        self.simulator.plot3D(X,Y,**kwargs)














    
