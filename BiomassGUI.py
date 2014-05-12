'''
Created on Apr 28, 2011

@author: root
'''
#!/usr/bin/python2.7


import pygtk
pygtk.require('2.0')
import gtk
import BiomassRunner_oo9 as BiomassRunner


class WidgetConstructor:
    """
         WidgetConstructor acts as a base class, but really, all of its functions could be module level.  
         It provides a set of functions that can be called by the GUI program to set up widgets or grab data
    """
    def makeBoxGroup(self,hargs = (),vargs=()):
        """
            makeBoxGroup packs the given vargs into a verticle box, and then packs that with the hargs into 
            a horizontal box
        """
        VBox = gtk.VBox(homogeneous = False,spacing = 10)
        HBox = gtk.HBox(homogeneous = False,spacing = 10)
        for value in vargs:
            VBox.pack_start(value,expand = True,fill = True,padding = 0)
        for value in hargs:
            HBox.pack_start(value,expand = False,fill = True,padding = 0)
        HBox.pack_start(VBox,expand = True,fill = True,padding = 0)
        VBox.show()
        HBox.show()
        return HBox
    
    def makeLabel(self,text):
        label = gtk.Label(text)
        label.show()
        return label
    
    #file selection callback
    def file_ok_sel(self,ok_button,*args):
        """
            This function sets an entry field after clicking ok on the file selection dialog button
        """
        filename = self.fileSelector.get_filename()  #Send this to the appropriate place
        for entry in args:
            entry.set_text(filename)
        
    def openFileSelector(self,widget,title,*args):
        self.fileSelector = gtk.FileSelection(title)
        #fileSelector.connect("destroy",lambda wid: gtk.main_quit())
        self.fileSelector.cancel_button.connect("clicked",lambda wid:self.fileSelector.destroy())
        self.fileSelector.ok_button.connect("clicked", self.file_ok_sel,*args)
        self.fileSelector.ok_button.connect("clicked", lambda wid:self.fileSelector.destroy())
        self.fileSelector.show()
    
    #make button callback
    def makeButton(self,type,**kwargs):
        button = getattr(gtk,"%s"%type)(**kwargs)  
        button.show()
        return button
    #make text entry box callback
    def makeTextEntry(self,text,editable):
        entry = gtk.Entry(max = 0)
        entry.set_text(text)
        entry.set_editable(editable)
        entry.show()
        return entry

    #make combo box callback
    def makeComboBox(self,strings):
        model = gtk.ListStore(str)
        combobox = gtk.ComboBox(model = None)
        cell = gtk.CellRendererText()
        for string in strings:
            model.append([string])
        combobox.set_model(model)
        combobox.pack_start(cell,True)
        combobox.add_attribute(cell, 'text',0)
        return combobox
    
    def setComboBox(self,combobox,strings):
        model = gtk.ListStore(str)
        for string in strings:
            model.append([string])
        combobox.set_model(model)
    
    def getComboBoxText(self,combobox):
        model = combobox.get_model()
        active = combobox.get_active()
        return model[active][0] 
    
    def getModelText(self,model,columns,collapse):
        """
            This Method returns the values in a model as a 2D list.  Each List contains the values for 
            a row in the model.
        """
        stringFormatFunction = collapse and (lambda s: " ".join(s.split(','))) or (lambda s: s)
        values = []
        for i in range(len(model)):
            iter = model.get_iter(i)
            values.append(list(model.get(iter,*columns)))
        newvalues = [[stringFormatFunction(value) for value in valueset] for valueset in values] 
        return newvalues

class GasifierGUI(WidgetConstructor):  
    """
        The GasifierGUI Class contains all the methods for setting up the GUI pages, grabbing the data from each page,
        and passing the parameters along to the Runner to launch the simulation
    """
    def runSimulation(self,widget,data = None):
        """
            Run the simulation by Instantiating a Runner and giving it parameters
        """
        self.getRunData()
        runner = BiomassRunner.SimulationRunner(self.runParametersDict)
        if self.runParametersDict['GeneralParameters']['UsePreviousSimulation'] == False:
            runner.runSimulation()
        runner.plotData()
        
    def getRunData(self):
        """
            Collect all of the Run Data from the GUI.  Defaults are used if information isn't provided in certain fields.
            This method interacts with all of the data members created when setting up the GUI pages
        """
        #Dictionary Containing parameters to be passed to the Runner.  Starts with Defaults.
        self.runParametersDict = {'GeneralParameters':{'MCSimulation':False,'UsePreviousSimulation':False,'LoadExperimentalData':False,
                        'UseParametersFromSpreadsheet':False,'ParameterSpreadsheetName':'/home/jordan/code/python/parameters.csv',
                        'ProfileDataFileName':'/home/jordan/code/python/profileResults.csv',
                        'ExitDataName':'/home/jordan/code/python/exitResults.csv',
                        'ExperimentalFileName':None,'PhaseType':'Gas','EnergyBalance':'Isothermal',
                        'GasFileName':'/usr/local/cantera/data/gri30.cti',
                        'SolidFileName':'/usr/local/cantera/data/gri30.cti',
                        'SaveFigureDirectory':'/home/jordan/code/python/Figures','AbsoluteTolerance':1E-12,
                        'RelativeTolerance':1E-12,'IntegrationSteps':100,'NormalizedSpecies':[]},
                        'MCParameters':{'MCqTolerance':1.0,'MCTempTolerance':1.0,'MCZSlices':15,'MCPhiSlices':1,
                        'MCRadialSlices':1.0,'MCBufferFraction':0.2,'MCInputFileName':'MCInput.csv',
                        'MCOutputFileName':'MCOutput.csv','MCFilePath':'/home/jordan/code/c++/TestMC'},
                        'RunParameters':{'Total Solids Flow [lb/hr]':[None]}}
        
        #Dictionary that is populated with all of the data in the GUI
        self.GUIDict = {}
        GeneralParametersDict = {'MCSimulation':self.MCButton_On.get_active(),'UsePreviousSimulation':self.useStoredData_On.get_active(),
                        'LoadExperimentalData':self.useExperimentalData_On.get_active(),
                        'UseParametersFromSpreadsheet':self.useParameterSpreadsheet_On.get_active(),
                        'ParameterSpreadsheetName':self.parameterDataEntry.get_text(),
                        'ProfileDataFileName':self.storedDataEntry.get_text(),
                        'ExitDataName':self.dataFileEntry.get_text(),
                        'ExperimentalFileName':self.experimentalDataEntry.get_text(),'PhaseType':self.getComboBoxText(self.phaseCombo),
                        'EnergyBalance':self.getComboBoxText(self.energyBalanceCombo),
                        'GasFileName':self.gasFileEntry.get_text(),'SolidFileName':self.solidFileEntry.get_text(),
                        'SaveFigureDirectory':self.saveFigureEntry.get_text(),'AbsoluteTolerance':float(self.atolEntry.get_text()),
                        'RelativeTolerance':float(self.rtolEntry.get_text()),
                        'IntegrationSteps':float(self.stepsEntry.get_text())}
        MCDict = {'MCqTolerance':float(self.qToleranceEntry.get_text()),'MCTempTolerance':float(self.tempToleranceEntry.get_text()),
                        'MCZSlices':int(self.zSlicesEntry.get_text()),'MCPhiSlices':int(self.phiSlicesEntry.get_text()),
                        'MCRadialSlices':int(self.radialSlicesEntry.get_text()),'MCBufferFraction':float(self.bufferFractionEntry.get_text()),
                        'MCInputFileName':self.MCInputFileEntry.get_text(),'MCOutputFileName':self.MCOutputFileEntry.get_text(),
                        'MCFilePath':self.MCFilePathEntry.get_text()}
        
        self.GUIDict['GeneralParameters'] = GeneralParametersDict
        self.GUIDict['MCParameters'] = MCDict
        #Get Data from the Simulator Run Parameters Model
        parameters = self.getModelText(self.parameterModel,(0,1),collapse = True)
        parametersDict = {}
        for value in parameters:
            parametersDict[value[0]] = [float(s) for s in value[1].split()]
        self.GUIDict['RunParameters'] = parametersDict
        #Get Solid Flow Rate Data    
        solidFlowRates = self.getModelText(self.solidFlowRateModel,(0,1),collapse = True)
        solidFlowRatesDict = {}
        for value in solidFlowRates:
            if value[0] != 'Select a Solid':
                solidFlowRatesDict['solid '+value[0]] = [float(s) for s in value[1].split()]
        self.GUIDict['SolidFlowRates'] = solidFlowRatesDict
        #Get Gas Flow Rate Data and Species that are Normalized for plotting
        gasFlowRates = self.getModelText(self.gasFlowRateModel,(0,1,2),collapse = True)
        normalizedSpecies = self.getModelText(self.gasFlowRateModel,(0,3), collapse = False)
        gasFlowRatesDict = {}
        gasTemperaturesDict = {}

        for value in gasFlowRates:
            if value[0] != 'Select a Gas':
                gasFlowRatesDict['gas '+value[0]] = [float(s) for s in value[1].split()]
                gasTemperaturesDict[value[0] + ' GasTemperature'] = [float(s) for s in value[2].split()]
        self.GUIDict['GasFlowRates'] = gasFlowRatesDict
        self.GUIDict['GasTemperatures'] = gasTemperaturesDict
        self.GUIDict['GeneralParameters']['NormalizedSpecies'] = [value[0] for value in normalizedSpecies if value[1] == True]
        #Get Tube Parameter Data
        tubeParameters = self.getModelText(self.tubeParametersModel,(0,1),collapse = True)
        tubeParametersDict = {}
        for value in tubeParameters:
            tubeParametersDict[value[0]] = [float(s) for s in value[1].split()]
        self.GUIDict['TubeParameters'] = tubeParametersDict
        
        #Get Plotting Profile Options
        plottingProfiles = self.getModelText(self.profileModel,(0,1),collapse = False)
        plottingProfilesDict = {}
        
        for value in plottingProfiles:
            plottingProfilesDict[value[0]] = value[1]
            
        self.GUIDict['PlottingProfiles'] = plottingProfilesDict
        
        #Get Selected Plotting Comparisons
        plottingComparisons = self.getModelText(self.comparisonModel,(0,1,2,3),collapse = True)
        print plottingComparisons
        plottingComparisonsDict = {}
        count = 0
        for value in plottingComparisons:
            if value[0] != 'X Data' and value[1]!= 'Y Data':
                plottingComparisonsDict[count]=(value[0],value[1],value[2],value[3])
                count += 1
        self.GUIDict['ComparisonProfiles'] = plottingComparisonsDict

        
        #Read in data to the run dictionary.  If a value is blank, a default will be used instead.
        for key1, value1 in self.GUIDict.items():
            if key1 not in self.runParametersDict.keys():
                self.runParametersDict[key1] = {}
            for key2, value2 in value1.items():    
                if self.GUIDict[key1][key2] != '' and self.GUIDict[key1][key2] != []:
                    self.runParametersDict[key1][key2] = value2
        
    def __init__(self):
        """
            The initializer sets up the window and connect the 'delete_event' signal to it.
            I also sets up the exterior table the notebook is attached to.
        """
        #Create window and notebook
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", lambda wid,event: gtk.main_quit())
        #window.connect("destroy", lambda wid: gtk.main_quit())
        #window.connect("delete_event", lambda a1,a2:gtk.main_quit())
        self.window.set_border_width(10)
        self.exteriorTable = gtk.Table(rows = 3,columns = 6,homogeneous = False)
        self.window.add(self.exteriorTable)

    def createNotebook(self):
        # Create a new notebook, place the position of the tabs
        notebook = gtk.Notebook()
        notebook.set_tab_pos(gtk.POS_LEFT)
        self.exteriorTable.attach(notebook, left_attach = 0,right_attach = 6,top_attach = 0,bottom_attach = 1)
        notebook.show()
        self.show_tabs = True
        self.show_border = True
        #Create page in notebook for each section of input
        self.frameDict = {}
        for i in range(5):
            pages = ["General","MC","Simulator","Tube","Plotting"]
            frameNames = ["General Information", "MC Information", "Simulator Data","Tube Data","Plotting Data"]
            frame = gtk.Frame(frameNames[i])
            self.frameDict[pages[i]] = frame
            frame.set_border_width(10)
            frame.set_size_request(1300, 700)
            frame.show()
            label = gtk.Label(pages[i])
            notebook.append_page(frame, label)
        
        # Create buttons for cycling through pages or closing
        button = gtk.Button("close")
        button.connect("clicked", lambda w: gtk.main_quit())
        self.exteriorTable.attach(button, 0,1,1,2)
        button.show()

        button = gtk.Button("next page")
        button.connect("clicked", lambda w: notebook.next_page())
        self.exteriorTable.attach(button, 1,2,1,2)
        button.show()

        button = gtk.Button("prev page")
        button.connect("clicked", lambda w: notebook.prev_page())
        self.exteriorTable.attach(button, 2,3,1,2)
        button.show()
        
        button = gtk.Button("Run Simulation")
        button.connect("clicked", self.runSimulation)
        self.exteriorTable.attach(button,3,4,1,2)
        button.show()

        self.exteriorTable.show()
        self.window.show()
        
    #GENERAL MENU OPTIONS
    def setupGeneralPage(self):
        def changedComboBox(combobox):
            if self.getComboBoxText(combobox) == "SolidGas":
                self.setComboBox(self.energyBalanceCombo,["ConstantWallTemp","Isothermal","ConstantHeatFlux"])
            elif self.getComboBoxText(combobox) =="Gas":
                self.setComboBox(self.energyBalanceCombo,["ConstantWallTemp","Isothermal","Adiabatic"])
        table = gtk.Table(rows = 5, columns = 5,homogeneous = False)
        self.frameDict["General"].add(table)
        table.show()
        #MC OPTIONS
        title = self.makeLabel("Monte Carlo On/Off?")
        table.attach(title,0,1,0,1)
        MClabel = self.makeLabel("MC Simulation")
        self.MCButton_On = self.makeButton(group = None, type = "RadioButton",label = "On")
        MCButton_Off = self.makeButton(group = self.MCButton_On,type = "RadioButton",label = "Off")
        MCButton_Off.set_active(1)
        MCBoxGroup = self.makeBoxGroup(hargs = (MClabel,),vargs = (self.MCButton_On,MCButton_Off))
        table.attach(MCBoxGroup,0,1,1,2)
        #MECHANISM DATA
        title = self.makeLabel("Mechanism Data")
        table.attach(title,3,4,0,1)
        #PHASE TYPE
        phaseLabel = self.makeLabel("Phase Type     ")
        self.phaseCombo = self.makeComboBox(["Gas","SolidGas"])
        self.phaseCombo.set_active(0)
        phaseBoxGroup = self.makeBoxGroup(hargs = (phaseLabel,self.phaseCombo))
        self.phaseCombo.connect('changed',changedComboBox)
        table.attach(phaseBoxGroup,3,4,1,2)
        #ENERGY BALANCE
        energyBalanceLabel = self.makeLabel("Energy Balance")
        self.energyBalanceCombo = self.makeComboBox(["Isothermal"])
        energyBalanceBoxGroup = self.makeBoxGroup(hargs = (energyBalanceLabel,self.energyBalanceCombo))
        table.attach(energyBalanceBoxGroup,3,4,2,3)
        #FILE NAMES
        gasFileLabel = self.makeLabel("Gas File Name ")
        gasFileButton = self.makeButton(type = "Button", label = "...")
        self.gasFileEntry = self.makeTextEntry(text = "/usr/local/cantera/data/completeMechanism_gas.cti", editable = False)
        gasFileButton.connect("clicked",self.openFileSelector,"Choose a gas mechanism file",self.gasFileEntry)
        gasFileBoxGroup = self.makeBoxGroup(hargs = (gasFileLabel,self.gasFileEntry,gasFileButton))
        table.attach(gasFileBoxGroup,3,4,3,4)
        solidFileLabel = self.makeLabel("Solid File Name")
        solidFileButton = self.makeButton(type = "Button", label = "...")
        self.solidFileEntry = self.makeTextEntry(text = "/usr/local/cantera/data/completeMechanism_solid.cti", editable = False)
        solidFileButton.connect("clicked", self.openFileSelector,"Choose a solid mechanism file",self.solidFileEntry)
        solidFileBoxGroup = self.makeBoxGroup(hargs = (solidFileLabel,self.solidFileEntry,solidFileButton))
        table.attach(solidFileBoxGroup,3,4,4,5)
        
        #DATA STORAGE AND RETRIEVAL
        title = self.makeLabel("Data Storage and Retrieval")
        table.attach(title,1,3,0,1)
        #USE OLD DATA FOR PLOTTING/DATA STORAGE
        useStoredDataLabel = self.makeLabel("Use Data From a \nPrevious Simulation? ")
        self.useStoredData_On = self.makeButton(group = None, type = "RadioButton",label = "Yes")
        useStoredData_Off = self.makeButton(group = self.useStoredData_On,type = "RadioButton",label = "No")
        useStoredData_Off.set_active(1)
        useStoredDataBoxGroup = self.makeBoxGroup(hargs = (useStoredDataLabel,),vargs = (self.useStoredData_On,useStoredData_Off))
        table.attach(useStoredDataBoxGroup,1,2,1,2,xoptions = gtk.FILL)
        storedDataLabel = self.makeLabel("Profile Data File Name          ")
        storedDataButton = self.makeButton(type = "Button", label = "...")
        self.storedDataEntry = self.makeTextEntry(text = "/home/jordan/code/python/profileResults.csv", editable = False)
        storedDataButton.connect("clicked",self.openFileSelector,"Choose the stored data file",self.storedDataEntry)
        storedDataBoxGroup = self.makeBoxGroup(hargs = (storedDataLabel,self.storedDataEntry,storedDataButton))
        table.attach(storedDataBoxGroup,2,3,1,2,xoptions = gtk.FILL)
        dataFileLabel = self.makeLabel("Exit Data File Name              ")
        dataFileButton = self.makeButton(type = "Button", label = "...")
        self.dataFileEntry = self.makeTextEntry(text = "/home/jordan/code/python/exitResults.csv", editable = False)
        dataFileButton.connect("clicked",self.openFileSelector,"Choose the Storage File Name",self.dataFileEntry)
        dataFileBoxGroup = self.makeBoxGroup(hargs = (dataFileLabel,self.dataFileEntry,dataFileButton))
        table.attach(dataFileBoxGroup,2,3,2,3,xoptions = gtk.FILL)
        
        #USE EXPERIMENTAL DATA
        useExperimentalDataLabel = self.makeLabel("Load In Experimental \nData?")
        self.useExperimentalData_On = self.makeButton(group = None, type = "RadioButton",label = "Yes")
        useExperimentalData_Off = self.makeButton(group = self.useExperimentalData_On,type = "RadioButton",label = "No")
        useExperimentalData_Off.set_active(1)
        useExperimentalDataBoxGroup = self.makeBoxGroup(hargs = (useExperimentalDataLabel,),vargs = (self.useExperimentalData_On,useExperimentalData_Off))
        table.attach(useExperimentalDataBoxGroup,1,2,3,4,xoptions = gtk.FILL)
        experimentalDataLabel = self.makeLabel("Experimental Data File Name")
        experimentalDataButton = self.makeButton(type = "Button", label = "...")
        self.experimentalDataEntry = self.makeTextEntry(text = "", editable = False)
        experimentalDataButton.connect("clicked",self.openFileSelector,"Choose the experimental data file",self.experimentalDataEntry)
        experimentalDataBoxGroup = self.makeBoxGroup(hargs = (experimentalDataLabel,self.experimentalDataEntry,experimentalDataButton))
        table.attach(experimentalDataBoxGroup,2,3,3,4,xoptions = gtk.FILL)
        
        #USE PARAMETERS FROM SPREADSHEET
        useParameterSpreadsheetLabel = self.makeLabel('Use Parameters from \nSpreadsheet?')
        self.useParameterSpreadsheet_On = self.makeButton(group = None, type = "RadioButton", label = "Yes")
        useParameterSpreadsheet_Off = self.makeButton(group = self.useParameterSpreadsheet_On, type = "RadioButton", label = "No")
        useParameterSpreadsheet_Off.set_active(1)
        useParameterSpreadsheetBoxGroup = self.makeBoxGroup(hargs = (useParameterSpreadsheetLabel,),vargs =(self.useParameterSpreadsheet_On,useParameterSpreadsheet_Off))
        table.attach(useParameterSpreadsheetBoxGroup,1,2,4,5)
        parameterSpreadsheetLabel = self.makeLabel('Parameter File Name             ')
        parameterDataButton = self.makeButton(type = "Button",label = "...")
        self.parameterDataEntry = self.makeTextEntry(text = '/home/jordan/code/python/parameters.csv', editable = False)
        parameterDataButton.connect("clicked",self.openFileSelector, "Select the Parameter Data File",self.parameterDataEntry)
        parameterDataBoxGroup = self.makeBoxGroup(hargs = (parameterSpreadsheetLabel,self.parameterDataEntry,parameterDataButton))
        table.attach(parameterDataBoxGroup,2,3,4,5)
        
        
        
        table.set_row_spacings(30)
        table.set_col_spacings(30)
    
    def setupMCPage(self):
        table = gtk.Table(rows = 5, columns = 5,homogeneous = False)
        self.frameDict["MC"].add(table)
        table.show()
        #MC CONVERGENCE PARAMETERS
        solutionConvergenceLabel = self.makeLabel("Convergence Parameters")
        table.attach(solutionConvergenceLabel,0,1,0,1)
        qToleranceLabel = self.makeLabel("q Tolerance                 ")
        self.qToleranceEntry = self.makeTextEntry(text = "0.05",editable = True)
        qToleranceBoxGroup = self.makeBoxGroup(hargs = (qToleranceLabel,self.qToleranceEntry))
        table.attach(qToleranceBoxGroup,0,1,1,2)
        tempToleranceLabel = self.makeLabel("Temperature Tolerance")
        self.tempToleranceEntry = self.makeTextEntry(text = "0.05",editable = True)
        tempToleranceBoxGroup = self.makeBoxGroup(hargs = (tempToleranceLabel,self.tempToleranceEntry))
        table.attach(tempToleranceBoxGroup,0,1,2,3)
        #MC Parameteres
        MCParametersLabel = self.makeLabel("Monte Carlo Parameters")
        table.attach(MCParametersLabel,1,2,0,1)
        zSlicesLabel = self.makeLabel("Z Slices          ")
        self.zSlicesEntry = self.makeTextEntry(text = "15", editable = True)
        zSlicesBoxGroup = self.makeBoxGroup(hargs = (zSlicesLabel,self.zSlicesEntry))
        phiSlicesLabel = self.makeLabel("Phi Slices        ")
        self.phiSlicesEntry = self.makeTextEntry(text = "1", editable = True)
        phiSlicesBoxGroup = self.makeBoxGroup(hargs = (phiSlicesLabel,self.phiSlicesEntry))
        radialSlicesLabel = self.makeLabel("Radial Slices   ")
        self.radialSlicesEntry = self.makeTextEntry(text = "1", editable = True)
        radialSlicesBoxGroup = self.makeBoxGroup(hargs = (radialSlicesLabel,self.radialSlicesEntry))
        bufferFractionLabel = self.makeLabel("Buffer Fraction")
        self.bufferFractionEntry = self.makeTextEntry(text = "0.2", editable = True)
        bufferFractionBoxGroup = self.makeBoxGroup(hargs = (bufferFractionLabel,self.bufferFractionEntry))
        table.attach(zSlicesBoxGroup,1,2,1,2)
        table.attach(phiSlicesBoxGroup,1,2,2,3)
        table.attach(radialSlicesBoxGroup,1,2,3,4)
        table.attach(bufferFractionBoxGroup,1,2,4,5)
        #MONTE CARLO FILE NAMES
        MCFileNamesLabel = self.makeLabel("File Names")
        table.attach(MCFileNamesLabel,2,3,0,1)
        MCInputFileLabel = self.makeLabel("MC Input File Name   ")
        self.MCInputFileEntry = self.makeTextEntry(text = "MCInput.csv", editable = True)
        MCInputFileBoxGroup = self.makeBoxGroup(hargs = (MCInputFileLabel,self.MCInputFileEntry))
        MCOutputFileLabel = self.makeLabel("MC Output File Name")
        self.MCOutputFileEntry = self.makeTextEntry(text = "MCOutput.csv", editable = True)
        MCOutputFileBoxGroup = self.makeBoxGroup(hargs = (MCOutputFileLabel,self.MCOutputFileEntry))
        MCFilePathLabel = self.makeLabel("MC File Path             ")
        MCFilePathButton = self.makeButton(type = "Button", label = "...")
        self.MCFilePathEntry = self.makeTextEntry(text = "/home/jordan/code/c++/radiation/MCTest", editable = False)
        MCFilePathButton.connect("clicked",self.openFileSelector,"Select the MC File",self.MCFilePathEntry)
        MCFilePathBoxGroup = self.makeBoxGroup(hargs = (MCFilePathLabel,self.MCFilePathEntry,MCFilePathButton))
        table.attach(MCInputFileBoxGroup,2,3,1,2)
        table.attach(MCOutputFileBoxGroup,2,3,2,3)
        table.attach(MCFilePathBoxGroup,2,3,3,4)
 
    def cellcombo_edited(self, cellrenderer, path, new_text, data):
        model, column = data
        iter = model.get_iter(path)
        model.set_value(iter, column, new_text)
 
    def edited_cb(self, cellrenderer, path, new_text, data):
        model, column = data 
        model[path][column] = new_text
        return

    def cell_toggled(self, cellrenderer, path, data):
        """
        Sets the toggled state on the toggle button to true or false.
        """
        model, column = data
        model[path][column] = not model[path][column]
        return

    def addRow(self,cellrenderer,path,data):
        model, row = data
        model.append(row)
    
    def removeRow(self,cellrenderer,path, model):
        iter = model.get_iter(path)
        model.remove(iter)
        
    def setupSimulatorPage(self):
        table = gtk.Table(rows = 5, columns = 5,homogeneous = False)
        self.frameDict["Simulator"].add(table)
        table.show()
        
        #Integration Parameters
        title = self.makeLabel('Integration Parameters')
        atolLabel = self.makeLabel('Absolute Tolerance')
        self.atolEntry = self.makeTextEntry('1E-12',editable = True)
        rtolLabel = self.makeLabel('Relative Tolerance')
        self.rtolEntry = self.makeTextEntry('1E-12',editable = True)
        stepsLabel = self.makeLabel('Integration Steps')
        self.stepsEntry = self.makeTextEntry('100',editable = True)
        boxgroup1 = self.makeBoxGroup(hargs = (atolLabel,self.atolEntry))
        boxgroup2 = self.makeBoxGroup(hargs = (rtolLabel,self.rtolEntry))
        boxgroup3 = self.makeBoxGroup(hargs = (stepsLabel,self.stepsEntry))
        integrationBoxGroup = self.makeBoxGroup(vargs = (title,boxgroup1,boxgroup2,boxgroup3))
        table.attach(integrationBoxGroup,0,1,1,2)
        
        
        #Tree Model
        #Create a model for the treeview.
        self.parameterModel = gtk.ListStore(str,str)
        self.solidFlowRateModel = gtk.ListStore(str,str,'gboolean','gboolean')
        self.gasFlowRateModel = gtk.ListStore(str,str,str,'gboolean','gboolean','gboolean')
        self.parameterModel.append(['Pressures [psig]','10',])
        #self.parameterModel.append(['Gas Temperatures [C]','800'])
        self.parameterModel.append(['Solid Temperatures [C]','25'])

        #self.parameterModel.append(['Total Solids Flow Rates',''])
        self.parameterModel.append(['Particle Diameters [microns]','100'])
        self.parameterModel.append(['Effectiveness Factor [phi]','0.5'])
        self.parameterModel.append(['Gas Thermal Conductivity [W/m - K]','0.050'])
        self.parameterModel.append(['Solid Density [kg/m^3]','1000'])
        self.parameterModel.append(['Total Solids Flow [lb/hr]',''])
        
        self.solidFlowRateModel.append(['CELL','0.370',True,True])
        self.solidFlowRateModel.append(['HCE','0.287',True,True])
        self.solidFlowRateModel.append(['LIGC','0.0528',True,True])
        self.solidFlowRateModel.append(['LIGH','0.0528',True,True])
        self.solidFlowRateModel.append(['LIGO','0.0528',True,True])
        self.solidFlowRateModel.append(['ASH','0.185',True,True])
        self.solidFlowRateModel.append(['Select a Solid','',True,True])
        
        self.gasFlowRateModel.append(['O2','25','800',False,True,True])
        self.gasFlowRateModel.append(['H2','5','800',False,True,True])
        self.gasFlowRateModel.append(['Select a Gas','','',False,True,True])
        #create a combobox columns of solids and gas species
        solidSpeciesList = ('CELL','HCE','LIGC','LIGH','LIGO','ASH')
        gasSpeciesList = ('N2','H2O','CH4','O2','H2','CO2')
        speciesList = ('CELL','HCE','LIGC','LIGH','LIGO','ASH','N2','H2O','CH4','O2','H2')
        self.solidSpeciesModel = gtk.ListStore(str)
        self.gasSpeciesModel = gtk.ListStore(str)
        self.speciesModel = gtk.ListStore(str)
        for species in solidSpeciesList:
            self.solidSpeciesModel.append([species])
        for species in gasSpeciesList:
            self.gasSpeciesModel.append([species])
        for species in speciesList:
            self.speciesModel.append([species])
            
       
                
                
        #Create cell renderers.  
        cell = gtk.CellRendererText()
        cell2 = gtk.CellRendererText()
        cell2.set_property('editable', True)
        cell2.connect('edited', self.edited_cb, (self.parameterModel,1))
        
        cellcombo3 = gtk.CellRendererCombo()
        cellcombo3.set_property("text-column", 0)      #The column in speciesModel to grab text from.  (The model only has 1 column)
        cellcombo3.set_property("editable", True)      #Make the cell editable
        cellcombo3.set_property("has-entry", False)  
        cellcombo3.set_property("model", self.solidSpeciesModel) #The model of values to use.  In this case, speciesModel
        cellcombo3.connect('edited',self.cellcombo_edited,(self.solidFlowRateModel,0))
        cell4 = gtk.CellRendererText()
        cell4.set_property('editable',True)
        cell4.connect('edited', self.edited_cb, (self.solidFlowRateModel,1))
        cell5 = gtk.CellRendererToggle()
        cell5.set_property('activatable', True)
        cell5.set_property('xalign',0.5)
        cell5.connect("toggled",self.addRow,(self.solidFlowRateModel,['Select a Solid','',True,True]))
        cell6 = gtk.CellRendererToggle()
        cell6.set_property('xalign', 0.5)
        cell6.connect('toggled',self.removeRow,self.solidFlowRateModel)
        
        cellcombo7 = gtk.CellRendererCombo()
        cellcombo7.set_property("text-column", 0)      #The column in speciesModel to grab text from.  (The model only has 1 column)
        cellcombo7.set_property("editable", True)      #Make the cell editable
        cellcombo7.set_property("has-entry", False)  
        cellcombo7.set_property("model", self.gasSpeciesModel) #The model of values to use.  In this case, speciesModel
        cellcombo7.connect('edited',self.cellcombo_edited,(self.gasFlowRateModel,0))
        
        cell8 = gtk.CellRendererText()
        cell8.set_property('editable',True)
        cell8.connect('edited', self.edited_cb, (self.gasFlowRateModel,1))
        
        cell9 = gtk.CellRendererText()
        cell9.set_property('editable',True)
        cell9.connect('edited', self.edited_cb, (self.gasFlowRateModel,2))
        
        cell10 = gtk.CellRendererToggle()
        cell10.set_property('activatable', True)
        cell10.connect('toggled', self.cell_toggled, (self.gasFlowRateModel,3))
        cell11 = gtk.CellRendererToggle()
        cell11.set_property('activatable', True)
        cell11.set_property('xalign',0.1)
        cell11.connect("toggled",self.addRow,(self.gasFlowRateModel,['Select a Gas','','',False,True,True]))
        cell12 = gtk.CellRendererToggle()
        cell12.set_property('xalign', 0.1)
        cell12.connect('toggled',self.removeRow,self.gasFlowRateModel)
        
        #Create the treeview and set the model as the ListStore created above
        treeview1 = gtk.TreeView(self.parameterModel)
        treeview2 = gtk.TreeView(self.solidFlowRateModel)
        treeview3 = gtk.TreeView(self.gasFlowRateModel)
        
        #Create columns.  First argument is column name, 2nd argument is the cell renderer, 
        #text argument is the column to take the text from in model (the ListStore).
        #NOTE: Packing multiple cell renderers into a column would require using column.pack_start(cell1,cell2,etc..)
        #NOTE2: Other renderer attributes can be set with column.set_attributes(cellRenderer, text, etc...)
        column = gtk.TreeViewColumn("Run Parameters",cell,text = 0)
        column2 = gtk.TreeViewColumn("Values",cell2,text = 1)
        column3 = gtk.TreeViewColumn("Solids",cellcombo3,text = 0)
        column4 = gtk.TreeViewColumn("Rates [lb/hr]",cell4,text = 1)
        column5 = gtk.TreeViewColumn("Add Row", cell5, radio = 2)
        column6 = gtk.TreeViewColumn("Remove Row", cell6, radio = 2,activatable = 3)
        
        column7 = gtk.TreeViewColumn("Gases",cellcombo7,text = 0)
        column8 = gtk.TreeViewColumn("Rates [lb/hr]",cell8,text = 1)
        column9 = gtk.TreeViewColumn("Temperature [C]",cell9,text = 2)
        column10 = gtk.TreeViewColumn("Normalize?", cell10, active = 3)
        column11 = gtk.TreeViewColumn("Add Row", cell11, radio = 4)
        column12 = gtk.TreeViewColumn("Remove Row", cell12, radio = 4,activatable = 5)


        #Add the columns to the treeview
        treeview1.append_column(column)
        treeview1.append_column(column2)
        treeview2.append_column(column3)
        treeview2.append_column(column4)
        treeview2.append_column(column5)
        treeview2.append_column(column6)
        treeview3.append_column(column7)
        treeview3.append_column(column8)
        treeview3.append_column(column9)
        treeview3.append_column(column10)
        treeview3.append_column(column11)
        treeview3.append_column(column12)

        def totalSolidFlowSet_cb(cellrenderer, path, new_text, data):
            model, columns = data
            iter = model.get_iter(path)
            values = (list(model.get(iter,*columns)))
            if values [0] == 'Total Solids Flow [lb/hr]' and values[1] != '':
                column4.set_title("Mass Fractions")
            elif values [0] == 'Total Solids Flow [lb/hr]' and values[1] == '':
                column4.set_title("Rates [lb/hr]")
                   
        cell2.connect('edited',totalSolidFlowSet_cb,(self.parameterModel,(0,1)))  #This callback checks if a value for total solid flow has been entered.
        
        
        table.attach(treeview1,0,1,0,1)
        table.attach(treeview2,1,2,0,1)
        table.attach(treeview3,1,2,1,2)
        self.window.show_all()
        
#        parameterTextEntry = self.makeTextEntry(text = "", editable = True)
#        parameterTextButton = self.makeButton(type = "Button", stock = gtk.STOCK_ADD)
#        parameterSetBoxGroup = self.makeBoxGroup(hargs = (parameterTextEntry,parameterTextButton))
#        table.attach(parameterSetBoxGroup,1,2,2,3)
        table.set_col_spacings(30)
        
    def setupTubePage(self):
        table = gtk.Table(rows = 5, columns = 5,homogeneous = False)
        self.frameDict["Tube"].add(table)
        table.show()
        #Tree Model
        #attributes = ['Pressure','Gas Temperature','Wall Temperature','Particle Diameter']
        self.tubeParametersModel = gtk.ListStore(str,str)
        self.tubeParametersModel.append(['Diameter [in]','3.5'])
        self.tubeParametersModel.append(['Length [ft]','7'])
        self.tubeParametersModel.append(['Wall Temperatures [C]','1400'])
        self.tubeParametersModel.append(['Wall Fluxes [W/m^2]','60'])
        self.tubeParametersModel.append(['Emissivity','0.9'])
        self.tubeParametersModel.append(['Overall Heat Transfer Coefficient [W/m^2 - K]','150'])
        #Create the treeview and set the model as the ListStore
        treeview = gtk.TreeView(self.tubeParametersModel)
        #Create a column for the treeview
        column = gtk.TreeViewColumn("Tube Parameters")
        column2 = gtk.TreeViewColumn("Values")
        #Create a cell renderer and add it to the column with the text attribute
        cell = gtk.CellRendererText()
        cell2 = gtk.CellRendererText()
        column.pack_start(cell)
        column2.pack_start(cell2)
        column.set_attributes(cell,text = 0)    #Sets text to uneditable
        #column2.set_attributes(cell2, text = 1)
        cell2.set_property('editable', True)
        #cell.set_property('editable',True)
        column2.set_attributes(cell2, text = 1) #Sets text to editable
        cell2.connect('edited', self.edited_cb, (self.tubeParametersModel,1))
    
        #Add the column to the treeview
        treeview.append_column(column)
        treeview.append_column(column2)
        treeview.get_selection().set_mode(gtk.SELECTION_MULTIPLE)  #Allow multiple selection

        #Connect the change signal to the print function above
        table.attach(treeview,0,1,0,5)
        self.window.show_all()
        
    def setupPlottingPage(self):
        table = gtk.Table(rows = 5, columns = 5,homogeneous = False)
        table.set_col_spacings(30)
        self.frameDict["Plotting"].add(table)
        table.show()

        self.profileModel = gtk.ListStore(str,'gboolean')
        self.comparisonModel = gtk.ListStore(str,str,str,str,'gboolean','gboolean','gboolean')
        self.profileModel.append(['Solids Temperature [C]',False])  
        self.profileModel.append(['Gas Temperature [C]',False])
        self.profileModel.append(['Phase Temperatures [C]', True]) #checked True to start
        self.profileModel.append(['Wall Temperature [C]',False])
        self.profileModel.append(['Carbon Conversion',False])
        self.profileModel.append(['Concentrations [Mole Fractions]',False])
        self.profileModel.append(['Molar Rates [mol/min]',False])
        self.profileModel.append(['Wall Flux [kW/m^2]',False])
        self.comparisonModel.append(['X Data','Y Data','','',False,False,False])
        
        attributeList = ('Length [ft]','Pressure [psig]','Exit Gas Temperature [C]','Exit Carbon Conversion','Heat Duty [W]','Solids Total Flow [lb/hr]',
                         'Molar Rate [mol/min]','Concentration [Mole Fraction]' )
        comboModel = gtk.ListStore(str)
        for species in attributeList:
            comboModel.append([species])
            
        #Create cell renderers.  
        cellcombo3 = gtk.CellRendererCombo()
        cellcombo3.set_property("text-column", 0)      #The column in comboModel to grab text from.  (The model only has 1 column)
        cellcombo3.set_property("editable", True)      #Make the cell editable
        cellcombo3.set_property("has-entry", False)  
        cellcombo3.set_property("model", comboModel) 
        cellcombo3.connect('edited',self.cellcombo_edited,(self.comparisonModel,0))
        
        cellcombo4 = gtk.CellRendererCombo()
        cellcombo4.set_property("text-column", 0)      #The column in speciesModel to grab text from.  (The model only has 1 column)
        cellcombo4.set_property("editable", True)      #Make the cell editable
        cellcombo4.set_property("has-entry", False)  
        cellcombo4.set_property("model", comboModel)   #The model of values to use.  In this case, speciesModel
        cellcombo4.connect('edited',self.cellcombo_edited,(self.comparisonModel,1))   #model, column number

        cellcombo5 = (gtk.CellRendererCombo())
        cellcombo5.set_property("text-column",0)
        cellcombo5.set_property("editable", False)
        cellcombo5.set_property("has-entry",False)
        cellcombo5.set_property("model",self.speciesModel)
        cellcombo5.connect('edited',self.cellcombo_edited,(self.comparisonModel,2))
        cellcombo6 = (gtk.CellRendererCombo())
        cellcombo6.set_property("text-column",0)
        cellcombo6.set_property("editable", False)
        cellcombo6.set_property("has-entry",False)
        cellcombo6.set_property("model",self.speciesModel)
        cellcombo6.connect('edited',self.cellcombo_edited,(self.comparisonModel,3))
        
        cell7 = gtk.CellRendererToggle()
        cell7.set_property('activatable', True)
        cell7.set_property('radio', True)
        cell7.set_property('xalign',0.1)
        cell7.connect("toggled",self.addRow,(self.comparisonModel,['X Data','Y Data','','',False,False,True]))
        cell8 = gtk.CellRendererToggle()
        cell8.set_property('xalign', 0.1)
        cell8.set_property('radio', True)
        cell8.connect('toggled',self.removeRow,self.comparisonModel)
        
        
        def conditional_edit(cellrenderer, path,new_text, data):
            model, column1, column2 = data
            if new_text == 'Molar Rate [mol/min]' or new_text == 'Concentration [Mole Fraction]':
                model[path][column1] = "Select a species"
                model[path][column2] = True
                
        cellcombo3.connect('edited',conditional_edit,(self.comparisonModel,2,4)) 
        cellcombo4.connect('edited',conditional_edit,(self.comparisonModel,3,5))  

        cell = gtk.CellRendererText()
        cell2 = gtk.CellRendererToggle()
        cell2.set_property('activatable', True)
        cell2.connect('toggled', self.cell_toggled, (self.profileModel,1))
        
        #Create treeview
        treeview = gtk.TreeView(self.profileModel)
        treeview2 = gtk.TreeView(self.comparisonModel)
        #Create Columns
        column = gtk.TreeViewColumn("Profile Name",cell,text = 0)
        column2 = gtk.TreeViewColumn("Check Profile",cell2, active = 1)
        column3 = gtk.TreeViewColumn("X Data",cellcombo3,text = 0)
        column4 = gtk.TreeViewColumn("Y Data",cellcombo4,text = 1)
        column5 = gtk.TreeViewColumn("X Species", cellcombo5, text = 2, editable = 4)
        column6 = gtk.TreeViewColumn("Y Species", cellcombo6,text = 3, editable = 5)
        column7 = gtk.TreeViewColumn("Add Row",cell7)
        column8 = gtk.TreeViewColumn("Remove Row",cell8,activatable = 6 )
        
        #Append columns to treeview
        treeview.append_column(column)
        treeview.append_column(column2)
        treeview2.append_column(column3)
        treeview2.append_column(column4)
        treeview2.append_column(column5)
        treeview2.append_column(column6)
        treeview2.append_column(column7)
        treeview2.append_column(column8)
        table.attach(treeview,0,1,0,4)
        table.attach(treeview2,1,2,0,4)
        
        #Save Figures
        saveFigureLabel = self.makeLabel("Save Figure Directory   ")
        self.saveFigureEntry = self.makeTextEntry(text = "/home/jordan/code/python/Figures/", editable = True)
        saveFigureBoxGroup = self.makeBoxGroup(hargs = (saveFigureLabel,self.saveFigureEntry))
        table.attach(saveFigureBoxGroup,0,2,4,5)
        self.window.show_all()

        

def main():
    gtk.main()
    return 0

if __name__ == "__main__":
    GUI = GasifierGUI()
    GUI.createNotebook()
    GUI.setupGeneralPage()
    GUI.setupMCPage()
    GUI.setupSimulatorPage()
    GUI.setupTubePage()
    GUI.setupPlottingPage()

    main()