"""
Title:      Mapping Code Wrapper for Batch Model Processing
Author:     James R. Beattie
Created:    May, 2019

Details:

This Python script is a wrapper for the Mappings code located at: https://miocene.anu.edu.au/mappings/

The aim is to call the Mappings code, populate it with the same values for each of the models, changing
either the density or the electron temperature, depending on the users command line arguments.

Currently the radiation source is set to the default.

"""

# Dependencies
#####################################################################################################

import subprocess               # for running the mapping code
from subprocess import PIPE     # for accessing the mapping output
import keyboard                 # using the keyboard, requires python to have permissions
import time                     # for waiting between events, otherwise Python is too fast for mapping
import argparse                 # command line arguments

# Command Line Arguments
#####################################################################################################

ap 			= argparse.ArgumentParser(description = 'Input arguments')
ap.add_argument('-numOfModels',required=False,default=1,help='the number of models that will be run.',type=int)
ap.add_argument('-initTemp',required=False,default=1e4,help='the initial temperature for the models',type=float)
ap.add_argument('-finalTemp',required=False,default=None,help='the final temperature for the models',type=float)
ap.add_argument('-initDens',required=False,default=1,help='the initial density for the models',type=float)
ap.add_argument('-finalDens',required=False,default=None,help='the final density for the models',type=float)
args 		= vars(ap.parse_args())

"""
#######################################################
Example Command Line Calls (from within iPython)
#######################################################

Using default values T = 1e4, and density = 1:

> run AutoMapping -numOfModels 50

Changing the density from 2 to 10 with constant temperature, T = 1e4:

> run AutoMapping -numOfModels 50 -initDens 2  -finalDens 10

Changing the temperature from 1e5 to 1e6 with constant density = 1:

> run AutoMapping -numOfModels 50 -initTemp 1e5  -finalTemp 1e6

"""

# Call the Mappings Code
#####################################################################################################

def RunMapping(command="Map51"):
    """
    RunMapping runs the mapping code and handles all of the mapping
    events as they arrive at command line.

    INPUTS:
    command - the command line call to the mapping software.

    OUTPUTS:
    the regular outputs from the mapping software.

    """

    class IPs:
        """
        A class for storing and updating some initial parameters / variables.

        """

        LineCounter = 0     # Counts lines after a particular state.
        State1      = None  # State1 controls some events.
        State2      = 0     # State2 controls some events.
        IterState   = 0     # Iteration state keeps track of the current iteration.

        # Declare the final iteration and the initial conditions
        FinalIter           = args['numOfModels']    # The total number of models that will be run.

        # The first set of initial conditions, which take command line arguments
        InitialConditions   = [args['initTemp'], args['initDens'], 1, 0.5, -1]

        # Create a log file for the mappings output
        MappingsLog         = open("MappingsLog.txt","w")

    def TypeAndPress(condition,write,press,wait):
        """
        TypeAndPress handles keyboard typing and keyboard pressing
        events.

        INPUTS:
        condition   - the read condition coming from mapping, e.g. 'Include cosmic ray heating? (Y/N):'
                    is read then press X and write Y.
        write       - what needs to be written from the keyboard
        press       - they key that is pressed e.g. 'enter'
        wait        - how long the wrapper needs to wait between events.

        OUTPUTS:
        a completed event.

        """

        if write is not None:
            if output.strip() == condition:
                time.sleep(wait)
                keyboard.write(write)
                time.sleep(wait)
                keyboard.press_and_release(press)
        else:
            if output.strip() == condition:
                time.sleep(wait)
                keyboard.press_and_release(press)

    def WaitAndEnterAndWait(InitialCondition):
        """
        A wait and press event handler for the initial conditions part of
        the Mappings code.

        INPUTS:
        InitialCondition - the initial conditions for each of the parameters
                            for the mappings code.

        OUTPUTS:
        a completed initial conditions event.

        """

        keyboard.write("{}".format(InitialCondition))
        time.sleep(0.01)
        keyboard.press_and_release('enter')
        time.sleep(0.01)

    def WaitAndWriteAndEnter(write):
        """
        A wait, write and enter key event handler.

        INPUTS:
        write   - the string input for what to write into the mappings code

        OUTPUTS:
        a completed wait, write and press enter event.

        """

        time.sleep(0.01)
        keyboard.write(write)
        keyboard.press_and_release('enter')

    def InitialiseParameters():
        """
        An initialise all parameters event handler for reseting the parameters
        after each finished model.

        INPUTS:
        None

        OUTPUTS:
        initialise the parameters event complete.

        """

        IPs.LineCounter = 0     # Counts lines after a particular state.
        IPs.State1      = None  # State1 controls some events.
        IPs.State2      = 0     # State2 controls some events.

    def UpdateICs():
        """
        An update the initial conditions event handler that picks up on which experiment to run based
        on the command line arguments.

        INPUTS:
        None

        OUTPUTS:
        updated initial conditions for the next set of models to run.

        """

        # increment the temperature for the temperature experiment
        if args['finalTemp'] is not None:
            IPs.InitialConditions[0] += (args['finalTemp'] - args['initTemp'])/(args['numOfModels']-1)
            PrintAndLogOutput("Temp")

        # increment the density for the density experiment.
        if args['finalDens'] is not None:
            IPs.InitialConditions[1] += (args['finalDens'] - args['initDens'])/(args['numOfModels']-1)
            PrintAndLogOutput("Dens")

    def PrintAndLogOutput(type):
        """
        Prints out the change in value of parameters or model and the logs it into the mappingslog variable.
        The output is MappingsLog.txt

        INPUT:
        type    - the type of input, either temp, density or the model itself

        OUPUT:
        a completed logging event of a change.
        """

        if type == "Temp":
            print("\n\n ################################################ \n\n")
            print(" The temperature has just been updated to {}".format(IPs.InitialConditions[0]))
            print("\n\n ################################################ \n\n")

            IPs.MappingsLog.write("\n\n ################################################ \n\n")
            IPs.MappingsLog.write(" The temperature has just been updated to {}".format(IPs.InitialConditions[0]))
            IPs.MappingsLog.write("\n\n ################################################ \n\n")

        elif type == "Dens":
            print("\n\n ################################################ \n\n")
            print(" The density has just been updated to {}".format(IPs.InitialConditions[1]))
            print("\n\n ################################################ \n\n")

            IPs.MappingsLog.write("\n\n ################################################ \n\n")
            IPs.MappingsLog.write("The density has just been updated to {}".format(IPs.InitialConditions[1]))
            IPs.MappingsLog.write("\n\n ################################################ \n\n")

        elif type == "Model":
            print("\n\n ################################################ \n\n")
            print(" Moving to model {} \n".format(IPs.IterState))
            print("\n\n ################################################ \n\n")

            IPs.MappingsLog.write("\n\n ################################################ \n\n")
            IPs.MappingsLog.write(" The model number is {}".format(IPs.IterState))
            IPs.MappingsLog.write("\n\n ################################################ \n\n")


    # Call the mappings code using a subprocess
    proc = subprocess.Popen(command,stdout=PIPE)

    # Print lines of the mapping code to command line until there is nothing to print
    # or if the process has been terminated .poll()
    while True:
        output = proc.stdout.readline()
        IPs.MappingsLog.write(output)

        # if nothing is being output or termination has occured, terminate
        if output == '' and proc.poll() is not None:
            break

        # if there exists an output
        if output:
            print(output.strip())

            # This is a random enter event that needs to be done to continue mappings
            if IPs.LineCounter == 3 and IPs.State1 == "EnterToStart" and IPs.State2 == 1:
                keyboard.press_and_release('enter')
                IPs.LineCounter = 0
                IPs.State1      = None

            # This is an exit condition event
            if output.strip() == "E,X,Q :  Exit":

                # if the current iterstate is less than the final iterstate, keep going
                # and reset the parameter states.
                if IPs.IterState < IPs.FinalIter:
                    # use the SINGLE SLAB model
                    WaitAndWriteAndEnter("SS")

                    # Reset the parameters
                    if IPs.IterState >= 1:
                        InitialiseParameters()

                        # update the initial coniditions for the next model
                        UpdateICs()

                    # Increase the iterstate by one.
                    IPs.IterState += 1
                    PrintAndLogOutput("Model")
                else:
                    # EXIT the mappings code after the iterstate > final iterstate
                    WaitAndWriteAndEnter("Exit")

            # This picks up one of the events where an enter needs to be pressed
            if output.strip() == "[ mu_neu:  1.2584       mu_ion: 0.60364      mu_h:  1.3668     ]":
                if IPs.State2 == 0:
                    IPs.LineCounter     += 1
                    IPs.State1           = "EnterToStart"
                    IPs.State2          += 1
                else:
                    IPs.LineCounter     = 0
                    IPs.State1          = None
                    IPs.State2          = 0

            # This presses N for the cosmic ray heating event
            if output.strip() == "Include cosmic ray heating? (Y/N):" and IPs.State2 == 0:
                IPs.State2 +=1
                WaitAndWriteAndEnter("N")

            # This begins the initial conditions event
            if output.strip() == "(Case A-B (H) taken as 0(A) <-> 1 (B), <0 auto)":
                IPs.LineCounter     +=1
                IPs.State1          = "InitialConditions"

            # This handles the initial conditions event
            if IPs.LineCounter == 2 and IPs.State1 == "InitialConditions":
                time.sleep(0.05)
                WaitAndEnterAndWait(IPs.InitialConditions[0])
                WaitAndEnterAndWait(IPs.InitialConditions[1])
                WaitAndEnterAndWait(IPs.InitialConditions[2])
                WaitAndEnterAndWait(IPs.InitialConditions[3])
                WaitAndEnterAndWait(IPs.InitialConditions[4])
                keyboard.press_and_release('enter')

            # Here are all of the rest of the events that can be handled as functionally the same event
            TypeAndPress("Change abundances (y/N) :","No","enter",0)
            TypeAndPress("Change Abundance Offsets (y/N)? :","No","enter",0)
            TypeAndPress("Dust is currently disabled.","No","enter",0)
            TypeAndPress("F :  Time dependent ionisation and temperature.","A","enter",0)
            TypeAndPress("X : eXit with current balance","X","enter",0)
            TypeAndPress("X  :   eXit with current source","X","enter",0)
            TypeAndPress("X   0.00       0.00   0.00   0.00   0.00   0.00      0.00      0.00",None,"enter",0)
            TypeAndPress("Spectrum printout required? (y/n)","y","enter",0)
            TypeAndPress("Give a name/code for this run:","Experiment_{}".format(IPs.IterState),"enter",0)

            # Start counting after the line counter has been set off
            if IPs.LineCounter >= 1:
                IPs.LineCounter += 1

            # Reset the line counter after it gets to a large number,
            # to resuse later.
            elif IPs.LineCounter > 10:
                IPs.LineCounter = 0

            # Reset one of the states if it gets larger than 2
            if IPs.State2 == 2:
                IPs.State2 = 0

            # wait a small amount of time per line print from mappings
            time.sleep(0.01)

    MapLine = proc.poll()
    IPs.MappingsLog.close()

    return MapLine

RunMapping()
