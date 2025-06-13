#!/bin/bash
###########################################################
# Create DVCS filelist from S3 folder                     #
###########################################################

###########################################################
# Help message                                            #
###########################################################
Help()
{
    # Display help
    echo "Create list of ePIC simulation campaign files by campaign date and beam settings."
    echo "NOTE: This script MUST be run within eic-shell to work."
    echo
    echo "Syntax:  getFilelist.sh [-b|c|e|h]"
    echo "options:"
    echo "b     Beam setting (hiAcc, hiDiv)."
    echo "c     Campaign to look at (23.08.0, 23.10.1, etc.)."
    echo "e     Beam energy configuration (\"5x41\", \"10x100\", \"18x275\" (physics beams), \"10x130\" (early science))."
    echo
    echo "Help message is run if no options are provided."
}

###########################################################
# Main script                                             #
###########################################################

# Print help message if no arguments are provided
if [[ -z $1 ]];
then
    Help
    exit
fi

# Set initial variables
BeamSet="X"
Camp="X"
Energy="X"


# Parse over all options
while getopts h:b:c:e: option; do
    case $option in
	h) # Print help message
	    Help
	    exit;;
	b) # Set beam setting
	    BeamSet=$OPTARG;;
	c) # Set campaign to use
	    Camp=$OPTARG;;
	e) # Set beam energy
	    Energy=$OPTARG;;
       \?) # Invalid options
	    echo "Error: Invalid flag"
	    exit;;
    esac
done

# Check for valid inputs: beam setting
if [ $BeamSet == "X" ]
then
    echo "Must declare beam setting."
    exit
elif [ $BeamSet != "hiDiv" ] && [ $BeamSet != "hiAcc" ]
then
    echo "Invalid beam setting."
    exit
fi
# Check for input: campaign
if [ $Camp == 'X' ]
then
    echo "Must declare campaign to use."
    exit
fi
# Check for valid inputs: beam energy
if [ $Energy == "X" ]
then
    echo "Must declare beam energy."
    exit
elif [ $Energy != "5x41" ] && [ $Energy != "10x100" ] && [ $Energy != "18x275" ] && [ $Energy != "10x130" ]
then
    echo "Invalid beam energy."
    exit
fi
# Only high acceptance detector setting available for early science beams - reject hiDiv & replace w/ hiAcc
if [ $Energy == "10x130" ] && [ $BeamSet == "hiDiv" ]
then
    echo "hiDiv setting not available for early science beams"
    echo "Using hiAcc instead"
    BeamSet="hiAcc"
fi

# Sep. '24 -> Simulation files now no longer written to S3. Need to convert purely to XRootD syntax
#xrdfs root://dtn-eic.jlab.org ls /work/eic2/EPIC/RECO/$Camp/epic_craterlake/EXCLUSIVE/DVCS_ABCONV/$Energy | grep $BeamSet | sed "s,/work,root://dtn-eic.jlab.org//work," > filelists/inputFileList_ePIC_"$Camp"_"$Energy"_"$BeamSet".list

# Jan. '25 -> Campaign files moved to a different location within JLab workdisk - for campaigns 25.01.1 and beyond
#xrdfs root://dtn-eic.jlab.org ls /volatile/eic/EPIC/RECO/$Camp/epic_craterlake/EXCLUSIVE/DVCS_ABCONV/$Energy | grep $BeamSet | sed "s,/volatile,root://dtn-eic.jlab.org//volatile," > filelists/inputFileList_ePIC_"$Camp"_"$Energy"_"$BeamSet".list

# May '25 -> include option for early science files (only available from campaign 25.04.1)
if [ $Energy == "10x130" ]
then
    xrdfs root://dtn-eic.jlab.org ls /volatile/eic/EPIC/RECO/$Camp/epic_craterlake/EXCLUSIVE/DVCS_ABCONV/EpIC_v1.1.6-1.0/$Energy/q2_1_100 | sed "s,/volatile,root://dtn-eic.jlab.org//volatile," > filelists/inputFileList_ePIC_"$Camp"_"$Energy"_"$BeamSet".list
    exit
else
    xrdfs root://dtn-eic.jlab.org ls /volatile/eic/EPIC/RECO/$Camp/epic_craterlake/EXCLUSIVE/DVCS_ABCONV/$Energy | grep $BeamSet | sed "s,/volatile,root://dtn-eic.jlab.org//volatile," > filelists/inputFileList_ePIC_"$Camp"_"$Energy"_"$BeamSet".list
    exit
fi
