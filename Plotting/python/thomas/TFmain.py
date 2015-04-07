#!/usr/bin/env python
"""
Thomas:

This is TFmain.py. This program runs the full chain of creating a matrix of transfer functions from a customized .root file. To create the customized .root file, see 'outputtree.py'

"""

########################################
# Imports
########################################

import pickle
import time
import datetime
import shutil

from config import Make_config
from readtree import Single_Bin_Fits
from fit_across_bins import Fit_Across_Bins
from draw_hists_and_fits import Draw_Hists_and_Fits

########################################
# Main
########################################

def main():

    print 'This is TFmain.py. This program runs the full chain of creating a matrix of transfer functions from a customized .root file. To create the customized .root file, see \'outputtree.py\'\n'

    print '\n###################################'
    print '# TFmain: START OF CONFIGURATION'
    print '###################################'
    # Create the config.dat file  
    Make_config()
    print '###################################'
    print '# TFmain: END OF CONFIGURATION'
    print '###################################'

    print '\n###################################'
    print '# TFmain: START OF SINGLE BIN FITTING'
    print '###################################'
    # Fit the individual E-eta bins
    Single_Bin_Fits()
    print '###################################'
    print '# TFmain: END OF SINGLE BIN FITTING'
    print '###################################'

    print '\n###################################'
    print '# TFmain: START OF FITTING ACROSS BINS'
    print '###################################'
    # Fit across the E-bins
    Fit_Across_Bins()
    print '###################################'
    print '# TFmain: END OF FITTING ACROSS BINS'
    print '###################################'

    pickle_f = open( 'config.dat', 'rb' )
    config = pickle.load( pickle_f )
    pickle_f.close()

    if config['Draw_ABF_in_single_bins']:
        print '\n###################################'
        print '# TFmain: START OF DRAWING ABF RESULTS'
        print '###################################'
        # Drawing the single bin histograms and the ABF results
        Draw_Hists_and_Fits()
        print '###################################'
        print '# TFmain: END OF DRAWING ABF RESULTS'
        print '###################################'

    print ''
    print config['info']

    ts = time.time()
    config['enddate'] = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    print ' Analysis end time:              {0}\n'.format( config['enddate'] )


    # Save latest version of config, and copy to the output folder
    pickle_f = open( 'config.dat', 'wb' )
    pickle.dump( config, pickle_f )
    pickle_f.close()
    shutil.copyfile( 'config.dat', '{0}/config.dat'.format( config['outputdir'] ) )

########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
