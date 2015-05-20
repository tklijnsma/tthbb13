#!/usr/bin/env python
"""
Thomas:

Retrieves the statistics from individual in job.stdout.gz files.

"""

########################################
# Imports
########################################

import os
import gzip


########################################
# Functions
########################################

def Create_Error_Report( workdir ):

    ########################################
    # Fill the dict
    ########################################

    n_specific_prints = 0

    dirs = os.listdir( 'work.{0}/output/'.format(workdir) )

    out_f = open( 'AllErrors-{0}.txt'.format(workdir) , 'w' )

    for jobdir in dirs:

        stdout_filename = 'work.{0}/output/{1}/job.stderr.gz'.format(workdir,jobdir)

        if not os.path.isfile( stdout_filename ):
            #print '{0} does not exist (yet)'.format( stdout_filename )
            continue

        out_f.write( '\nERROR IN {0}\n=================\n'.format(jobdir) )

        f = gzip.open(stdout_filename, 'rb')
        out_f.write( f.read() )
        f.close()

    out_f.close()



########################################
# Main
########################################

def main():

    Create_Error_Report( 'sig' )
    #Create_Error_Report( 'bkg' )
    




########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
