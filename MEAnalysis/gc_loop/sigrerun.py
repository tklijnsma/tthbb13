#!/usr/bin/env python
"""
Thomas:

"""

########################################
# Imports
########################################

import subprocess
import os
import shutil
import re

from CreateStatisticsReport import Create_sig_Report

from time import sleep


########################################
# Functions
########################################

def Change_cfg_heppy( gc_loop_dir, Cut_criteria ):

    fn_cfg_temp = '{0}/MEAnalysis_cfg_heppy_template.py'.format(gc_loop_dir)
    fn_cfg = '{0}/MEAnalysis_cfg_heppy.py'.format(gc_loop_dir)

    f_cfg_temp = open( fn_cfg_temp, 'r' )
    cfg_temp = f_cfg_temp.read()
    f_cfg_temp.close()

#    Cut_criteria = \
#        """
#            ( 'pt'  , '>', '200.0' ),
#            ( 'mass', '>', '120.0' ),
#            ( 'mass', '<', '220.0' ),
#            ( 'fW'  , '<', '0.175' )
#        """

    cfg = cfg_temp.replace( '###INSERTCUTCRITERIA###', Cut_criteria )

    f_cfg = open( fn_cfg, 'w' )
    f_cfg.write( cfg )
    f_cfg.close()

    ( dst, junk ) = os.path.split( gc_loop_dir )
    dst = dst + '/python/MEAnalysis_cfg_heppy.py'

    shutil.copy( fn_cfg, dst )



########################################
# Main
########################################

def main():

    gc_loop_dir = os.getcwd()

    gc_rel_dir = '../gc'
    os.chdir( gc_rel_dir )
    gc_dir = os.getcwd()

    print os.getcwd()

    gopy = 'grid-control/go.py'

    submitcmd   = [ gopy , 'confs/sig-psi.conf', '-q' ]
    statuscmd   = [ gopy , 'confs/sig-psi.conf', '-qs' ]
    retrievecmd = [ gopy , 'confs/sig-psi.conf', '-r' ]

    var_list = [ '80.0', '100.0', '120.0', '130.0', '135.0', '140.0', '145.0', '150.0', '155.0', '160.0', '180.0' ]


    # Start of loop
    for i_iter in range(len(var_list)):

        # Delete currently existing jobs
        subprocess.call( [ 'qdel', '-u', 'tklijnsm' ] )

        # Remove the work.sig directory
        if os.path.isdir('work.sig'):
            print 'Removing work.sig'
            shutil.rmtree('work.sig')

        # Set the cut criteria
        Cut_criteria = \
            """
                ( 'pt'  , '>', '200.0' ),
                ( 'mass', '>', {0} ),
                ( 'mass', '<', '250.0' ),
                ( 'fW'  , '<', '0.175' )
            """.format( var_list[i_iter] )

        # Write the cfg file and copy it to MEAnalysis/python/
        Change_cfg_heppy( gc_loop_dir, Cut_criteria )

        # Submit the jobs
        print 'Submitting jobs'
        subprocess.call( submitcmd, stdout=open(os.devnull, 'wb') )

        # Check status repeatedly - detect 100% success rate to stop repeating
        n_limit_checks = 100
        for i_check in range(n_limit_checks):

            print 'Checking status (call {0})'.format(i_check)

            subprocess.call( statuscmd, stdout=open(os.devnull, 'wb'))

            output = subprocess.Popen( retrievecmd, stdout=subprocess.PIPE ).communicate()[0]

            match = re.search( r'FAILED:\s*\d+\s*(\d+)', output )
            if int(match.group(1)) > 0:
                print 'There was a failed job. Check work.*/output/job*/stderr.gz to debug'
            
            match = re.search( r'SUCCESS:\s*(\d+)\s*(\d+)', output )
            n_success = match.group(1)
            p_success = int(match.group(2))

            match = re.search( r'RUNNING:\s*(\d+)\s*(\d+)', output )
            n_running = match.group(1)

            match = re.search( r'QUEUED:\s*(\d+)\s*(\d+)', output )
            n_queued = match.group(1)

            print '    Running: {0:4s} , Queued: {1:4s} , Finished: {2:4s}'.format(
                n_running, n_queued, n_success )

            if p_success == 100:
                print 'Jobs completed'
                break
            else:
                sleep(15)

        print 'Creating statistics reports'
        Create_sig_Report( '{0}/reporttest'.format(gc_loop_dir), i_iter )




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
