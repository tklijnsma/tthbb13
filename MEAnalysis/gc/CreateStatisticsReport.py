#!/usr/bin/env python
"""
Thomas:

Retrieves the statistics from individual in job.stdout.gz files.

"""

########################################
# Imports
########################################

import re
import os
import gzip

from operator import add

########################################
# Functions
########################################

def Create_Statistics_Report( workdir, outdir ):


    ########################################
    # Fill the dict
    ########################################

    dirs = os.listdir( 'work.{0}/output/'.format(workdir) )

    for jobdir in dirs:

        stdout_filename = 'work.{0}/output/{1}/job.stdout.gz'.format(workdir,jobdir)

        if not os.path.isfile( stdout_filename ):
            print '{0} does not exist (yet)'.format( stdout_filename )
            continue

        f = gzip.open(stdout_filename, 'rb')
        full_out = f.read()
        f.close()


        begin_match = re.search( r'Statistics\n==========', full_out )
        end_match = re.search( r'==========\nEnd of Statistics', full_out )
        
        if not begin_match or not end_match:
            print 'Could not find Statistics'
            return

        begin_text = begin_match.end()
        end_text = end_match.start()

        report = full_out[ begin_text : end_text ]

        # Get the key_list, runs just once
        if 'key_list' not in locals():

            # Create the key_list
            key_matches = re.finditer( r'([\w\+]+)\s+=\s\d+', report )

            key_list = []
            for key_match in key_matches:

                key = key_match.group(1)

                if report[key_match.end()] == ')':
                    key_list.append('#' + key)
                else:
                    key_list.append(key)

            # Initialize the stat_dict
            stat_dict = {}
            for key in key_list:
                if key[0] == '#':
                    stat_dict[key[1:]] = 0
                else:
                    stat_dict[key] = 0


        key_matches = re.findall( r'([\w\+]+)\s+=\s(\d+)', report )

        for (key, value) in key_matches:
            stat_dict[key] += int(value)
        

    ########################################
    # Write the dict
    ########################################

    out_f = open( '{0}/StatisticsReport-{1}.txt'.format(outdir, workdir) , 'w' )

    for key in key_list:

        if key[0] == '#':
            print '  ({0:25s} = {1})'.format( key[1:], stat_dict[key[1:]] )
            out_f.write( '  ({0:25s} = {1})\n'.format( key[1:], stat_dict[key[1:]]))
        else:
            print '{0:30s} = {1}'.format( key, stat_dict[key] )
            out_f.write( '{0:30s} = {1}\n'.format( key, stat_dict[key] ) )

    out_f.close()


def Create_Match_Tables_Report( workdir, outdir ):

    ########################################
    # Fill the dict
    ########################################

    dirs = os.listdir( 'work.{0}/output/'.format(workdir) )

    total_match_table_b = [ 0 for i in range(8) ]
    total_match_table_l = [ 0 for i in range(8) ]


    for jobdir in dirs:

        stdout_filename = 'work.{0}/output/{1}/job.stdout.gz'.format(workdir,jobdir)

        if not os.path.isfile( stdout_filename ):
            print '{0} does not exist (yet)'.format( stdout_filename )
            continue

        f = gzip.open(stdout_filename, 'rb')
        full_out = f.read()
        f.close()

        begin_match = re.search( r'Match Tables\n==========', full_out )
        end_match = re.search( r'==========\nEnd of Match Tables', full_out )
        
        if not begin_match or not end_match:
            print 'Could not find Match Tables'
            return

        begin_text = begin_match.end()
        end_text = end_match.start()

        report = full_out[ begin_text : end_text ]

        count_matches = re.findall( r'\D\d+\D', report )

        clist = [ int(count_match) for count_match in count_matches ]
        
        match_table_b = clist[:8]
        match_table_l = clist[8:]

        total_match_table_b = map(add, match_table_b, total_match_table_b)
        total_match_table_l = map(add, match_table_l, total_match_table_l)

    Print_Match_Table( 'b', total_match_table_b )
    Print_Match_Table( 'l', total_match_table_l )

    Write_Match_Table( 'b', total_match_table_b, workdir, outdir )
    Write_Match_Table( 'l', total_match_table_l, workdir, outdir )

        



# Prints the matching tables
def Print_Match_Table( particle, match_table ):

    print '\n{0:23s}| {1:23s}| {2:23s}|'.format(
        'Match table: {0}'.format(particle),
        'suc. {0}-quark to {0}-jet'.format(particle),
        'fail {0}-quark to {0}-jet'.format(particle) )

    print '-----------------------|------------------------|------------------------|'

    print '{0:23s}| {1:23s}| {2:23s}|'.format(
        'suc. {0}-quark to {0}-sj'.format(particle),
        str(match_table[0]) ,
        str(match_table[1]) )

    print '{0:23s}| {1:23s}| {2:23s}|'.format(
        '   suc. jet-sj',
        str(match_table[2]) ,
        str(match_table[3]) )

    print '-----------------------|------------------------|------------------------|'

    print '{0:23s}| {1:23s}| {2:23s}|'.format(
        'fail {0}-quark to {0}-sj'.format(particle),
        str(match_table[4]) ,
        str(match_table[5]) )

    print '{0:23s}| {1:23s}| {2:23s}|'.format(
        '   suc. jet-sj',
        str(match_table[6]) ,
        str(match_table[7]) )

    print '--------------------------------------------------------------------------'


# Writes the matching tables
def Write_Match_Table( particle, match_table, workdir, outdir ):

    if os.path.isfile('{0}/StatisticsReport-{1}.txt'.format(outdir,workdir) ):
        f = open( '{0}/StatisticsReport-{1}.txt'.format(outdir,workdir) , 'a' )
    else:
        f = open( '{0}/StatisticsReport-{1}.txt'.format(outdir,workdir) , 'w' )


    f.write( '\n{0:23s}| {1:23s}| {2:23s}|\n'.format(
        'Match table: {0}'.format(particle),
        'suc. {0}-quark to {0}-jet'.format(particle),
        'fail {0}-quark to {0}-jet'.format(particle) ) )

    f.write( '-----------------------|------------------------|------------------------|\n' )

    f.write( '{0:23s}| {1:23s}| {2:23s}|\n'.format(
        'suc. {0}-quark to {0}-sj'.format(particle),
        str(match_table[0]) ,
        str(match_table[1]) ) )

    f.write( '{0:23s}| {1:23s}| {2:23s}|\n'.format(
        '   suc. jet-sj',
        str(match_table[2]) ,
        str(match_table[3]) ) )

    f.write( '-----------------------|------------------------|------------------------|\n' )

    f.write( '{0:23s}| {1:23s}| {2:23s}|\n'.format(
        'fail {0}-quark to {0}-sj'.format(particle),
        str(match_table[4]) ,
        str(match_table[5]) ) )

    f.write( '{0:23s}| {1:23s}| {2:23s}|\n'.format(
        '   suc. jet-sj',
        str(match_table[6]) ,
        str(match_table[7]) ) )

    f.write( '--------------------------------------------------------------------------\n' )


    f.close()



def Create_Specific_Print_Report( workdir, outdir ):

    ########################################
    # Fill the dict
    ########################################

    dirs = os.listdir( 'work.{0}/output/'.format(workdir) )

    out_f = open( '{0}/SpecificPrintReport-{1}.txt'.format(outdir,workdir) , 'w' )

    for jobdir in dirs:

        stdout_filename = 'work.{0}/output/{1}/job.stdout.gz'.format(workdir,jobdir)

        if not os.path.isfile( stdout_filename ):
            print '{0} does not exist (yet)'.format( stdout_filename )
            continue

        out_f.write( 'Specific Prints in {0}\n====================\n'.format(jobdir) )

        f = gzip.open(stdout_filename, 'rb')
        full_out = f.read()
        f.close()

        begin_match = re.search( r'Specific Print\n==========', full_out )
        end_match = re.search( r'==========\nEnd of Specific Print', full_out )
        
        if not begin_match or not end_match:
            out_f.write( 'Could not find Specific Print\n\n' )
            continue

        begin_text = begin_match.end()
        end_text = end_match.start()

        specific_print = full_out[ begin_text : end_text ]

        out_f.write( specific_print )
        out_f.write( '\n' )

        



########################################
# Main
########################################

def main():

    outdir = 'Reports'
    
    if not os.path.isdir( 'Reports'):
        os.makedirs( 'Reports' )

    Create_Statistics_Report( 'bkg', outdir )
    Create_Match_Tables_Report( 'bkg', outdir )
    Create_Specific_Print_Report( 'bkg', outdir )

    Create_Statistics_Report( 'sig', outdir )
    Create_Match_Tables_Report( 'sig', outdir )
    Create_Specific_Print_Report( 'sig', outdir )




########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
