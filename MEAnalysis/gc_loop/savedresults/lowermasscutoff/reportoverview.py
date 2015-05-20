#!/usr/bin/env python
"""
Thomas:

"""

########################################
# Imports
########################################

import re
import os
import ROOT


########################################
# Functions
########################################

def Return_Number_From_Str( string ):

    match = re.search( r'\d+', string )
    return int(match.group())


def Extract_Key_From_Statfile( key, statdir ):

    n_list = []

    statfiles = sorted( os.listdir(statdir), key=Return_Number_From_Str )

    for statfile in statfiles:
        f = open( statdir+statfile, 'r' )
        report = f.read()
        f.close()

        match = re.search( r'{0}\s*=\s*(\d+)'.format(key) , report )
        n_list.append( match.group(1) )

    return n_list


def Extract_Number_From_Treefile( i, treedir ):

    n_list = []

    treefiles = sorted( os.listdir(treedir), key=Return_Number_From_Str )

    for treefile in treefiles:
        f = open( treedir+treefile, 'r' )
        report = f.read()
        f.close()

        matches = re.findall( r'\D(\d+)\D' , report )
        n_list.append( matches[i] )

    return n_list

    
def ROOT_Line_Plot( var_list, n_list, label_tup, hf, c1 ):

    #c1 = ROOT.TCanvas("c1","c1",500,400)
    #c1.SetGrid()

    ( title, xlabel, ylabel ) = label_tup

    gr = ROOT.TGraph( len(n_list) )

    for i in range(len(n_list)):
        gr.SetPoint( i, float(var_list[i]), float(n_list[i]) )

    gr.Draw('A*')

    gr.SetTitle("{0};{1};{2}".format(title,xlabel,ylabel));

    c1.Update()

    filename = 'plots/' + title

    # Write pdf
    c1.Print( filename, "pdf" )

    # Write png
    print 'Writing {0}.png'.format(filename)
    img = ROOT.TImage.Create()
    img.FromPad(c1)
    img.WriteImage( '{0}.png'.format(filename) )

    # Write line to html
    hf.write('<a href="{0}"><img width="700" src="{0}.png"></a>\n'.format(filename) )




########################################
# Main
########################################

def main():

    ########################################
    # Extraction
    ########################################

    var_list = [ '80.0', '100.0', '120.0', '130.0', '135.0', '140.0', '145.0', '150.0', '155.0', '160.0', '180.0' ]

    reportdir = 'reporttest'

    statdir = '{0}/Statistics/'.format(reportdir)
    treedir = '{0}/Match_Trees/'.format(reportdir)
    printdir = '{0}/Specific_Prints/'.format(reportdir)


    Numbers_to_get = [

        [ 'stat' , 'n_survivedcut' ] ,
        [ 'stat' , 'n_0_cand'      ] ,
        [ 'tree' , 'n_suc_jet_fail_sj', 6 ] ,
        [ 'ratio', ( 'tree', 'n_suc_jet_fail_sj', 6 ), ( 'stat', 'n_survivedcut' ) ],


        ]


    ########################################
    # Retrieving and Drawing
    ########################################

    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

    # Don't display standard statistics
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)

    # Open an html-file
    hf = open( 'reportoverview.html', 'w' )
    hf.write( '<html><body>\n<h1>Run Summary:\n</h1>\n<br>\n<hr />' )
    hf.write( '<h2>Title</h2>' )

    # Setup the canvas
    c1 = ROOT.TCanvas("c1","c1",500,400)
    c1.SetGrid()

    for arglist in Numbers_to_get:
        
        if arglist[0] == 'stat':
            n_list = Extract_Key_From_Statfile( arglist[1], statdir )
        elif arglist[0] == 'tree':
            n_list = Extract_Number_From_Treefile( arglist[2] , treedir )
        n_name = arglist[1]

        if arglist[0] == 'ratio':

            tup1 = arglist[1]
            if tup1[0] == 'stat':
                n1_list = Extract_Key_From_Statfile( tup1[1], statdir )
            if tup1[0] == 'tree':
                n1_list = Extract_Number_From_Treefile( tup1[2] , treedir )
            n1_name = tup1[1]

            tup2 = arglist[2]
            if tup2[0] == 'stat':
                n2_list = Extract_Key_From_Statfile( tup2[1], statdir )
            if tup2[0] == 'tree':
                n2_list = Extract_Number_From_Treefile( tup2[2] , treedir )
            n2_name = tup2[1]

            n_name = '{0}_OVER_{1}'.format( n1_name, n2_name )
            n_list = [ str(float(n1)/float(n2)) for (n1,n2) in zip( n1_list, n2_list ) ]


        label_tup = ( n_name, 'lower mass cutoff (GeV)', n_name )
        ROOT_Line_Plot( var_list, n_list, label_tup, hf, c1 )




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
