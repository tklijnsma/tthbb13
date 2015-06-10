#!/usr/bin/env python
"""
Thomas:

Examples on how to use TFMatrix.dat

2 plots are created: 1 transfer function, and 1 cumulative distribution function

"""

########################################
# Imports
########################################

import pickle
import ROOT
from math import pi


########################################
# Main
########################################

def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
    ROOT.gStyle.SetOptFit(1011)


    # Open TFMatrix.dat to TFmat
    pickle_f = open( '../V11_full_jets/TFMatrix.dat', 'rb' )
    TFmat = pickle.load( pickle_f )
    pickle_f.close()

    # Choose a TF; As an example the TF for a b in the first eta-bin is chosen.
    #   (b is fitted with double Gaussians, l with single Gaussians)
    myTF = TFmat['b'][0]


    ########################################
    # Plot an example of a response function
    ########################################

    # Determine Make_Formula output:
    # - True:  [0] = reconstr., x = mc/gen/quark
    # - False: [0] = mc/gen/quark, x = reconstr.
    # No argument is treated as True
    set_reconstructed_eval_gen = False

    f1 = myTF.Make_Formula( set_reconstructed_eval_gen )

    # Set gen pt to 58 GeV
    f1.SetParameter( 0, 58.0 )

    print 'Specifics of the TF:'
    print f1.GetTitle()
    for i in range( f1.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f1.GetParameter(i) )


    print '\n====================='
    print 'RMS_G1 function:'
    f_RMS_G1 = myTF.AcrossBinFuncs[2].Get_TF1_with_set_parameters()
    print f_RMS_G1.GetTitle()
    for i in range( f_RMS_G1.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f_RMS_G1.GetParameter(i) )
    RMS_G1 = myTF.AcrossBinFuncs[2].Evaluate(58.0)
    print '    Evaluated at 58.0: RMS = {0}'.format( RMS_G1 )


    print '\n====================='
    print 'RMS_G2 function:'
    f_RMS_G2 = myTF.AcrossBinFuncs[4].Get_TF1_with_set_parameters()
    print f_RMS_G2.GetTitle()
    for i in range( f_RMS_G2.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f_RMS_G2.GetParameter(i) )
    RMS_G2_addition = myTF.AcrossBinFuncs[4].Evaluate(58.0)
    print '    Evaluated at 58.0: RMS addition = {0}'.format( RMS_G2_addition )
    RMS_G2 = RMS_G1 + RMS_G2_addition
    print '    RMS = RMS1 + addition           = {0}'.format( RMS_G2 )


    expected_integral = pow( 2.0*pi , 0.5 ) * ( 0.7*RMS_G1 + 0.3*RMS_G2 )
    print '\nExpected unnorm. integral: sqrt(2*pi)*(0.7*RMS1+0.3*RMS2)'
    print '     = {0}'.format( expected_integral )

    actual_integral = f1.Integral( -40.0, 150.0 )
    print '\nActual integral of TF1:'
    print '     = {0}'.format( actual_integral )




    # Drawing
    f1.SetRange( 0.0 , 100.0)
    c1 = ROOT.TCanvas("c1","c1",500,400)
    c1.SetGrid()
    f1.Draw()
    c1.Update()
    c1.Print("exampleplot_TF","pdf")


    """
    ########################################
    # Plot an example of a CDF
    ########################################

    # Create the cumulative distribution function
    #  - Only works for single or double Gaussian functions
    #  - [0] = mc/gen/quark pt, x = pt_cutoff
    f2 = myTF.Make_CDF()

    # Set the mc/gen/quark pt
    f2.SetParameter( 0, 58.0 )

    print 'Specifics of the CDF:'
    print f2.GetTitle()
    for i in range( f2.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f2.GetParameter(i) )

    # Drawing
    f2.SetRange( 0.0 , 100.0)
    f2.Draw()
    c1.Update()
    c1.Print("exampleplot_CDF","pdf")
    """


########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
