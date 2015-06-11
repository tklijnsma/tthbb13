#!/usr/bin/env python
"""
Thomas:

Examples on how to use TFMatrix.dat

2 plots are created: 1 transfer function, and 1 cumulative distribution function

"""

########################################
# Imports
########################################

import os
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

    # Open directory for output plots
    outputdir = 'TFtest_plots'
    if not os.path.isdir( outputdir ):
        os.makedirs( outputdir )

    # Open TFMatrix.dat to TFmat
    pickle_f = open( '../V11_full_subjets_REVERSED_BRANCHES/TFMatrix.dat', 'rb' )
    #pickle_f = open( '../V11_full_subjets/TFMatrix.dat', 'rb' )
    TFmat = pickle.load( pickle_f )
    pickle_f.close()

    # Choose a TF; As an example the TF for a b in the first eta-bin is chosen.
    #   (b is fitted with double Gaussians, l with single Gaussians)
    myTF = TFmat['b'][0]

    # Choose a pt value to test with
    test_pt = 78.0


    ########################################
    # Plot an example of a response function
    ########################################

    # Determine Make_Formula output:
    # - True:  [0] = reconstr., x = mc/gen/quark
    # - False: [0] = mc/gen/quark, x = reconstr.
    # No argument is treated as True
    set_reconstructed_eval_gen = False

    f1 = myTF.Make_Formula( set_reconstructed_eval_gen )

    # Set gen pt to [test_pt] GeV
    f1.SetParameter( 0, test_pt )

    # Print some specifics of the created TF
    print '\nSpecifics of the TF:'
    print f1.GetTitle()
    for i in range( f1.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f1.GetParameter(i) )

    # Retrieve the RMS's of the two Gaussians ( [2] and [2]+[4] )
    RMS_G1 = myTF.AcrossBinFuncs[2].Evaluate(test_pt)

    RMS_G2_addition = myTF.AcrossBinFuncs[4].Evaluate(test_pt)
    RMS_G2 = RMS_G1 + RMS_G2_addition

    expected_integral = pow( 2.0*pi , 0.5 ) * ( 0.7*RMS_G1 + 0.3*RMS_G2 )
    print '\nTF normalized by: sqrt(2*pi)*(0.7*RMS1+0.3*RMS2)'
    print '     = {0}'.format( expected_integral )

    actual_integral = f1.Integral( -40.0, 500.0 )
    print '\nNumerical integral of TF:'
    print '     = {0}'.format( actual_integral )


    print '\n====================='
    print 'RMS_G1 function:'
    f_RMS_G1 = myTF.AcrossBinFuncs[2].Get_TF1_with_set_parameters()
    print f_RMS_G1.GetTitle()
    for i in range( f_RMS_G1.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f_RMS_G1.GetParameter(i) )
    print '    Evaluated at {1}: RMS = {0}'.format( RMS_G1, test_pt )

    print '\n====================='
    print 'RMS_G2 function:'
    f_RMS_G2 = myTF.AcrossBinFuncs[4].Get_TF1_with_set_parameters()
    print f_RMS_G2.GetTitle()
    for i in range( f_RMS_G2.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f_RMS_G2.GetParameter(i) )
    print '    Evaluated at {1}: RMS addition = {0}'.format(RMS_G2_addition,test_pt)
    print '    RMS = RMS1 + RMS addition       = {0}'.format( RMS_G2 )


    # Drawing
    f1.SetRange( 0.0 , 2.5*test_pt )
    c1 = ROOT.TCanvas("c1","c1",500,400)
    c1.SetGrid()
    f1.Draw()
    c1.Update()
    c1.Print("{0}/exampleplot_TF".format(outputdir),"pdf")
    img = ROOT.TImage.Create()
    img.FromPad(c1)
    img.WriteImage("{0}/exampleplot_TF.png".format(outputdir))


    # Example plot withe pt_reco = 58.0, pt_gen = x
    # ======================================

    f_rev = myTF.Make_Formula( True )

    # Set gen pt to [test_pt] GeV
    f_rev.SetParameter( 0, test_pt )
    
    print '\nSpecifics of TF ( [0]:reco_pt, x:gen_pt ) :'
    print f_rev.GetTitle()
    for i in range( f_rev.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f_rev.GetParameter(i) )

    rev_integral = f_rev.Integral( -40.0, 500.0 )
    print '\nNumerical integral of TF:'
    print '     = {0}'.format( rev_integral )

    # Drawing
    f_rev.SetRange( 0.0 , 2.5*test_pt )
    f_rev.Draw()
    c1.Update()
    c1.Print("{0}/exampleplot_TF_set_reco".format(outputdir),"pdf")
    img = ROOT.TImage.Create()
    img.FromPad(c1)
    img.WriteImage("{0}/exampleplot_TF_set_reco.png".format(outputdir))

    ########################################
    # Plot an example of a CDF
    ########################################

    # Create the cumulative distribution function
    #  - Only works for single or double Gaussian functions
    #  - [0] = mc/gen/quark pt, x = pt_cutoff
    f2 = myTF.Make_CDF()

    # Set the mc/gen/quark pt
    f2.SetParameter( 0, test_pt )

    print '\nSpecifics of the CDF:'
    print f2.GetTitle()
    for i in range( f2.GetNumberFreeParameters() ):
        print '    [{0}] = {1}'.format( i, f2.GetParameter(i) )

    # Drawing
    f2.SetRange( 0.0 , 2.0*test_pt )
    f2.Draw()
    c1.Update()
    c1.Print("{0}/exampleplot_CDF".format(outputdir),"pdf")
    img = ROOT.TImage.Create()
    img.FromPad(c1)
    img.WriteImage("{0}/exampleplot_CDF.png".format(outputdir))

    hf = open( '{0}/0_overview.html'.format(outputdir), 'w' )
    hf.write( '<html><body>\n<h1>Test of Transfer Function at pt = {0} GeV:\n</h1>\n<br>\n<hr />'.format( test_pt ) )

    hf.write( \
        """
        <table style="width:100%">
          <tr>
            <td>
                gen_pt={5} set, looking at reco_pt response
                <br>
                <a href="{0}"><img width="450" src="{0}.png"></a>
                <br>
                integral over whole function: {3}
            </td>
            <td>
                reco_pt={5} set, looking at gen_pt response
                <br>
                <a href="{1}"><img width="450" src="{1}.png"></a>
                <br>
                integral over whole function: {4}
            </td>
            <td>
                pt_cutoff={5} set, looking at lost fraction of events
                <br>
                <a href="{2}"><img width="450" src="{2}.png"></a>
            </td>
          </tr>
        </table>
        """.format( 'exampleplot_TF',
                    'exampleplot_TF_set_reco',
                    'exampleplot_CDF',
                    actual_integral,
                    rev_integral,
                    test_pt,
                    ))

    #for fn in [ 'exampleplot_TF', 'exampleplot_TF_set_reco', 'exampleplot_CDF' ]:
    #    hf.write('\n<a href="{0}"><img width="450" src="{0}.png"></a>\n'.format(fn) )
    hf.close()



########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
