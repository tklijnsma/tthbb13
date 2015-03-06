#!/usr/bin/env python
"""
Thomas:
Step 3 towards a transfer function.

 - Specify input particles (which quarks/jets)
 - Specify which extra variables should be passed to the new N_Tuple
 - Loop:
   - Retrieve all data (4-vector & extra variables) from 1 event, cut-off pt
   - Link jets to quarks, cut-off delta R
   - Write all links to new .root file
   - Repeat

ADDED:
 - Manually calculates mcE and outputs it

"""

########################################
# Definitions
########################################

input_root_file_name = '/scratch/tklijnsm/VHBBHeppyV10_TTbarH_M-125_13TeV_PU20bx25.root'
input_tree_name = 'tree'

Pt_CutOff = 20  # Only work with (pseudo-)particles with pt > 30 GeV
R_CutOff = 0.1  # Only write links with R < 0.3

########################################
# Imports
########################################

import ROOT
import TTH.TTHNtupleAnalyzer.AccessHelpers as AH

########################################
# Functions
########################################

def Get_min_2D( Mat ):

    # Mat should be a list of lists: [ [...], [...], ...]

    if len(Mat)==0: print "Mat is empty!"

    Min = Mat[0][0]
    i_index = 0
    j_index = 0

    for i in range( len(Mat) ):
        for j in range( len(Mat[0]) ):
            if Mat[i][j] < Min:
                Min = Mat[i][j]
                i_index = i
                j_index = j

    return [ Min , i_index , j_index ]


def LinkJettoQuark( tl_jets, tl_quarks, R_CutOff ):

    # Count number of entries in both inputs
    if len(tl_jets) > len(tl_quarks):
        number_of_links = len(tl_quarks)
    else:
        number_of_links = len(tl_jets)

    # General method:
    #  - Compute matrix of delR values (rows = jets, columns = quarks)
    #  - Loop until out of quarks or jets:
    #    - Find minimum
    #    - Write the corresponding jet, quark and delR value to lists
    #    - Delete the row and column from the matrix
    #    - Repeat

    links = []
    tl_quarks_out = []
    tl_jets_out = []
    delR = []

    # Computing matrix of delR values
    Rmat = [[ (tl_jets[i].DeltaR( tl_quarks[j] )) for j in range(len(tl_quarks))] for i in range(len(tl_jets)) ]

    # Find the links
    for k in range(number_of_links):

        # Find minimum and indices of minimum in delR matrix
        [ Rmin, Jet_min_index, Quark_min_index ] = Get_min_2D( Rmat )
        
        # Don't use del R values > 0.3
        if Rmin > R_CutOff:
            break

        if False:
            # Print which quarks got linked to which jets
            print 'Minimum delR: ', Rmat[ Jet_min_index ][ Quark_min_index ], Rmin
            print 'linked Jet:'
            print [tl_jets[Jet_min_index].Pt(),
                tl_jets[Jet_min_index].Eta(),
                tl_jets[Jet_min_index].Phi(),
                tl_jets[Jet_min_index].M()]
            print 'to Quark:'
            print [tl_quarks[Quark_min_index].Pt(),
                tl_quarks[Quark_min_index].Eta(),
                tl_quarks[Quark_min_index].Phi(),
                tl_quarks[Quark_min_index].M()]
            print ''

        # Write corresponding jet, quark and delR-value to lists.
        # Also deletes these jets and quarks from the input lists, which is necessary
        # to find the right minimum indices.
        tl_jets_out.append( tl_jets.pop(Jet_min_index) )
        tl_quarks_out.append( tl_quarks.pop(Quark_min_index) )
        delR.append( Rmat[ Jet_min_index ][ Quark_min_index ] )
        
        # Delete jet and quark from delR matrix
        Rmat.pop(Jet_min_index)
        for row in Rmat:
            row.pop(Quark_min_index)

    return [tl_jets_out, tl_quarks_out, delR]


def Get_TLorentz( input_tree, particletype, is_jet, extra_vars, Pt_CutOff):

    val_ex_output = []
    tl_output = []
    
    # Get the variables for the TLorentz vector
    Pt =  AH.getter(input_tree, particletype + '_pt')
    Eta =  AH.getter(input_tree, particletype + '_eta')
    Phi =  AH.getter(input_tree, particletype + '_phi')
    Mass =  AH.getter(input_tree, particletype + '_mass')


    # Get the variables to also construct the TLorentz vector for the mc values
    # This is only needed to compute mcE, which is not in the input root file.
    if is_jet:
        mcPt =  AH.getter(input_tree, particletype + '_mcPt')
        mcEta =  AH.getter(input_tree, particletype + '_mcEta')
        mcPhi =  AH.getter(input_tree, particletype + '_mcPhi')
        mcM =  AH.getter(input_tree, particletype + '_mcM')


    # Get the values for the extra variables (beyond what is necessary for TLorentz)
    for var in extra_vars:
        val_ex_output.append( AH.getter(input_tree, particletype + '_' + var) )

    # Remove duplicates: 2 particles are often repeated in the WZQuark entries
    removedupl = 0
    if len(Pt)>2:
        if Pt[0]==Pt[ len(Pt)-2 ] and Pt[1]==Pt[ len(Pt)-1 ]:
            removedupl = 2

    # Cut off pt<30, write pt>30 to TLorentzVector
    for i in range( len(Pt) - removedupl ):
        if Pt[i] > Pt_CutOff:

            y = ROOT.TLorentzVector()
            y.SetPtEtaPhiM( Pt[i] , Eta[i] , Phi[i] , Mass[i] )

            # Fill in the extra variables
            for var, extra_val_list in zip(extra_vars, val_ex_output):
                setattr(y, var, extra_val_list[i])

            # Fill in mcE
            if is_jet:
                x = ROOT.TLorentzVector()
                x.SetPtEtaPhiM( mcPt[i] , mcEta[i] , mcPhi[i] , mcM[i] )
                setattr(y, 'mcE', x.E() )

            tl_output.append( y )
                
    return tl_output




########################################
# Main
########################################

def main():

    ROOT.gROOT.SetBatch(True)
    
    # Linking switch: If Just_Jets=True, .root file will simply contain (unlinked)
    # jets.
    Just_Jets = False

    output_root_file_name = 'Jets_with_Quarks.root'

    quarktypes = ['GenBQuarkFromTop', 'GenBQuarkFromH', 'GenWZQuark' ]
    jettypes = ['aJets']

    standard_vars = [ 'pt', 'eta', 'phi', 'mass', 'E' ]

    extra_quark_vars = [ 'pdgId', 'charge', 'status' ]
    extra_jet_vars = [ 'mcPt', 'mcEta', 'mcPhi', 'mcM',
        'mcFlavour', 'mcMatchId', 'btagCSV' ]

    # Only for output purposes
    separate_vars = ['delR', 'Jet_mcE']


    ########################################
    # Setup I/O
    ########################################

    # Input tree
    input_root_file = ROOT.TFile(input_root_file_name)
    input_tree = input_root_file.Get(input_tree_name)

    # Output tree
    output_root_file = ROOT.TFile(output_root_file_name,'RECREATE')
    output_tree = ROOT.TTree('tree','My test tree')

    # Define branches in output tree
    branches = []

    branches.extend( [ 'Jet_' + var for var in standard_vars ] )
    branches.extend( [ 'Jet_' + var for var in extra_jet_vars ] )

    branches.extend( [ 'Quark_' + var for var in standard_vars ] )
    branches.extend( [ 'Quark_' + var for var in extra_quark_vars ] )

    branches.extend( [ var for var in separate_vars ] )

    # Create dicitionaries to hold the information that will be
    # written as new branches
    variables      = {}
    variable_types = {}

    # Setup the output branches for the true object
    AH.addScalarBranches(variables,
                         variable_types,
                         output_tree,
                         branches,
                         datatype = 'float')


    ########################################
    # Event loop
    ########################################
    
    n_entries = input_tree.GetEntries()
    print "Processing {0} events".format(n_entries)

    for i_event in range(n_entries):
    #for i_event in range(5000):

        #print '******* Event', i_event
        if not i_event % 1000:
            print "{0:.1f}%".format( 100.*i_event /n_entries)

        input_tree.GetEntry( i_event )

        ########################################
        # Get quark and jet data
        #  - list of TLorentzVectors
        #  - each TLorentzVector has extra variables as attributes
        ########################################

        tl_quarks = []

        for quarktype in quarktypes:
            tl_quarks.extend( Get_TLorentz( input_tree,
                quarktype,
                False,
                extra_quark_vars,
                Pt_CutOff ) )

        tl_jets = []

        for jettype in jettypes:
            if Just_Jets == True:
                # Set Pt_CutOff to 20, fixed
                tl_jets.extend( Get_TLorentz( input_tree,
                    jettype,
                    True,
                    extra_jet_vars,
                    20 ) )
            else:
                tl_jets.extend( Get_TLorentz( input_tree,
                    jettype,
                    True,
                    extra_jet_vars,
                    Pt_CutOff ) )


        ########################################
        # delR combinatorics
        ########################################

        if Just_Jets == True:
            # Don't perform linking if 'Just_Jets' parameter is true
            # In that case, just write all found jet data
            for jet in tl_jets:

                variables['Jet_pt'][0] = jet.Pt()
                variables['Jet_eta'][0] = jet.Eta()
                variables['Jet_phi'][0] = jet.Phi()
                variables['Jet_mass'][0] = jet.M()

                variables['Jet_E'][0] = jet.E()
                variables['Jet_mcE'][0] = jet.mcE

                # Retrieve the extra variables from set attributes
                for var in extra_jet_vars:
                    variables['Jet_'+var][0] = getattr( jet, var )

                output_tree.Fill()
            continue
        
        # Otherwise, proceed with linking quarks and jets    
        [linked_jets, linked_quarks, delRs] = LinkJettoQuark( tl_jets,
            tl_quarks,
            R_CutOff )

        # linked_jets and linked_quarks are (ordered) lists of TLorentzVectors. The
        # extra variables are contained in attributes.
        
        ########################################
        # Write to file
        ########################################

        for jet, quark, delR in zip( linked_jets, linked_quarks, delRs):

            variables['Jet_pt'][0] = jet.Pt()
            variables['Jet_eta'][0] = jet.Eta()
            variables['Jet_phi'][0] = jet.Phi()
            variables['Jet_mass'][0] = jet.M()

            variables['Jet_E'][0] = jet.E()
            variables['Jet_mcE'][0] = jet.mcE

            variables['Quark_pt'][0] = quark.Pt()
            variables['Quark_eta'][0] = quark.Eta()
            variables['Quark_phi'][0] = quark.Phi()
            variables['Quark_mass'][0] = quark.M()

            variables['Quark_E'][0] = quark.E()

            variables['delR'][0] = delR

            # Retrieve the extra variables from set attributes
            for var in extra_quark_vars:
                variables['Quark_'+var][0] = getattr( quark, var )
                
            for var in extra_jet_vars:
                variables['Jet_'+var][0] = getattr( jet, var )

            output_tree.Fill()


    output_root_file.Write()
    output_root_file.Close()


########################################
# End of main
########################################   


if __name__ == "__main__":
    main()
