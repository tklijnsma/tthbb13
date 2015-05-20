import ROOT
import itertools
from PhysicsTools.HeppyCore.framework.analyzer import Analyzer
import copy
import math
from TTH.MEAnalysis.VHbbTree import *

#Load the MEM integrator libraries
# ROOT.gSystem.Load("libFWCoreFWLite")
# ROOT.gROOT.ProcessLine('AutoLibraryLoader::enable();')
# ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gSystem.Load("libCintex")
ROOT.gROOT.ProcessLine('ROOT::Cintex::Cintex::Enable();')
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")

from ROOT import MEM
o = MEM.MEMOutput

#Pre-define shorthands for permutation and integration variable vectors
CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")
CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")


def lvec(self):
    """
    Converts an object with pt, eta, phi, mass to a TLorentzVector
    """
    lv = ROOT.TLorentzVector()
    lv.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.mass)
    return lv


# Very basic tree node class
class TreeNode:

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.children = []

    def Print(self):

        print 'Class TreeNode: {0}'.format( self.split_name )
        print '    n = {0}'.format( self.n )

        if len(self.children) > 0:
            child_names = [ child.name for child in self.children ]
            print '    children = {0}'.format( child_names )
        else:
            print '    children = []'

    def Get_child( self, *input_args ):

        if len( input_args ) == 0:
            print "Format: <TreeNode>.Get_child( <1/0>, <1/0>, ... )"
            return 0

        t = self

        for i in input_args:

            if len( t.children ) == 0:
                print 'Too many input arguments - '\
                      'Tree does not have that many children'
                return 0

            t = t.children[i]

        return t
#--------------------------------------#


class FilterAnalyzer(Analyzer):
    """
    A generic analyzer that may filter events.
    Counts events the number of processed and passing events.
    """
    def beginLoop(self, setup):
        super(FilterAnalyzer, self).beginLoop(setup)
        self.counters.addCounter("processing")
        self.counters["processing"].register("processed")
        self.counters["processing"].register("passes")


class SubjetAnalyzer(FilterAnalyzer):
    """
    Subjet analyzer by Thomas
    """

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(SubjetAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf

        self.verbose = False
        self.require_quark_match = True

        self.R_cut = 0.3
        self.top_mass = 172.04

        if hasattr( self.conf, 'httCandidatecut' ):
            self.Cut_criteria = self.conf.httCandidatecut
            print 'Using httCandidate cut criteria from configuration file'
        else:
            self.Cut_criteria = [
                ( 'pt'  , '>', '200.0' ),
                ( 'mass', '>', '120.0' ),
                ( 'mass', '<', '220.0' ),
                ( 'fW'  , '<', '0.175' ) ]


        self.Statistics_keylist = [
            'n_processed',
            'n_cat1+btag',

            'n_enough_initial_jets',
            '#n_not_enough_wcand',
            '#n_not_enough_btagged',
            '#n_not_enough_httCand',

            'n_survivedcut',
            '#n_0_cand',
            '#n_1_cand',
            '#n_2+_cand',

            'n_right_nr_of_quarks',
            '#n_too_few_WZ',
            '#n_too_few_B',
            '#n_too_many_WZ',
            '#n_too_many_B',

            'n_passed_to_MEM',

            ]

        self.Statistics = {}
        for key in self.Statistics_keylist:
            if key[0] == '#':
                self.Statistics[key[1:]] = 0
            else:
                self.Statistics[key] = 0


        # Create match tree for b-particles

        self.match_tree_b = TreeNode( 'b' )

        # 1st level
        self.match_tree_b.children = [ TreeNode('Jet_Qrk_fail'),
                                       TreeNode('Jet_Qrk_success') ]

        # 2nd level
        b = self.match_tree_b.Get_child( 0 )
        b.children = [ TreeNode('Sj_Qrk_fail'), TreeNode('Sj_Qrk_success') ]

        b = self.match_tree_b.Get_child( 1 )
        b.children = [ TreeNode('Sj_Qrk_fail'), TreeNode('Sj_Qrk_success') ]

        # 3rd level
        b = self.match_tree_b.Get_child( 0, 0 )
        b.children = [ TreeNode('Jet_Sj_fail'), TreeNode('Jet_Sj_success') ]

        b = self.match_tree_b.Get_child( 0, 1 )
        b.children = [ TreeNode('Jet_Sj_fail'), TreeNode('Jet_Sj_success') ]

        b = self.match_tree_b.Get_child( 1, 0 )
        b.children = [ TreeNode('Jet_Sj_fail'), TreeNode('Jet_Sj_success') ]

        b = self.match_tree_b.Get_child( 1, 1 )
        b.children = [ TreeNode('Jet_Sj_fail'), TreeNode('Jet_Sj_success') ]


        # Create match tree for l-particles

        self.match_tree_l = copy.deepcopy( self.match_tree_b )
        self.match_tree_l.name = 'l'



    def beginLoop(self, setup):
        super(SubjetAnalyzer, self).beginLoop(setup)


    def endLoop(self, setup):

        print '\nStatistics'
        print '=========='

        for key in self.Statistics_keylist:

            if key[0] == '#':
                print '  ({0:25s} = {1})'.format(key[1:], self.Statistics[key[1:]] )
            else:
                print '{0:30s} = {1}'.format( key, self.Statistics[key] )

        print '=========='
        print 'End of Statistics\n'


        print '\nMatch Trees'
        print '=========='

        self.Print_Tree( self.match_tree_b )
        self.Print_Tree( self.match_tree_l )

        print '=========='
        print 'End of Match Trees\n'



    def process(self, event):

        self.Statistics['n_processed'] += 1

        #if self.verbose:
        #    print 'Printing from SubjetAnalyzer! iEv = {0}'.format(event.iEv)

        print 'Printing from SubjetAnalyzer! iEv = {0}'.format(event.iEv)


        ########################################
        # Check event suitability
        ########################################

        # Check if event is category 1
        if not ( event.cat == 'cat1' and event.cat_btag == 'H' ):
            return 0
        self.Statistics['n_cat1+btag'] += 1

        # Check if there are at least 2 wquark_candidate_jets
        if len( event.wquark_candidate_jets ) < 2:
            self.Statistics['n_not_enough_wcand'] += 1
            return 0

        # Check if there is at least 1 btagged_jet
        if len( event.btagged_jets ) == 0:
            self.Statistics['n_not_enough_btagged'] += 1
            return 0

        # Check if there is an httCandidate
        if len( event.httCandidate ) == 0:
            self.Statistics['n_not_enough_httCand'] += 1
            return 0

        self.Statistics['n_enough_initial_jets'] += 1

        # Apply the cuts
        tops = []
        for candidate in event.httCandidate:
            if self.Apply_Cut_criteria( candidate ):
                tops.append( copy.deepcopy(candidate) )

        other_top_present = False

        # Check if any candidates survived the cutoff criteria
        if len(tops) == 0:
            if self.verbose: print 'No candidate survived the cuts'
            self.Statistics['n_0_cand'] += 1
            return 0

        elif len(tops) == 1:
            # SortTops sets the delR_lepton attribute, so also use it for 1 top
            tops = self.SortTops( event, tops )
            top = tops[0]
            self.Statistics['n_1_cand'] += 1

        # If more than 1 candidate survived the cutoff criteria, choose the
        # one with a mass closest to top mass

        else:
            #tops = sorted( tops, key=lambda x: abs(x.mass - self.top_mass) )
            tops = self.SortTops( event, tops )
            other_top_present = True
            top = tops[0]
            other_top = tops[1]
            self.Statistics['n_2+_cand'] += 1

        self.Statistics['n_survivedcut'] += 1

        # Write delRmin to event
        setattr( event, 'httCandidate_delRmin', abs(top.Rmin-top.RminExpected) )
        print 'httCandidate_delRmin = {0}'.format( abs(top.Rmin-top.RminExpected) )

        # Necessary to remove duplicates in GenWZQuark branch
        self.CompareWZQuarks( event )


        ########################################
        # Get the lists of particles: quarks, jets and subjets
        ########################################

        # Note: if self.require_quark_match is set to False, this part is only
        # useful for getting some statistics on quark matching. No information is
        # used further.

        # Get a list of the 3 generated quarks that should correspond to a
        # top candidate
        ( tl_genquarks, tl_leptonicb ) = self.Get_tl_genquarks( event )

        # If there is an error in getting the quarks, the function returns 0
        if tl_genquarks == 0 and self.require_quark_match: return 0


        # Get list of btagged_jets

        tl_btagged_jets = []

        for (i_jet, jet) in enumerate(event.btagged_jets):
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM( jet.pt, jet.eta, jet.phi, jet.mass )
            setattr( x, 'i_jet', i_jet )
            tl_btagged_jets.append( x )


        # Get list of wquark_candidate_jets

        tl_wquark_candidate_jets = []

        for (i_jet, jet) in enumerate(event.wquark_candidate_jets):
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM( jet.pt, jet.eta, jet.phi, jet.mass )
            setattr( x, 'i_jet', i_jet )
            tl_wquark_candidate_jets.append( x )

        # Get list of subjets for the found httCandidate
        tl_subjets = []
        prefixes = [ 'sjW1', 'sjW2', 'sjNonW' ]
        for (i_subjet, prefix) in enumerate( prefixes ):
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM(
                getattr( top, prefix + 'pt' ),
                getattr( top, prefix + 'eta' ),
                getattr( top, prefix + 'phi' ),
                getattr( top, prefix + 'mass' ) )
            setattr( x, 'i_subjet', i_subjet )
            tl_subjets.append( x )

        # For printing purposes
        tl_subjets_backup = copy.deepcopy( tl_subjets )

        # If another top was found in the event, also get the subjets from that
        if other_top_present:
            tl_subjets_other_top = []
            prefixes = [ 'sjW1', 'sjW2', 'sjNonW' ]
            for (i_subjet, prefix) in enumerate( prefixes ):
                x = ROOT.TLorentzVector()
                x.SetPtEtaPhiM(
                    getattr( other_top, prefix + 'pt' ),
                    getattr( other_top, prefix + 'eta' ),
                    getattr( other_top, prefix + 'phi' ),
                    getattr( other_top, prefix + 'mass' ) )
                setattr( x, 'i_subjet', i_subjet )
                tl_subjets_other_top.append( x )


        
        # Create the attributes in the event class for the subjet case
        #   The subjets will be popped out of the *_jets_minus_sj lists        

        btagged_jets_minus_sj = copy.deepcopy( event.btagged_jets )
        btagged_subjet = []

        # Note: event.wquark_candidate_jets is of type 'set' for some reason
        wquark_candidate_jets_minus_sj = copy.deepcopy(
            list( event.wquark_candidate_jets ) )

        wquark_candidate_subjets = []


        ########################################
        # Perform combinatorics and calculate delR
        ########################################

        ### Match subjet to a GenBQuark and GenWZQuark

        # Copy the subjet list (in case a b-subjet needs to be popped)
        tl_subjets_copy = copy.deepcopy( tl_subjets )

        # Do the matching
        (i_sj, i_genquark ) = self.Link_smallest_delR(  tl_subjets_copy,
                                                        [tl_genquarks[0]] )

        if i_sj == 'No link':
            succes_bsj_bquark = False
        else:
            succes_bsj_bquark = True
            # Pop the subjet so it can't be matched again to a light quark
            tl_subjets_copy.pop(i_sj)

        # Do the matching
        ( i_sj1, i_sj2, i_q1, i_q2 ) = self.Link_2smallest_delR(
                                                        tl_subjets_copy,
                                                        tl_genquarks[1:] )

        if i_sj == 'No link':
            succes_lsj_lquark = False
        else:
            succes_lsj_lquark = True


        ### Match b-jets and w-jets to a GenBQuark and GenWZQuark

        # Do the matching
        (i_jet, i_genquark ) = self.Link_smallest_delR( tl_btagged_jets,
                                                        [tl_genquarks[0]] )

        if i_jet == 'No link':
            succes_bjet_bquark = False
        else:
            succes_bjet_bquark = True

        # Do the matching
        ( i_jet1, i_jet2, i_q1, i_q2 ) = self.Link_2smallest_delR(
                                                        tl_wquark_candidate_jets,
                                                        tl_genquarks[1:] )

        if i_jet1 == 'No link':
            succes_ljet_lquark = False
        else:
            succes_ljet_lquark = True


        ### Match a subjet to a btagged_jet

        # Do the matching
        (i_sj, i_bt_j ) = self.Link_smallest_delR( tl_subjets, tl_btagged_jets )

        # Treat potential errors        
        if i_sj == 'No link':
            if self.verbose: print 'Could not match a subjet to a btagged_jet'
            succes_sj_bjet = False

        else:
            succes_sj_bjet = True

            # Set the btagFlag and the TFs as attributes (needed in mem code)
            setattr(tl_subjets[i_sj],
                    'btagFlag',
                    btagged_jets_minus_sj[i_bt_j].btagFlag )

            setattr(tl_subjets[i_sj],
                    'tf_sj_b',
                    copy.deepcopy( btagged_jets_minus_sj[i_bt_j].tf_sj_b ) )

            setattr(tl_subjets[i_sj],
                    'tf_sj_l',
                    copy.deepcopy( btagged_jets_minus_sj[i_bt_j].tf_sj_l ) )

            # Remove the jet from btagged_jets_minus_sj that is the subjet
            btagged_jets_minus_sj.pop( i_bt_j )

            # Remove the b_tagged_subjet from the tl_subjets
            # CAREFUL: This is a TLorentzVector object, not a Jet object!
            btagged_subjet.append( tl_subjets.pop( i_sj ) )


        ### Match subjets to a wquark_candidate_jets

        # Do the matching
        ( i_sj1, i_sj2, i_wj1, i_wj2 ) = self.Link_2smallest_delR(
            tl_subjets, tl_wquark_candidate_jets )

        # Treat potential errors
        if i_sj1 == 'No link':
            if self.verbose:
                print 'Could not match a subjet to a wquark_candidate_jet'
            succes_sj_ljet = False

        else:
            succes_sj_ljet = True

            # First light quark match

            # Set the btagFlag and the TFs as attributes (needed in mem code)
            setattr(tl_subjets[i_sj1],
                    'btagFlag',
                    wquark_candidate_jets_minus_sj[i_wj1].btagFlag )

            setattr(tl_subjets[i_sj1],
                    'tf_sj_b',
                    copy.deepcopy(wquark_candidate_jets_minus_sj[i_wj1].tf_sj_b) )

            setattr(tl_subjets[i_sj1],
                    'tf_sj_l',
                    copy.deepcopy(wquark_candidate_jets_minus_sj[i_wj1].tf_sj_l) )

            # Pop the jet matched to a subjet from the list
            wquark_candidate_jets_minus_sj.pop( i_wj1 )

            # Fill in the appropiate subjet (this is a TLorentzVector object)
            wquark_candidate_subjets.append( tl_subjets.pop( i_sj1 ) )


            # Second light quark match

            # Set the btagFlag and the TFs as attributes (needed in mem code)
            setattr(tl_subjets[i_sj2],
                    'btagFlag',
                    wquark_candidate_jets_minus_sj[i_wj2].btagFlag )

            setattr(tl_subjets[i_sj2],
                    'tf_sj_b',
                    copy.deepcopy( wquark_candidate_jets_minus_sj[i_wj2].tf_sj_b ) )

            setattr(tl_subjets[i_sj2],
                    'tf_sj_l',
                    copy.deepcopy( wquark_candidate_jets_minus_sj[i_wj2].tf_sj_l ) )

            # Pop the jet matched to a subjet from the list
            wquark_candidate_jets_minus_sj.pop( i_wj2 )

            # Fill in the appropiate subjet (this is a TLorentzVector object)
            wquark_candidate_subjets.append( tl_subjets.pop( i_sj2 ) )

        
        # Fill the match tree

        self.match_tree_b.Get_child( succes_bjet_bquark ).n += 1

        self.match_tree_b.Get_child( succes_bjet_bquark,
                                     succes_bsj_bquark ).n += 1

        self.match_tree_b.Get_child( succes_bjet_bquark,
                                     succes_bsj_bquark,
                                     succes_sj_bjet, ).n += 1


        self.match_tree_l.Get_child( succes_ljet_lquark ).n += 1

        self.match_tree_l.Get_child( succes_ljet_lquark,
                                     succes_lsj_lquark ).n += 1

        self.match_tree_l.Get_child( succes_ljet_lquark,
                                     succes_lsj_lquark,
                                     succes_sj_ljet, ).n += 1


        """
        specific_print = (
            succes_bsj_bquark and \
            succes_lsj_lquark and \
            not succes_bjet_bquark and \
            succes_ljet_lquark and \
            succes_sj_bjet and \
            succes_sj_ljet )
        """

        #specific_print = other_top_present and succes_bjet_bquark and not succes_bsj_bquark
        specific_print = succes_bjet_bquark and not succes_bsj_bquark

        if specific_print:

            print '\nSpecific Print'
            print '=========='

            print 'Successful bjet to bquark, Unsuccessful bsubjet to bquark'

            print 'Input jets:'
            self.Print_particle_lists(
                ( event.good_leptons, 'Class', 'event.good_leptons'),
                ( event.btagged_jets, 'Class', 'event.btagged_jets'),
                ( event.wquark_candidate_jets, 'Class',
                    'event.wquark_candidate_jets'),
                )

            print 'Chosen httCandidate:'
            print 'delR with lepton: {0}'.format( top.delR_lepton )
            self.Print_particle_lists(
                ( tl_subjets_backup, 'TL', 'tl_subjets' ),
                )

            if other_top_present:
                print 'Another httCandidate was found:'
                print 'delR with lepton: {0}'.format( other_top.delR_lepton )
                self.Print_particle_lists(
                    ( tl_subjets_other_top, 'TL', 'tl_subjets_other_top' )
                    )

            print 'Quarks:'
            self.Print_particle_lists(
                ( [tl_genquarks[0]], 'TL', 'Hadronic b-quark'),
                ( [tl_leptonicb], 'TL', 'Leptonic b-quark'),
                ( tl_genquarks[1:], 'TL', 'Light quarks'),
                ( event.GenBQuarkFromH, 'Class', 'Higgs b-quark'),
                )

            print 'Output particles:'
            self.Print_particle_lists(

                ( btagged_jets_minus_sj, 'Class', 'btagged_jets_minus_sj'),
                ( btagged_subjet, 'TL', 'btagged_subjet'),

                ( wquark_candidate_jets_minus_sj, 'Class',
                    'wquark_candidate_jets_minus_sj' ),
                ( wquark_candidate_subjets, 'TL', 'wquark_candidate_subjets'),

                )

            print '=========='
            print 'End of Specific Print\n'


        # Only continue if ALL matches were successful
        pass_to_MEM = (
            succes_bsj_bquark and \
            succes_lsj_lquark and \
            succes_bjet_bquark and \
            succes_ljet_lquark and \
            succes_sj_bjet and \
            succes_sj_ljet )

        if not pass_to_MEM:
            return 0
        else:
            print 'Passing to MEM Analyzer'
            self.Statistics['n_passed_to_MEM'] += 1
            


        ########################################
        # Write the lists to the event class
        ########################################

        setattr(event,
                'btagged_jets_minus_sj',
                btagged_jets_minus_sj )

        setattr(event,
                'wquark_candidate_jets_minus_sj',
                wquark_candidate_jets_minus_sj )

        setattr(event,
                'btagged_subjet',
                btagged_subjet )

        setattr(event,
                'wquark_candidate_subjets',
                wquark_candidate_subjets )

        #print '\n    Printing found particles:'
        #self.Print_found_particles( event, tl_subjets_backup )

        #print '\n    Printing matched particles:'
        #self.Print_matched_particles( event )



    ########################################
    # Functions
    ########################################

    # Applies the cut criteria - returns True (survived) or False (did not survive)
    def Apply_Cut_criteria( self, candidate ):

        for ( attr, operator, cut_off ) in self.Cut_criteria:

            if not eval( '{0}{1}{2}'.format(
                getattr( candidate, attr ),
                operator,
                cut_off ) ):

                return False

        return True
    #--------------------------------------#

    # Sorts tops - criterium to be tested
    def SortTops( self, event, tops ):

        l = event.good_leptons[0]

        # Create TLorentzVector for lepton
        tl_lepton = ROOT.TLorentzVector()
        tl_lepton.SetPtEtaPhiM( l.pt, l.eta, l.phi, l.mass )

        # Create TLorentzVector for tops
        tl_tops = []
        for top in tops:
            tl_top = ROOT.TLorentzVector()
            tl_top.SetPtEtaPhiM( top.pt, top.eta, top.phi, top.mass )
            setattr( top, 'delR_lepton', tl_top.DeltaR(tl_lepton) )

        tops = sorted( tops, key=lambda x: -x.delR_lepton )

        return tops
    #--------------------------------------#


    # Prints a match tree
    def Print_Tree( self, t ):

        print '\nPrinting match tree {0}'.format( t.name )
        print '===================='
        print '  Total count           = {0}'.format( t.n )
        print '  Total passable to MEM = {0}\n'.format( t.Get_child(1,1,1).n )

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            'Jet to Qrk',
            '',
            'Subj to Qrk',
            '',
            'Jet to Subj',
            '',
            )

        print '|-------------+----------+-------------+----------+-------------+----------|'

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            'Success',
            str( t.Get_child(1).n ),
            'Success',
            str( t.Get_child(1,1).n ),
            'Success',
            str( t.Get_child(1,1,1).n ),
            )


        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            '',
            '',
            '',
            '',
            'Failed',
            str( t.Get_child(1,1,0).n ),
            )

        print '|             |          |-------------+----------+-------------+----------|'

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            '',
            '',
            'Failed',
            str( t.Get_child(1,0).n ),
            'Success',
            str( t.Get_child(1,0,1).n ),
            )

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            '',
            '',
            '',
            '',
            'Failed',
            str( t.Get_child(1,0,0).n ),
            )

        print '|-------------+----------+-------------+----------+-------------+----------|'

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            'Failed',
            str( t.Get_child(0).n ),
            'Success',
            str( t.Get_child(0,1).n ),
            'Success',
            str( t.Get_child(0,1,1).n ),
            )


        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            '',
            '',
            '',
            '',
            'Failed',
            str( t.Get_child(0,1,0).n ),
            )

        print '|             |          |-------------+----------+-------------+----------|'

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            '',
            '',
            'Failed',
            str( t.Get_child(0,0).n ),
            'Success',
            str( t.Get_child(0,0,1).n ),
            )

        print '|{0:13s}|{1:10s}|{2:13s}|{3:10s}|{4:13s}|{5:10s}|'.format(
            '',
            '',
            '',
            '',
            'Failed',
            str( t.Get_child(0,0,0).n ),
            )

        print '----------------------------------------------------------------------------'
    #--------------------------------------#

    
    # Print a non-predefined number of particle lists
    # Input is a list of tuples: input_args = [ (tuple), (tuple), ... ]
    # Every tuple is structured as follows:
    # (tuple) = ( particle list, TLorentz/Class label, name of the particle )
    def Print_particle_lists( self, *input_args ):

        for ( p_list, mode, name ) in input_args:

            print '\n    Printing {0}:'.format( name )

            if mode == 'TL':

                for q in p_list:
                    print '    [ {0:9s} | {1:9s} | {2:9s} | {3:9s} ]'.format(
                        '{0:.3f}'.format(q.Pt()),
                        '{0:.3f}'.format(q.Eta()),
                        '{0:.3f}'.format(q.Phi()),
                        '{0:.3f}'.format(q.M()) )

            if mode == 'Class':

                for q in p_list:
                    print '    [ {0:9s} | {1:9s} | {2:9s} | {3:9s} ]'.format(
                        '{0:.3f}'.format(q.pt),
                        '{0:.3f}'.format(q.eta),
                        '{0:.3f}'.format(q.phi),
                        '{0:.3f}'.format(q.mass) )

        print ''
    #--------------------------------------#

    # Deletes duplicate WZ Quarks
    def CompareWZQuarks(self, event ):

        if len(event.GenWZQuark) < 4:
            return 0

        quarks = event.GenWZQuark

        Is_Duplicate = ( quarks[-1].pt==quarks[1].pt and \
            quarks[-2].pt==quarks[0].pt )

        if Is_Duplicate:
            event.GenWZQuark = event.GenWZQuark[0:-2]

        return 0
    #--------------------------------------#

    # Gets a list of 3 quarks
    #  - The first quark is a Gen B quark. This should be the hadronic B quark.
    #    Which quark is hadronic is determined by adding the light quarks to the
    #    B quarks, and seeing which combined mass comes closer to the top mass
    #  - Output looks like: [ BQuark, lightQuark1, lightQuark2 ], where the list
    #    entries are TLorentzVector objects. 
    def Get_tl_genquarks(self, event ):

        # Check if right amount of quarks was generated
        if len(event.GenWZQuark)<2:
            self.Statistics['n_too_few_WZ'] += 1
            return (0,0)
        if len(event.GenBQuarkFromTop)<2:
            self.Statistics['n_too_few_B'] += 1
            return (0,0)
        elif len(event.GenWZQuark)>2:
            self.Statistics['n_too_many_WZ'] += 1
            return (0,0)
        elif len(event.GenBQuarkFromTop)>2:
            self.Statistics['n_too_many_B'] += 1
            return (0,0)

        self.Statistics['n_right_nr_of_quarks'] += 1

        # Make list of TLorentzVector objects for the 2 light quarks
        tl_GenWZQuarks = []
        for l in event.GenWZQuark:
            tl_l = ROOT.TLorentzVector()
            tl_l.SetPtEtaPhiM( l.pt, l.eta, l.phi, l.mass )
            tl_GenWZQuarks.append( tl_l )

        # Make list for the 2 B quarks, and a list of B quarks + light quarks
        tl_GenBQuarks = []
        tl_Combined = []
        for (b_i, b) in enumerate(event.GenBQuarkFromTop):

            tl_b = ROOT.TLorentzVector()
            tl_b.SetPtEtaPhiM( b.pt, b.eta, b.phi, b.mass )

            tl_GenBQuarks.append( tl_b )
            tl_Combined.append( tl_b + tl_GenWZQuarks[0] + tl_GenWZQuarks[1] )

        # Calculate mass difference from top mass
        delmass0 = abs(tl_Combined[0].M() - self.top_mass)
        delmass1 = abs(tl_Combined[1].M() - self.top_mass)

        # Save the mass difference with top mass
        setattr( event.GenBQuarkFromTop[0], 'delmass_top', delmass0 )
        setattr( event.GenBQuarkFromTop[1], 'delmass_top', delmass1 )

        # Make sure the B quark with lowest del mass to top mass is at index 0
        # (both for the tl list and in the event)
        if delmass1 < delmass0:

            tl_GenBQuarks = [ tl_GenBQuarks[1], tl_GenBQuarks[0] ]

            event.GenBQuarkFromTop = [
                event.GenBQuarkFromTop[1],
                event.GenBQuarkFromTop[0] ]

        setattr( tl_GenBQuarks[0], 'is_hadr', 1 )
        setattr( event.GenBQuarkFromTop[0], 'is_hadr', 1 )
        setattr( event.GenBQuarkFromTop[1], 'is_hadr', 0 )

        # Create the definitive list of 3 quarks
        tl_GenQuarks = []
        tl_GenQuarks.append( tl_GenBQuarks[0] )

        # There should be only 2 generated WZQuarks:
        tl_GenQuarks.append( tl_GenWZQuarks[0] )
        tl_GenQuarks.append( tl_GenWZQuarks[1] )
        
        return ( tl_GenQuarks, tl_GenBQuarks[1] )
    #--------------------------------------#

    # Simple algorithm that matches the smallest delta R for two lists of TL vectors
    def Link_smallest_delR( self, tl_genquarks, tl_jets ):

        n_jets = len(tl_jets)
        n_quarks = len(tl_genquarks)

        Rmat = [[ (tl_genquarks[i].DeltaR( tl_jets[j] )) \
            for j in range(n_jets)] for i in range(n_quarks) ]

        Rmin = 9999.0
        
        for (r, row) in enumerate(Rmat):
            for (c, ele) in enumerate(row):
                if ele < Rmin and ele < self.R_cut:
                    Rmin = ele
                    r_min = r
                    c_min = c

        if Rmin == 9999.0: return ( 'No link', 0)

        return (r_min, c_min)


    def Link_2smallest_delR( self, tl_quarks, tl_jets ):

        ( i_q1, i_j1 ) = self.Link_smallest_delR( tl_quarks, tl_jets )

        if i_q1 == 'No link': return ( 'No link', 0, 0, 0 )

        # Essentially a '.pop' operation, without destroying the initial list
        sec_tl_quarks = tl_quarks[:i_q1]
        sec_tl_quarks.extend( tl_quarks[i_q1+1:] )

        sec_tl_jets = tl_jets[:i_j1]
        sec_tl_jets.extend( tl_jets[i_j1+1:] )

        ( i_q2, i_j2 ) = self.Link_smallest_delR( sec_tl_quarks, sec_tl_jets )

        if i_q2 == 'No link': return ( 'No link', 0, 0, 0 )

        return ( i_q1, i_q2, i_j1, i_j2 )


    def Link_3smallest_delR( self, tl_quarks_orig, tl_jets_orig ):

        tl_quarks = copy.deepcopy( tl_quarks_orig )
        tl_jets = copy.deepcopy( tl_jets_orig )

        ( i_q1, i_j1 ) = self.Link_smallest_delR( tl_quarks, tl_jets )
        if i_q1 == 'No link': return ( 'No link', 1, 0, 0, 0, 0)
        tl_quarks.pop(i_q1)
        tl_jets.pop(i_j1)

        ( i_q2, i_j2 ) = self.Link_smallest_delR( tl_quarks, tl_jets )
        if i_q2 == 'No link': return ( 'No link', 2, 0, 0, 0, 0 )
        tl_quarks.pop(i_q2)
        tl_jets.pop(i_j2)

        ( i_q3, i_j3 ) = self.Link_smallest_delR( tl_quarks, tl_jets )
        if i_q3 == 'No link': return ( 'No link', 3, 0, 0, 0, 0 )

        return ( i_q1, i_q2, i_q3, i_j1, i_j2, i_j3 )


    # Matches quarks with specified jets
    # - The function just takes two lists of TLorentzVector objects and matches
    #   elements from both lists based on minimum SUMMED delta R value.
    # - The tl_genquarks list must have 3 elements. Either the list of quarks or
    #   the list of subjets can be inserted here.
    def Do_delR_combinatorics( self, tl_genquarks, tl_jets ):

        n_jets = len(tl_jets)
        n_quarks = len(tl_genquarks)

        # Create delR matrix:
        Rmat = [[ (tl_genquarks[i].DeltaR( tl_jets[j] )) \
            for j in range(n_jets)] for i in range(n_quarks) ]

        """        
        print '\ndelR matrix:'
        for row in Rmat:
            for j in row:
                sys.stdout.write( '{0:.5f}'.format(j) + ' ')
            print ''
        """

        links_per_quark = [ [] for i in range(n_quarks) ]

        for i in range(n_quarks):
            for j in range(n_jets):
                if Rmat[i][j] < self.R_cut:
                    links_per_quark[i].append(j)
        
        
        # Perform some checks: see if there are 2 quarks linked to only 1 jet, and
        # check whether every quark has at least 1 jet it can be linked to

        for i in range(n_quarks):

            if len(links_per_quark[i]) == 0:
                return ('quark_has_no_jet',0)

            # Check if 2 quarks can only be linked to the same jet
            for j in range(n_quarks):
                if i==j: continue
                
                if len(links_per_quark[i])==1 and len(links_per_quark[j])==1 \
                    and links_per_quark[i]==links_per_quark[j]:
                    return ('2_quarks_need_1_jet',0)


        # Combinatorics: find the lowest sum of delR values
        
        sumR = 100000.0
        final_links = 0

        for (i_ind, i) in enumerate(links_per_quark[0]):
            for (j_ind, j) in enumerate(links_per_quark[1]):
                for (k_ind, k) in enumerate(links_per_quark[2]):

                    # Check if i,j,k gets a unique combination
                    if len( set( [ i, j, k ] ) ) == 3:

                        # Replace the sumR if it's smaller than the previous minimum
                        if Rmat[0][i] + Rmat[1][j] + Rmat[2][k] < sumR:
                            sumR = Rmat[0][i] + Rmat[1][j] + Rmat[2][k]
                            final_links = [ i, j, k ]

        if final_links == 0:
            print links_per_quark
            print 'No unique permutation found'
            return ('no_unique_permutation',0)

        final_delR_list = [ Rmat[i][final_links[i]] for i in range(n_quarks) ]

        return ( final_links , final_delR_list )
    #--------------------------------------#




class EventIDFilterAnalyzer(FilterAnalyzer):
    """
    """

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(EventIDFilterAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf
        self.event_whitelist = self.conf.general.get("eventWhitelist", None)

    def beginLoop(self, setup):
        super(EventIDFilterAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        passes = True
        if not self.event_whitelist is None:
            passes = False
            if (event.input.run, event.input.lumi, event.input.evt) in self.event_whitelist:
                print "IDFilter", (event.input.run, event.input.lumi, event.input.evt)
                passes = True

        if passes and "eventboundary" in self.conf.general["verbosity"]:
            print "---", event.input.run, event.input.lumi, event.input.evt
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class LeptonAnalyzer(FilterAnalyzer):
    """
    Analyzes leptons and applies single-lepton and di-lepton selection.

    Relies on TTH.MEAnalysis.VHbbTree.EventAnalyzer for inputs.

    Configuration:
    Conf.leptons[channel][cuttype] where channel=mu,ele, cuttype=tight,loose,(+veto)
    the lepton cuts must specify pt, eta and isolation cuts.

    Returns:
    event.good_leptons (list of VHbbTree.selLeptons): contains the leptons that pass the SL XOR DL selection.
        Leptons are ordered by flavour and pt.
    event.is_sl, is_dl (bool): specifies if the event passes SL or DL selection.

    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(LeptonAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf

    def beginLoop(self, setup):
        super(LeptonAnalyzer, self).beginLoop(setup)

        self.counters["processing"].register("sl")
        self.counters["processing"].register("dl")
        self.counters["processing"].register("slanddl")

        self.counters.addCounter("leptons")
        self.counters["leptons"].register("any")
        for l in ["mu", "el"]:
            for a in ["tight", "loose"]:
                for b in ["", "_veto"]:
                    lt = l + "_" + a + b
                    self.counters["leptons"].register(lt)

    def process(self, event):
        self.counters["processing"].inc("processed")
        self.counters["leptons"].inc("any", len(event.selLeptons))

        event.mu = filter(
            lambda x: abs(x.pdgId) == 13,
            event.selLeptons,
        )
        event.el = filter(
            lambda x: abs(x.pdgId) == 11,
            event.selLeptons,
        )

        for a in ["tight", "loose"]:
            for b in ["", "_veto"]:
                sumleps = []
                for l in ["mu", "el"]:
                    lepcuts = self.conf.leptons[l][a+b]
                    incoll = getattr(event, l)

                    leps = filter(
                        lambda x: (
                            x.pt > lepcuts["pt"]
                            and abs(x.eta) < lepcuts["eta"]
                            and abs(getattr(x, self.conf.leptons[l]["isotype"])) < lepcuts["iso"]
                        ), incoll
                    )

                    if b == "_veto":
                        good = getattr(event, "{0}_{1}".format(l, a))
                        leps = filter(lambda x: x not in good, leps)

                    if a == "tight":
                        leps = filter(
                            lambda x: x.tightId,
                            leps
                        )
                    elif a == "loose":
                        leps = filter(
                            lambda x: x.looseIdPOG,
                            leps
                        )
                    lep = sorted(leps, key=lambda x: x.pt, reverse=True)
                    sumleps += leps
                    lt = l + "_" + a + b
                    setattr(event, lt, leps)
                    setattr(event, "n_"+lt, len(leps))
                    self.counters["leptons"].inc(lt, len(leps))

                setattr(event, "lep_{0}".format(a+b), sumleps)
                setattr(event, "n_lep_{0}".format(a+b), len(sumleps))


        event.is_sl = (event.n_lep_tight == 1 and event.n_lep_tight_veto == 0)
        event.is_dl = (event.n_lep_loose == 2 and event.n_lep_loose_veto == 0)

        if event.is_sl:
            self.counters["processing"].inc("sl")
            event.good_leptons = event.mu_tight + event.el_tight
        if event.is_dl:
            self.counters["processing"].inc("dl")
            event.good_leptons = event.mu_loose + event.el_loose

        passes = event.is_sl or event.is_dl
        if event.is_sl and event.is_dl:
            #print "pathological SL && DL event: {0}".format(event)
            self.counters["processing"].inc("slanddl")
            passes = False

        if passes:
            self.counters["processing"].inc("passes")
        return passes



class JetAnalyzer(FilterAnalyzer):
    """
    Performs jet selection and b-tag counting.
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(JetAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf

    def beginLoop(self, setup):
        super(JetAnalyzer, self).beginLoop(setup)
        self.counters.addCounter("jets")
        self.counters["jets"].register("any")
        self.counters["jets"].register("good")
        for (btag_wp_name, btag_wp) in self.conf.jets["btagWPs"].items():
            self.counters["jets"].register(btag_wp_name)


    def process(self, event):
        self.counters["processing"].inc("processed")
        self.counters["jets"].inc("any", len(event.Jet))

        #pt-descending input jets
        if "input" in self.conf.general["verbosity"]:
            for j in event.Jet:
                print "ijet", j.pt, j.eta, j.phi, j.mass, j.btagCSV, j.mcFlavour

        event.good_jets = sorted(
            filter(
                lambda x: (
                    x.pt > self.conf.jets["pt"]
                    and abs(x.eta) < self.conf.jets["eta"]
                ), event.Jet
            ),
            key=lambda x: x.pt, reverse=True
        )


        #Assing jet transfer functions
        for jet in event.good_jets:
            jet_eta_bin = 0
            if abs(jet.eta)>1.0:
                jet_eta_bin = 1

            #If True, TF [0] - reco, x - gen
            #If False, TF [0] - gen, x - reco
            eval_gen = False
            jet.tf_b = self.conf.tf_matrix['b'][jet_eta_bin].Make_Formula(eval_gen)
            jet.tf_l = self.conf.tf_matrix['l'][jet_eta_bin].Make_Formula(eval_gen)
            #If [0] - gen, x - pt cutoff
            jet.tf_b_lost = self.conf.tf_matrix['b'][jet_eta_bin].Make_CDF()
            jet.tf_l_lost = self.conf.tf_matrix['l'][jet_eta_bin].Make_CDF()

            #Set jet pt threshold for CDF
            jet.tf_b_lost.SetParameter(0, self.conf.jets["pt"])
            jet.tf_l_lost.SetParameter(0, self.conf.jets["pt"])


            # Transfer functions for subjets

            #If True, TF [0] - reco, x - gen
            #If False, TF [0] - gen, x - reco
            eval_gen = False
            jet.tf_sj_b = self.conf.tf_sj_matrix['b'][jet_eta_bin].Make_Formula(eval_gen)
            jet.tf_sj_l = self.conf.tf_sj_matrix['l'][jet_eta_bin].Make_Formula(eval_gen)
            #If [0] - gen, x - pt cutoff
            jet.tf_sj_b_lost = self.conf.tf_sj_matrix['b'][jet_eta_bin].Make_CDF()
            jet.tf_sj_l_lost = self.conf.tf_sj_matrix['l'][jet_eta_bin].Make_CDF()

            #Set jet pt threshold for CDF
            jet.tf_sj_b_lost.SetParameter(0, self.conf.jets["pt"])
            jet.tf_sj_l_lost.SetParameter(0, self.conf.jets["pt"])



        event.numJets = len(event.good_jets)
        self.counters["jets"].inc("good", len(event.good_jets))

        event.btagged_jets_bdisc = {}
        event.buntagged_jets_bdisc = {}
        for (btag_wp_name, btag_wp) in self.conf.jets["btagWPs"].items():
            algo, wp = btag_wp
            event.btagged_jets_bdisc[btag_wp_name] = filter(
                lambda x: getattr(x, algo) > wp,
                event.good_jets
            )
            event.buntagged_jets_bdisc[btag_wp_name] = filter(
                lambda x: getattr(x, algo) <= wp,
                event.good_jets
            )
            self.counters["jets"].inc(btag_wp_name,
                len(event.btagged_jets_bdisc[btag_wp_name])
            )
            setattr(event, "nB"+btag_wp_name, len(event.btagged_jets_bdisc[btag_wp_name]))

        #Find jets that pass/fail the specified default b-tagging algo/working point
        event.buntagged_jets_bdisc = event.buntagged_jets_bdisc[self.conf.jets["btagWP"]]
        event.btagged_jets_bdisc = event.btagged_jets_bdisc[self.conf.jets["btagWP"]]

        #Find how many of these tagged jets are actually true b jets
        event.n_tagwp_tagged_true_bjets = 0
        for j in event.btagged_jets_bdisc:
            if abs(j.mcFlavour) == 5:
                event.n_tagwp_tagged_true_bjets += 1

        #Require at least 4 good jets in order to continue analysis
        passes = len(event.good_jets) >= 4
        if passes:
            self.counters["processing"].inc("passes")

        corrMet_px = event.met[0].px
        corrMet_py = event.met[0].py
        sum_dEx = 0
        sum_dEy = 0
        for jet in event.good_jets:
            Prec = lvec(jet)
            Pgen = lvec(jet)
            Pgen.SetPtEtaPhiM(jet.mcPt, jet.mcEta, jet.mcPhi, jet.mcM)
            Erec = Prec.E()
            Egen = Pgen.E()
            dEx = (Erec-Egen) * Prec.Px()/Prec.P()
            dEy = (Erec-Egen) * Prec.Py()/Prec.P()
            #print Erec, Egen
            sum_dEx += dEx
            sum_dEy += dEy
        corrMet_px += sum_dEx
        corrMet_py += sum_dEy
        #print (sum_dEx, sum_dEy), (corrMet_px, event.met[0].px), (corrMet_py, event.met[0].py)
        event.met_jetcorr = [MET(px=corrMet_px, py=corrMet_py)]
        event.met_gen = [MET(pt=event.met[0].genPt, phi=event.met[0].genPhi)]

        return passes


class BTagLRAnalyzer(FilterAnalyzer):
    """
    Performs b-tag likelihood ratio calculations
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(BTagLRAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf
        self.bTagAlgo = self.conf.jets["btagAlgo"]
        self.cplots_old = ROOT.TFile(self.conf.general["controlPlotsFileOld"])
        self.cplots = ROOT.TFile(self.conf.general["controlPlotsFile"])

        cplots_fmt = self.conf.general.get("controlPlotsFormat", "8tev")
        self.csv_pdfs_old = {
        }
        for x in ["b", "c", "l"]:
            for b in ["Bin0", "Bin1"]:
                self.csv_pdfs_old[(x, b)] = self.cplots_old.Get(
                    "csv_{0}_{1}__csv_rec".format(x, b)
                )

        self.csv_pdfs = {
        }
        for x in ["b", "c", "l"]:
            for b in ["Bin0", "Bin1"]:
                self.csv_pdfs[(x, b)] = self.cplots.Get(
                    "csv_{0}_{1}__csv_rec".format(x, b)
                )
                self.csv_pdfs[(x, b)].Scale(1.0 / self.csv_pdfs[(x, b)].Integral())
            self.csv_pdfs[(x, "pt_eta")] = self.cplots.Get(
                "csv_{0}_pt_eta".format(x)
            )
            self.csv_pdfs[(x, "pt_eta")].Scale(1.0 / self.csv_pdfs[(x, "pt_eta")].Integral())

    def get_pdf_prob(self, flavour, pt, eta, csv, kind):

        _bin = "Bin1" if abs(eta)>1.0 else "Bin0"

        if kind == "old":
            h = self.csv_pdfs_old[(flavour, _bin)]
        elif kind == "new_eta_1bin":
            h = self.csv_pdfs[(flavour, _bin)]
        elif kind == "new_pt_eta_bin_3d":
            h = self.csv_pdfs[(flavour, "pt_eta")]

        assert h != None, "flavour={0} kind={1}".format(flavour, kind)

        if csv < 0:
            csv = 0.0
        if csv > 1.0:
            csv = 1.0

        if kind == "old" or kind == "new_eta_1bin":
            nb = h.FindBin(csv)
            #if csv = 1 -> goes into overflow and pdf = 0.0
            #as a solution, take the next-to-last bin
            if nb >= h.GetNbinsX():
                nb = nb - 1
            ret = h.GetBinContent(nb)
        elif kind == "new_pt_eta_bin_3d":
            nb = h.FindBin(pt, abs(eta), csv)
            ret = h.GetBinContent(nb)
        return ret

    def beginLoop(self, setup):
        super(BTagLRAnalyzer, self).beginLoop(setup)

    def evaluate_jet_prob(self, pt, eta, csv, kind):
        return (
            self.get_pdf_prob("b", pt, eta, csv, kind),
            self.get_pdf_prob("c", pt, eta, csv, kind),
            self.get_pdf_prob("l", pt, eta, csv, kind)
        )

    def btag_likelihood(self, probs, nB, nC):

        perms = itertools.permutations(range(len(probs)))

        P = 0.0
        max_p = -1.0
        nperms = 0
        best_perm = None

        np = len(probs)
        for perm in perms:
            p = 1.0

            for i in range(0, nB):
                p *= probs[perm[i]][0]
            for i in range(nB, min(nB + nC, np)):
                p *= probs[perm[i]][1]
            for i in range(nB + nC, np):
                p *= probs[perm[i]][2]

            #print nperms, p, perm, max_p, best_perm
            if p > max_p:
                best_perm = perm
                max_p = p

            P += p
            nperms += 1
        P = P / float(nperms)
        assert nperms > 0
        return P, best_perm
        #end permutation loop

    def process(self, event):
        self.counters["processing"].inc("processed")

        #Take first 6 most b-tagged jets for btag LR
        jets_for_btag_lr = sorted(
            event.good_jets,
            key=lambda x: getattr(x, self.bTagAlgo),
            reverse=True,
        )[0:6]

        jet_probs = {
            kind: [
                self.evaluate_jet_prob(j.pt, j.eta, getattr(j, self.bTagAlgo), kind)
                for j in jets_for_btag_lr
            ]
            for kind in [
            "old", "new_eta_1bin",
            "new_pt_eta_bin_3d"
            ]
        }

        # for nj, j in enumerate(jets_for_btag_lr):
        #     print j.btagCSV, j.mcFlavour, jet_probs["old"][nj], jet_probs["new_eta_1bin"][nj], jet_probs["new_pt_eta_bin_3d"][nj]
        #
        jet_csvs = [
            getattr(j, self.bTagAlgo)
            for j in event.good_jets
        ]

        ph = None
        best_4b_perm = 0
        best_2b_perm = 0
        event.btag_lr_4b_old, ph = self.btag_likelihood(jet_probs["old"], 4, 0)
        event.btag_lr_2b_old, v = self.btag_likelihood(jet_probs["old"], 2, 0)

        event.btag_lr_4b, best_4b_perm = self.btag_likelihood(jet_probs["new_eta_1bin"], 4, 0)
        event.btag_lr_4b_1c, ph = self.btag_likelihood(jet_probs["new_eta_1bin"], 4, 1)
        event.btag_lr_2b_2c, ph = self.btag_likelihood(jet_probs["new_eta_1bin"], 2, 2)

        event.btag_lr_2b, best_2b_perm = self.btag_likelihood(jet_probs["new_eta_1bin"], 2, 1)
        event.btag_lr_2b_1c, best_2b_perm = self.btag_likelihood(jet_probs["new_eta_1bin"], 2, 1)

        event.btag_lr_4b_alt, best_4b_perm_alt = self.btag_likelihood(jet_probs["new_pt_eta_bin_3d"], 4, 0)
        event.btag_lr_2b_alt, best_2b_perm_alt = self.btag_likelihood(jet_probs["new_pt_eta_bin_3d"], 2, 0)

        def lratio(l1, l2):
            if l1+l2>0:
                return l1/(l1+l2)
            else:
                return 0.0

        event.btag_LR_4b_2b_old = lratio(event.btag_lr_4b_old, event.btag_lr_2b_old)
        event.btag_LR_4b_2b = lratio(event.btag_lr_4b, event.btag_lr_2b)
        event.btag_LR_4b_2b_alt = lratio(event.btag_lr_4b_alt, event.btag_lr_2b_alt)
        #event.btag_LR_4b_2b_alt = 0

        event.buntagged_jets_by_LR_4b_2b = [jets_for_btag_lr[i] for i in best_4b_perm[4:]]
        event.btagged_jets_by_LR_4b_2b = [jets_for_btag_lr[i] for i in best_4b_perm[0:4]]

        for i in range(len(event.good_jets)):
            event.good_jets[i].btagFlag = 0.0

        #Jets are untagged according to the b-tagging likelihood ratio permutation
        if self.conf.jets["untaggedSelection"] == "btagLR":
            event.buntagged_jets = event.buntagged_jets_by_LR_4b_2b
            event.btagged_jets = event.btagged_jets_by_LR_4b_2b
        #Jets are untagged according to b-discriminatr
        elif self.conf.jets["untaggedSelection"] == "btagCSV":
            event.buntagged_jets = event.buntagged_jets_bdisc
            event.btagged_jets = event.btagged_jets_bdisc

        #Take first 4 most b-tagged jets
        btagged = sorted(event.btagged_jets, key=lambda x: x.btagCSV, reverse=True)[0:4]
        #Set these jets to be used as b-quarks in the MEM
        #We don't want to use more than 4 b-quarks in the hypothesis
        for jet in btagged:
            idx = event.good_jets.index(jet)
            event.good_jets[idx].btagFlag = 1.0

        passes = True
        if passes:
            self.counters["processing"].inc("passes")

        return passes

class MECategoryAnalyzer(FilterAnalyzer):
    """
    Performs ME categorization
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(MECategoryAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.cat_map = {"NOCAT":-1, "cat1": 1, "cat2": 2, "cat3": 3, "cat6":6}
        self.btag_cat_map = {"NOCAT":-1, "L": 0, "H": 1}

    def beginLoop(self, setup):
        super(MECategoryAnalyzer, self).beginLoop(setup)

        for c in ["NOCAT", "cat1", "cat2", "cat3", "cat6"]:
            self.counters["processing"].register(c)

    def process(self, event):
        self.counters["processing"].inc("processed")

        cat = "NOCAT"
        pass_btag_lr = (self.conf.jets["untaggedSelection"] == "btagLR" and
            event.btag_LR_4b_2b > self.conf.mem["btagLRCut"][event.cat]
        )
        pass_btag_csv = (self.conf.jets["untaggedSelection"] == "btagCSV" and
            len(event.btagged_jets) >= 4
        )
        cat_btag = "NOCAT"

        if pass_btag_lr or pass_btag_csv:
            cat_btag = "H"

        if event.is_sl:

            #at least 6 jets, if 6, Wtag in [60,100], if more Wtag in [72,94]
            if ((len(event.good_jets) == 6 and event.Wmass >= 60 and event.Wmass < 100) or
               (len(event.good_jets) > 6 and event.Wmass >= 72 and event.Wmass < 94)):
               cat = "cat1"
               #W-tagger fills wquark_candidate_jets
            #at least 6 jets, no W-tag
            elif len(event.good_jets) >= 6:
                cat = "cat2"
            #one W daughter missing
            elif len(event.good_jets) == 5:
                event.wquark_candidate_jets = event.buntagged_jets
                cat = "cat3"
        elif event.is_dl and len(event.good_jets)>=4:
            event.wquark_candidate_jets = []
            cat = "cat6"

        self.counters["processing"].inc(cat)
        event.cat = cat
        event.cat_btag = cat_btag
        event.catn = self.cat_map.get(cat, -1)
        event.cat_btag_n = self.btag_cat_map.get(cat_btag, -1)

        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class WTagAnalyzer(FilterAnalyzer):
    """
    Performs W-mass calculation on pairs of untagged jets.

    Jets are considered untagged according to the b-tagging permutation which
    gives the highest likelihood of the event being a 4b+Nlight event.
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(WTagAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

    def beginLoop(self, setup):
        super(WTagAnalyzer, self).beginLoop(setup)

    def pair_mass(self, j1, j2):
        """
        Calculates the invariant mass of a two-particle system.
        """
        lv1, lv2 = [lvec(j) for j in [j1, j2]]
        tot = lv1 + lv2
        return tot.M()

    def find_best_pair(self, jets):
        """
        Finds the pair of jets whose invariant mass is closest to mW=80 GeV.
        Returns the sorted vector of [(mass, jet1, jet2)], best first.
        """
        ms = []

        #Keep track of index pairs already calculated
        done_pairs = set([])

        #Double loop over all jets
        for i in range(len(jets)):
            for j in range(len(jets)):

                #Make sure we haven't calculated this index pair yet
                if (i,j) not in done_pairs and i!=j:
                    m = self.pair_mass(jets[i], jets[j])
                    ms += [(m, jets[i], jets[j])]

                    #M(i,j) is symmetric, hence add both pairs
                    done_pairs.add((i,j))
                    done_pairs.add((j,i))
        ms = sorted(ms, key=lambda x: abs(x[0] - 80.0))
        return ms

    def process(self, event):
        self.counters["processing"].inc("processed")

        event.Wmass = 0.0

        #we keep a set of the Q quark candidate jets
        event.wquark_candidate_jets = set([])

        #Need at least 2 untagged jets to calculate W mass
        if len(event.buntagged_jets)>=2:
            bpair = self.find_best_pair(event.buntagged_jets)

            #Get the best mass
            event.Wmass = bpair[0][0]

            #All masses
            event.Wmasses = [bpair[i][0] for i in range(len(bpair))]

            #Add at most 2 best pairs two W quark candidates
            for i in range(min(len(bpair), 2)):
                event.wquark_candidate_jets.add(bpair[i][1])
                event.wquark_candidate_jets.add(bpair[i][2])

            if "reco" in self.conf.general["verbosity"]:
                print "Wmass", event.Wmass, event.good_jets.index(bpair[1]), event.good_jets.index(bpair[2])

        #If we can't calculate W mass, untagged jets become the candidate
        else:
            for jet in event.buntagged_jets:
                event.wquark_candidate_jets.add(jet)


        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class GenRadiationModeAnalyzer(FilterAnalyzer):
    """
    Performs B/C counting in order to classify heavy flavour / light flavour events.

    We count the number of reconstructed jets which are matched to b/c quarks by CMSSW (ghost clustering).
    From this, the jets matched to b quarks from tops are subtracted.

    Therefore, nMatchSimB == 2 corresponds to 2 additional gluon radiation b quarks
    which are reconstructed as good jets.
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(GenRadiationModeAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

    def beginLoop(self, setup):
        super(GenRadiationModeAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        event.nMatchSimB = 0
        event.nMatchSimC = 0
        lv_bs = map(lvec, event.GenBQuarkFromTop)
        for jet in event.good_jets:
            lv_j = lvec(jet)

            if (lv_j.Pt() > 20 and abs(lv_j.Eta()) < 2.5):
                if any([lv_b.DeltaR(lv_j) < 0.5 for lv_b in lv_bs]):
                    continue
                absid = abs(jet.mcFlavour)
                if absid == 5:
                    event.nMatchSimB += 1
                if absid == 4:
                    event.nMatchSimC += 1

        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class GenTTHAnalyzer(FilterAnalyzer):
    """
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(GenTTHAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

    def beginLoop(self, setup):
        super(GenTTHAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        #Somehow, the GenWZQuark distribution is duplicated
        event.l_quarks_w = event.GenWZQuark[0:len(event.GenWZQuark)/2]
        event.b_quarks_t = event.GenBQuarkFromTop
        event.b_quarks_h = event.GenBQuarkFromH
        event.lep_top = event.GenLepFromTop
        event.nu_top = event.GenNuFromTop

        event.cat_gen = None
        event.n_cat_gen = -1

        if (len(event.lep_top) == 1 and
            len(event.nu_top) == 1 and
            len(event.l_quarks_w) == 2 and
            len(event.b_quarks_t) == 2):
            event.cat_gen = "sl"
            event.n_cat_gen = 0
        elif (len(event.lep_top) == 2 and
            len(event.nu_top) == 2 and
            len(event.l_quarks_w) == 0 and
            len(event.b_quarks_t) == 2):
            event.cat_gen = "dl"
            event.n_cat_gen = 1
        elif (len(event.lep_top) == 0 and
            len(event.nu_top) == 0 and
            len(event.l_quarks_w) == 4 and
            len(event.b_quarks_t) == 2):
            event.cat_gen = "fh"
            event.n_cat_gen = 2

        #Get the total MET from the neutrinos
        spx = 0
        spy = 0
        for nu in event.nu_top:
            p4 = lvec(nu)
            spx += p4.Px()
            spy += p4.Py()
        event.tt_met = [MET(px=spx, py=spy)]


        #Get the total ttH visible pt at gen level
        spx = 0
        spy = 0
        for p in (event.l_quarks_w + event.b_quarks_t +
            event.b_quarks_h + event.lep_top):
            p4 = lvec(p)
            spx += p4.Px()
            spy += p4.Py()
        event.tth_px_gen = spx
        event.tth_py_gen = spy
        
        #Calculate tth recoil
        #rho = -met - tth_matched
        event.tth_rho_px_gen = -event.met_gen[0].px - event.tth_px_gen
        event.tth_rho_py_gen = -event.met_gen[0].py - event.tth_py_gen


        if "gen" in self.conf.general["verbosity"]:
            for j in event.l_quarks_w:
                print "q(W)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.b_quarks_t:
                print "b(t)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.lep_top:
                print "l(t)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.nu_top:
                print "n(t)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.b_quarks_h:
                print "b(h)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            print "gen cat", event.cat_gen, event.n_cat_gen

        #Store for each jet, specified by it's index in the jet
        #vector, if it is matched to any gen-level quarks
        matched_pairs = {}

        def match_jets_to_quarks(jetcoll, quarkcoll, label, label_numeric):
            for ij, j in enumerate(jetcoll):
                for iq, q in enumerate(quarkcoll):
                    l1 = lvec(q)
                    l2 = lvec(j)
                    dr = l1.DeltaR(l2)
                    if dr < 0.3:
                        if matched_pairs.has_key(ij):
                            if matched_pairs[ij][1] > dr:
                                matched_pairs[ij] = (label, iq, dr, label_numeric)
                        else:
                            matched_pairs[ij] = (label, iq, dr, label_numeric)
        #print "GEN", len(event.GenWZQuark), len(event.GenBQuarkFromTop), len(event.GenBQuarkFromH)
        match_jets_to_quarks(event.good_jets, event.l_quarks_w, "wq", 0)
        match_jets_to_quarks(event.good_jets, event.b_quarks_t, "tb", 1)
        match_jets_to_quarks(event.good_jets, event.b_quarks_h, "hb", 2)

        #Number of reco jets matched to quarks from W, top, higgs
        event.nMatch_wq = 0
        event.nMatch_tb = 0
        event.nMatch_hb = 0
        #As above, but also required to be tagged/untagged for b/light respectively.
        event.nMatch_wq_btag = 0
        event.nMatch_tb_btag = 0
        event.nMatch_hb_btag = 0

        for ij, jet in enumerate(event.good_jets):

            jet.tth_match_label = None
            jet.tth_match_index = -1
            jet.tth_match_dr = -1
            jet.tth_match_label_numeric = -1

            if not matched_pairs.has_key(ij):
                continue
            mlabel, midx, mdr, mlabel_num = matched_pairs[ij]

            jet.tth_match_label = mlabel
            jet.tth_match_index = midx
            jet.tth_match_dr = mdr
            jet.tth_match_label_numeric = mlabel_num

            if mlabel == "wq":
                event.nMatch_wq += 1

                #If this jet is considered to be un-tagged (CSV or LR)
                if jet.btagFlag < 0.5:
                    event.nMatch_wq_btag += 1
            if mlabel == "tb":
                event.nMatch_tb += 1

                #If this jet is considered to be b-tagged (CSV or LR)
                if jet.btagFlag >= 0.5:
                    event.nMatch_tb_btag += 1

            if mlabel == "hb":
                event.nMatch_hb += 1
                if jet.btagFlag >= 0.5:
                    event.nMatch_hb_btag += 1

        if "matching" in self.conf.general["verbosity"]:
            matches = {"wq":event.l_quarks_w, "tb": event.b_quarks_t, "hb":event.b_quarks_h}

            for ij, jet in enumerate(event.good_jets):
                if not matched_pairs.has_key(ij):
                    continue
                mlabel, midx, mdr = matched_pairs[ij]
                print "jet match", ij, mlabel, midx, mdr, jet.pt, matches[mlabel][midx].pt

        #reco-level tth-matched system
        spx = 0.0
        spy = 0.0
        for jet in event.good_jets:
            if not (jet.tth_match_label is None):
                p4 = lvec(jet)
                spx += p4.Px()
                spy += p4.Py()
        for lep in event.good_leptons:
            p4 = lvec(lep)
            match = False
            for glep in event.lep_top:
                p4g = lvec(glep)
                if p4g.DeltaR(p4) < 0.3:
                    match = True
                    break
            if match:
                spx += p4.Px()
                spy += p4.Py()

        event.tth_px_reco = spx
        event.tth_py_reco = spy

        #Calculate tth recoil
        #rho = -met - tth_matched
        event.tth_rho_px_reco = -event.met[0].px - event.tth_px_reco
        event.tth_rho_py_reco = -event.met[0].py - event.tth_px_reco

        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class MEAnalyzer(FilterAnalyzer):
    """
    Performs ME calculation using the external integrator.
    It supports multiple MEM algorithms at the same time, configured via the
    self.configs dictionary. The outputs are stored in a vector in the event.

    The ME algorithms are run only in case the njet/nlep/Wtag category (event.cat)
    is in the accepted categories specified in the config.
    Additionally, we require the b-tagging category (event.cat_btag) to be "H" (high).

    For each ME configuration on each event, the jets which are counted to be b-tagged
    in event.btagged_jets are added as the candidates for t->b (W) or h->bb.
    These jets must be exactly 4, otherwise no permutation is accepted (in case
    using BTagged/QUntagged assumptions).

    Any additional jets are assumed to come from the hadronic W decay. These are
    specified in event.wquark_candidate_jets.

    Based on the event njet/nlep/Wtag category, if a jet fmor the W is counted as missing,
    it is integrated over using additional variables set by self.vars_to_integrate.

    The MEM top pair hypothesis (di-leptonic or single leptonic top pair) is chosen based
    on the reconstructed lepton multiplicity (event.good_leptons).

    The algorithm is shortly as follows:
    1. check if event passes event.cat and event.cat_btag
    2. loop over all MEM configurations i=[0...Nmem)
        2a. add all 4 b-tagged jets to integrator
        2b. add all 0-3 untagged jets to integrator
        2c. add all leptons to integrator
        2d. decide SL/DL top pair hypo based on leptons
        2e. based on event.cat, add additional integration vars
        2f. run ME integrator for both tth and ttbb hypos
        2g. save output in event.mem_output_tth[i] (or ttbb)
        2i. clean up event in integrator

    Relies on:
    event.good_jets, event.good_leptons, event.cat, event.input.met_pt

    Produces:
    mem_results_tth (MEMOutput): probability for the tt+H(bb) hypothesis
    mem_results_ttbb (MEMOutput): probability for the tt+bb hypothesis

    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(MEAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

        self.configs = {
            "default": MEM.MEMConfig(),
            "MissedWQ": MEM.MEMConfig(),
            "oldTF": MEM.MEMConfig(),
            "NumPointsDouble": MEM.MEMConfig(),
            "NumPointsHalf": MEM.MEMConfig(),
            "NoJacobian": MEM.MEMConfig(),
            "NoDecayAmpl": MEM.MEMConfig(),
            "NoPDF": MEM.MEMConfig(),
            "NoScattAmpl": MEM.MEMConfig(),
            "QuarkEnergy98": MEM.MEMConfig(),
            "QuarkEnergy10": MEM.MEMConfig(),
            "NuPhiRestriction": MEM.MEMConfig(),
            "JetsPtOrder": MEM.MEMConfig(),
            "JetsPtOrderIntegrationRange": MEM.MEMConfig(),
            "Recoil": MEM.MEMConfig(),
            "Sudakov": MEM.MEMConfig(),
            "Minimize": MEM.MEMConfig(),
        }
        for k in self.configs.keys():
            self.configs[k].transfer_function_method = MEM.TFMethod.External

        self.memkeys = self.conf.mem["methodsToRun"]

        self.configs["default"].defaultCfg()
        self.configs["oldTF"].defaultCfg()
        self.configs["MissedWQ"].defaultCfg()
        self.configs["NumPointsDouble"].defaultCfg(2.0)
        self.configs["NumPointsHalf"].defaultCfg(0.5)
        self.configs["NoJacobian"].defaultCfg()
        self.configs["NoDecayAmpl"].defaultCfg()
        self.configs["NoPDF"].defaultCfg()
        self.configs["NoScattAmpl"].defaultCfg()
        self.configs["QuarkEnergy98"].defaultCfg()
        self.configs["QuarkEnergy10"].defaultCfg()
        self.configs["NuPhiRestriction"].defaultCfg()
        self.configs["JetsPtOrder"].defaultCfg()
        self.configs["JetsPtOrderIntegrationRange"].defaultCfg()
        self.configs["Recoil"].defaultCfg()
        self.configs["Sudakov"].defaultCfg()
        self.configs["Minimize"].defaultCfg()

        self.configs["oldTF"].transfer_function_method = MEM.TFMethod.Builtin
        self.configs["NoJacobian"].int_code &= ~ MEM.IntegrandType.Jacobian
        self.configs["NoDecayAmpl"].int_code &= ~ MEM.IntegrandType.DecayAmpl
        self.configs["NoPDF"].int_code &= ~ MEM.IntegrandType.PDF
        self.configs["NoScattAmpl"].int_code &=  ~ MEM.IntegrandType.ScattAmpl
        self.configs["QuarkEnergy98"].j_range_CL = 0.98
        self.configs["QuarkEnergy98"].b_range_CL = 0.98
        self.configs["QuarkEnergy10"].j_range_CL = 0.10
        self.configs["QuarkEnergy10"].b_range_CL = 0.10
        self.configs["NuPhiRestriction"].m_range_CL = 99
        self.configs["JetsPtOrder"].highpt_first  = 0
        self.configs["JetsPtOrderIntegrationRange"].highpt_first  = 0
        self.configs["JetsPtOrderIntegrationRange"].j_range_CL = 0.99
        self.configs["JetsPtOrderIntegrationRange"].b_range_CL = 0.99
        self.configs["Recoil"].int_code |= MEM.IntegrandType.Recoil
        self.configs["Sudakov"].int_code |= MEM.IntegrandType.Sudakov
        self.configs["Minimize"].do_minimize = 1
        self.configs["Minimize"].int_code = 0

        for cfn, cfg in self.configs.items():
            cfg.mem_assumptions = set([])
            #A function Event -> boolean which returns true if this ME should be calculated
            cfg.do_calculate = lambda x: True
            cfg.enabled = True
        
        self.configs["default"].do_calculate = (
            lambda x: len(x.wquark_candidate_jets) >= 2
        )
        self.configs["MissedWQ"].mem_assumptions.add("missed_wq")

        #Can't integrate in dilepton
        self.configs["MissedWQ"].do_calculate = (
            lambda x: len(x.good_leptons) == 1 and
            len(x.wquark_candidate_jets) >= 1
        )
        #only in 6J SL
        self.configs["Sudakov"].do_calculate = (
            lambda x: len(x.good_jets) == 6 and
            len(x.good_leptons) == 1
        )

        #Create the ME integrator.
        #Arguments specify the verbosity
        self.integrator = MEM.Integrand(
            #0,
            MEM.output,
            self.configs["default"]
        )

        #Create an emtpy std::vector<MEM::Permutations::Permutations>
        self.permutations = CvectorPermutations()

        #Assume that only jets passing CSV>0.5 are b quarks
        self.permutations.push_back(MEM.Permutations.BTagged)

        #Assume that only jets passing CSV<0.5 are l quarks
        self.permutations.push_back(MEM.Permutations.QUntagged)

        self.integrator.set_permutation_strategy(self.permutations)

        #Pieces of ME to calculate
        # self.integrator.set_integrand(
        #     MEM.IntegrandType.Constant
        #     |MEM.IntegrandType.ScattAmpl
        #     |MEM.IntegrandType.DecayAmpl
        #     |MEM.IntegrandType.Jacobian
        #     |MEM.IntegrandType.PDF
        #     |MEM.IntegrandType.Transfer
        # )
        #self.integrator.set_sqrts(13000.);

        #Create an empty vector for the integration variables
        self.vars_to_integrate = CvectorPSVar()

    def add_obj(self, objtype, **kwargs):
        """
        Add an event object (jet, lepton, MET) to the ME integrator.

        objtype: specifies the object type
        kwargs: p4s: spherical 4-momentum (pt, eta, phi, M) as a tuple
                obsdict: dict of additional observables to pass to MEM
                tf_dict: Dictionary of MEM.TFType->TF1 of transfer functions
        """
        if kwargs.has_key("p4s"):
            pt, eta, phi, mass = kwargs.pop("p4s")
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(pt, eta, phi, mass);
        elif kwargs.has_key("p4c"):
            v = ROOT.TLorentzVector(*kwargs.pop("p4c"))
        obs_dict = kwargs.pop("obs_dict", {})
        tf_dict = kwargs.pop("tf_dict", {})

        o = MEM.Object(v, objtype)

        #Add observables from observable dictionary
        for k, v in obs_dict.items():
            o.addObs(k, v)
        for k, v in tf_dict.items():
            o.addTransferFunction(k, v)
        self.integrator.push_back_object(o)

    def beginLoop(self, setup):
        super(MEAnalyzer, self).beginLoop(setup)

    def configure_mem(self, event, mem_cfg):
        self.integrator.set_cfg(mem_cfg)
        self.vars_to_integrate.clear()
        self.integrator.next_event()
        mem_cfg.enabled = True

        missed_wq_cat = event.cat in ["cat2", "cat3"]
        can_integrate_wq = event.cat in ["cat1", "cat2", "cat3"]
        #One quark from W missed, integrate over its direction if possible
        if "missed_wq" in mem_cfg.mem_assumptions:
            if can_integrate_wq:
                self.vars_to_integrate.push_back(MEM.PSVar.cos_qbar1)
                self.vars_to_integrate.push_back(MEM.PSVar.phi_qbar1)
            else:
                mem_cfg.enabled = False

        Use_subjets = True

        if not Use_subjets:

            #Add heavy flavour jets that are assumed to come from top/higgs decay
            for jet in event.btagged_jets:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.pt, jet.eta, jet.phi, jet.mass),
                    obs_dict={MEM.Observable.BTAG: jet.btagFlag},
                    tf_dict={
                        MEM.TFType.bReco: jet.tf_b, MEM.TFType.qReco: jet.tf_l,
                        MEM.TFType.bLost: jet.tf_b, MEM.TFType.qLost: jet.tf_l
                    }
                )

            #Add light jets that are assumed to come from hadronic W decay
            for jet in event.wquark_candidate_jets:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.pt, jet.eta, jet.phi, jet.mass),
                    obs_dict={MEM.Observable.BTAG: jet.btagFlag},
                    tf_dict={
                        MEM.TFType.bReco: jet.tf_b, MEM.TFType.qReco: jet.tf_l,
                        MEM.TFType.bLost: jet.tf_b, MEM.TFType.qLost: jet.tf_l
                    }
                )

        else:
    
            #Add heavy flavour jets that are assumed to come from top/higgs decay
            for jet in event.btagged_jets_minus_sj:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.pt, jet.eta, jet.phi, jet.mass),
                    obs_dict={MEM.Observable.BTAG: jet.btagFlag},
                    tf_dict={
                        MEM.TFType.bReco: jet.tf_b, MEM.TFType.qReco: jet.tf_l,
                        MEM.TFType.bLost: jet.tf_b, MEM.TFType.qLost: jet.tf_l
                    }
                )

            #Add heavy flavour subjet that is assumed to come from top/higgs decay
            for jet in event.btagged_subjet:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.Pt(), jet.Eta(), jet.Phi(), jet.M()),
                    obs_dict={MEM.Observable.BTAG: jet.btagFlag},
                    tf_dict={
                        MEM.TFType.bReco: jet.tf_sj_b,
                        MEM.TFType.qReco: jet.tf_sj_l,
                        MEM.TFType.bLost: jet.tf_sj_b,
                        MEM.TFType.qLost: jet.tf_sj_l
                    }
                )

            #Add light jets that are assumed to come from hadronic W decay
            for jet in event.wquark_candidate_jets_minus_sj:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.pt, jet.eta, jet.phi, jet.mass),
                    obs_dict={MEM.Observable.BTAG: jet.btagFlag},
                    tf_dict={
                        MEM.TFType.bReco: jet.tf_b, MEM.TFType.qReco: jet.tf_l,
                        MEM.TFType.bLost: jet.tf_b, MEM.TFType.qLost: jet.tf_l
                    }
                )

            #Add light subjets that are assumed to come from hadronic W decay
            for jet in event.wquark_candidate_subjets:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.Pt(), jet.Eta(), jet.Phi(), jet.M()),
                    obs_dict={MEM.Observable.BTAG: jet.btagFlag},
                    tf_dict={
                        MEM.TFType.bReco: jet.tf_sj_b,
                        MEM.TFType.qReco: jet.tf_sj_l,
                        MEM.TFType.bLost: jet.tf_sj_b,
                        MEM.TFType.qLost: jet.tf_sj_l
                    }
                )



        for lep in event.good_leptons:
            self.add_obj(
                MEM.ObjectType.Lepton,
                p4s=(lep.pt, lep.eta, lep.phi, lep.mass),
                obs_dict={MEM.Observable.CHARGE: lep.charge},
            )
        self.add_obj(
            MEM.ObjectType.MET,
            #MET is caused by massless object
            p4s=(event.input.met_pt, 0, event.input.met_phi, 0),
        )

    #Check if event.nMatch_label >= conf.mem[cat][label]
    def require(self, required_match, label, event):
        nreq = required_match.get(label, None)
        if nreq is None:
            return True

        nmatched = getattr(event, "nMatch_"+label)

        #In case event did not contain b form higgs (e.g. ttbb)
        if "hb" in label and len(event.GenBQuarkFromH) < nreq:
            return True

        passes = (nmatched >= nreq)
        return passes

    def process(self, event):
        self.counters["processing"].inc("processed")

        #Clean up any old MEM state
        self.vars_to_integrate.clear()
        self.integrator.next_event()

        #Initialize members for tree filler
        event.mem_results_tth = []
        event.mem_results_ttbb = []

        #jets = sorted(event.good_jets, key=lambda x: x.pt, reverse=True)
        leptons = event.good_leptons
        met_pt = event.input.met_pt
        met_phi = event.input.met_phi

        if "reco" in self.conf.general["verbosity"]:
            for j in jets:
                print "jet", j.pt, j.eta, j.phi, j.mass, j.btagCSV, j.btagFlag, j.mcFlavour
            for l in leptons:
                print "lep", l.pt, l.eta, l.phi, l.mass, l.charge

        #Check if event passes reco-level requirements to calculate ME
        if event.cat_btag == "H":
            print "MEM RECO PASS", (event.input.run, event.input.lumi, event.input.evt,
                event.cat, event.btag_LR_4b_2b, len(event.btagged_jets),
                len(event.wquark_candidate_jets), len(event.good_leptons),
                len(event.btagged_jets), len(event.buntagged_jets)
            )
        else:
            print 'Not calculating Matrix Element'
            print 'event.cat_btag should be "H", but is "{0}"'.format(
                event.cat_btag )
            #Don't calculate ME
            return True


        #Here we optionally restrict the ME calculation to only matched events
        #Get the conf dict specifying which matches we require
        required_match = self.conf.mem.get("requireMatched", {}).get(event.cat, {})

        #Calculate all the match booleans
        #match label -> matched boolean
        passd = {
            p: self.require(required_match, p, event) for p in
            ["wq", "wq_btag", "tb", "tb_btag", "hb", "hb_btag"]
        }

        #Fail if any fails
        for k, v in passd.items():
            if not v:
                #print "Failed to match", k
                return True

        fstate = MEM.FinalState.TTH
        if len(event.good_leptons) == 2:
            fstate = MEM.FinalState.LL
        if len(event.good_leptons) == 1:
            fstate = MEM.FinalState.LH

        res = {}
        for hypo in [MEM.Hypothesis.TTH, MEM.Hypothesis.TTBB]:
            for confname in self.memkeys:
                mem_cfg = self.configs[confname]

                #Run MEM if we did not explicitly disable it
                if (self.conf.mem["calcME"] and
                        mem_cfg.do_calculate(event) and mem_cfg.enabled
                    ):
                    print "MEM started", ("hypo", hypo), ("conf", confname)
                    self.configure_mem(event, mem_cfg)
                    r = self.integrator.run(
                        fstate,
                        hypo,
                        self.vars_to_integrate
                    )
                    print "MEM done", ("hypo", hypo), ("conf", confname)

                    res[(hypo, confname)] = r
                else:
                    r = MEM.MEMOutput()
                    res[(hypo, confname)] = r

        if "default" in self.memkeys:
            p1 = res[(MEM.Hypothesis.TTH, "default")].p
            p2 = res[(MEM.Hypothesis.TTBB, "default")].p

            #In case of an erroneous calculation, print a message
            if self.conf.mem["calcME"] and (p1<=0 or p2<=0 or (p1 / (p1+0.02*p2))<0.0001):
                print "MEM BADPROB", p1, p2

        #print out full MEM result dictionary
        #print "RES", [(k, res[k].p) for k in sorted(res.keys())]

        event.mem_results_tth = [res[(MEM.Hypothesis.TTH, k)] for k in self.memkeys]
        event.mem_results_ttbb = [res[(MEM.Hypothesis.TTBB, k)] for k in self.memkeys]
