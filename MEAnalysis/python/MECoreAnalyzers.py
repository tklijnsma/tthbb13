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

        self.R_cut = 0.3
        self.top_mass = 172.04

        self.bTagAlgo = self.conf.jets["btagAlgo"]

        if hasattr( self.conf, 'httCandidatecut' ):
            self.Cut_criteria = self.conf.httCandidatecut
            print 'Using httCandidate cut criteria from configuration file'
        else:
            self.Cut_criteria = [
                ( 'pt'  , '>', '200.0' ),
                ( 'mass', '>', '120.0' ),
                ( 'mass', '<', '220.0' ),
                ( 'fW'  , '<', '0.175' ) ]

    def beginLoop(self, setup):
        super(SubjetAnalyzer, self).beginLoop(setup)


    def endLoop(self, setup):

        print 'Running endLoop'


    def process(self, event):

        print 'Printing from SubjetAnalyzer! iEv = {0}'.format(event.iEv)

        ########################################
        # Minimal event suitability:
        #  - Number of quarks needs to be correct (2 b from top, 2 light)
        #  - Needs to be single leptonic
        #  - At least 1 httCandidate
        ########################################

        # Check if the event is single leptonic
        if not event.is_sl:
            return 0

        # Necessary to remove duplicates in GenWZQuark branch
        # This step can be removed for V12 samples
        self.CompareWZQuarks( event )

        # Check if event contains the right number of quarks
        genquarks_cat1_present = True
        genhiggs_present = True
        if len(event.GenWZQuark)<2:
            genquarks_cat1_present = False
        if len(event.GenWZQuark)>2:
            genquarks_cat1_present = False
        if len(event.GenBQuarkFromTop)<2:
            genquarks_cat1_present = False
        if len(event.GenBQuarkFromTop)>2:
            genquarks_cat1_present = False
        if len(event.GenBQuarkFromH)<2:
            genhiggs_present = False
        if len(event.GenBQuarkFromH)>2:
            genhiggs_present = False

        # Just run normal mem if there is no httCandidate present
        # Check if there is an httCandidate
        if len( event.httCandidate ) == 0:
            return True

        # Apply the cuts on the httCandidate
        tops = []
        for candidate in event.httCandidate:
            if self.Apply_Cut_criteria( candidate ):
                tops.append( copy.deepcopy(candidate) )

        # Just run normal mem if there is no httCandidate surviving the cuts
        # Check if any candidates survived the cutoff criteria
        if len(tops) == 0:
            return True

        # Calculates the delR with the (single) lepton that was found
        self.Set_DelRwithLepton( event, tops )

        # If exactly 1 survived, simply continue with that candidate
        if len(tops) == 1:
            top = tops[0]
            other_top_present = False

        # If more than 1 candidate survived the cutoff criteria, choose the
        # one whose delR with the lepton was biggest
        else:
            tops = sorted( tops, key=lambda x: -x.delR_lepton )
            other_top_present = True
            top = tops[0]
            other_top = tops[1]

        # Write delRmin to event
        setattr( event, 'httCandidate_delRmin', abs(top.Rmin-top.RminExpected) )
        print 'httCandidate_delRmin = {0}'.format( abs(top.Rmin-top.RminExpected) )

        # event.wquark_candidate_jets is a set instead of a list (not sure why)
        event.wquark_candidate_jets = list( event.wquark_candidate_jets )

        ########################################
        # Get the lists of particles: quarks, jets and subjets
        ########################################

        if genquarks_cat1_present:
            # Get the hadronic & leptonic b-quark and the two light quarks
            ( tl_HadronicBQuark,
              tl_GenWZQuark1,
              tl_GenWZQuark2,
              tl_LeptonicBQuark )  = self.Get_tl_genquarks( event )

            # Some convenient lists of quarks
            tl_QuarksFromHadrTop = [
                tl_HadronicBQuark, tl_GenWZQuark1, tl_GenWZQuark2 ]
            tl_LightQuarks = [ tl_GenWZQuark1, tl_GenWZQuark2 ]

        if genhiggs_present:
            # Get the 2 b-quarks from Higgs
            tl_BQuarksFromH = []
            for q in event.GenBQuarkFromH:
                tl_q = ROOT.TLorentzVector()
                tl_q.SetPtEtaPhiM( q.pt, q.eta, q.phi, q.mass )
                tl_BQuarksFromH.append( tl_q )

        # Get the subjets from the httCandidate
        # Also sets transfer functions as attributes, and a btagFlag based on kin.
        tl_subjets = self.Get_Subjets( top )
        if other_top_present:
            tl_subjets_other_top = self.Get_Subjets( other_top )

        # Get the list of btagged_jets
        tl_btagged_jets = []
        for jet in event.btagged_jets:
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM( jet.pt, jet.eta, jet.phi, jet.mass )
            setattr( x, 'origin_jet', jet )
            tl_btagged_jets.append( x )

        # Get list of wquark_candidate_jets
        tl_wquark_candidate_jets = []
        for jet in event.wquark_candidate_jets:
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM( jet.pt, jet.eta, jet.phi, jet.mass )
            setattr( x, 'origin_jet', jet )
            tl_wquark_candidate_jets.append( x )


        ########################################
        # Matching
        ########################################

        # Match subjet to a bjet
        Match_subjet_bjet = self.Match_two_tl_lists(
            tl_subjets, 'subjet',
            tl_btagged_jets, 'bjet' )

        # Match subjet to a ljet
        Match_subjet_ljet = self.Match_two_tl_lists(
            tl_subjets , 'subjet',
            tl_wquark_candidate_jets, 'ljet' )

        # In case of double matching, choose the match with lowest delR
        # (This is not expected to happen often)
        for tl_subjet in tl_subjets:
            if hasattr( tl_subjet, 'bjet_match' ) and \
                hasattr( tl_subjet, 'ljet_match' ) :
                print 'Double match detected'
                if tl_subjet.bjet_match_delR < tl_subjet.ljet_match_delR:
                    del tl_subjet.ljet_match
                    del tl_subjet.ljet_match_delR
                    Match_subjet_ljet -= 1
                else:
                    del tl_subjet.bjet_match
                    del tl_subjet.bjet_match_delR
                    Match_subjet_bjet -= 1


        if genquarks_cat1_present and genhiggs_present:

            # Quark matching
            # ==============

            # Quark matching is only done to compare events to truth level.
            # Events are not selected based on successful quark matching.

            # Matching quarks with subjets

            Match_hadr_bquark_subjet = self.Match_two_tl_lists(
                tl_HadronicBQuark, 'hadr_bquark',
                tl_subjets, 'subjet' )

            Match_lquark1_subjet = self.Match_two_tl_lists(
                tl_GenWZQuark1, 'lquark1',
                tl_subjets, 'subjet' )

            Match_lquark2_subjet = self.Match_two_tl_lists(
                tl_GenWZQuark2, 'lquark2',
                tl_subjets, 'subjet' )

            Match_lept_bquark_subjet = self.Match_two_tl_lists(
                tl_LeptonicBQuark, 'lept_bquark',
                tl_subjets, 'subjet' )

            Match_bquark_higgs1_subjet = self.Match_two_tl_lists(
                tl_BQuarksFromH[0], 'bquark_higgs1',
                tl_subjets, 'subjet' )

            Match_bquark_higgs2_subjet = self.Match_two_tl_lists(
                tl_BQuarksFromH[1], 'bquark_higgs2',
                tl_subjets, 'subjet' )

            # Matching quarks with jets

            Match_hadr_bquark_jet = self.Match_two_tl_lists(
                tl_HadronicBQuark, 'hadr_bquark',
                tl_btagged_jets, 'jet' )

            Match_lquark1_jet = self.Match_two_tl_lists(
                tl_GenWZQuark1, 'lquark1',
                tl_wquark_candidate_jets, 'jet' )

            Match_lquark2_jet = self.Match_two_tl_lists(
                tl_GenWZQuark2, 'lquark2',
                tl_wquark_candidate_jets, 'jet' )

            Match_lept_bquark_jet = self.Match_two_tl_lists(
                tl_LeptonicBQuark, 'lept_bquark',
                tl_btagged_jets, 'jet' )

            Match_bquark_higgs1_jet = self.Match_two_tl_lists(
                tl_BQuarksFromH[0], 'bquark_higgs1',
                tl_btagged_jets, 'jet' )

            Match_bquark_higgs2_jet = self.Match_two_tl_lists(
                tl_BQuarksFromH[1], 'bquark_higgs2',
                tl_btagged_jets, 'jet' )

            setattr( event, 'QMatching_sj_hadr_bquark'  , Match_hadr_bquark_subjet )
            setattr( event, 'QMatching_sj_lquark1'      , Match_lquark1_subjet )
            setattr( event, 'QMatching_sj_lquark2'      , Match_lquark2_subjet )
            setattr( event, 'QMatching_sj_lept_bquark'  , Match_lept_bquark_subjet )
            setattr( event, 'QMatching_sj_bquark_higgs1',
                Match_bquark_higgs1_subjet )
            setattr( event, 'QMatching_sj_bquark_higgs2',
                Match_bquark_higgs2_subjet )

            setattr( event, 'QMatching_jet_hadr_bquark'  , Match_hadr_bquark_jet )
            setattr( event, 'QMatching_jet_lquark1'      , Match_lquark1_jet )
            setattr( event, 'QMatching_jet_lquark2'      , Match_lquark2_jet )
            setattr( event, 'QMatching_jet_lept_bquark'  , Match_lept_bquark_jet )
            setattr( event, 'QMatching_jet_bquark_higgs1', Match_bquark_higgs1_jet )
            setattr( event, 'QMatching_jet_bquark_higgs2', Match_bquark_higgs2_jet )


        ########################################
        # Logic
        ########################################

        # Goal of this section is make sure exactly 1 subjet has btagFlag = 1.0,
        # and exactly 2 subjets have btagFlag = 0.0.

        # Regarding the strategy: a successful match is always prioritized.

        # b | l | Strategy
        # --+---+-------------------------------------------------------------------
        # 0 | 0 | Trust btagFlags based on prefix ('sjNonW'=b) (Strategy 1)
        # 0 | 3 |
        # --+---+-------------------------------------------------------------------
        # 0 | 1 | If the matched subjet was 'sjW1' or 'sjW2', trust the btagFlags
        #   |   | based on prefix. (Strategy 1)
        #   |   | Else, trust the matching, and tag one of the remaining subjets as
        #   |   | the b, based on highest pt (to be tested) (Strategy 2)
        # --+---+-------------------------------------------------------------------
        # 0 | 2 | Tag the remaining unmatched subjet as the b (Strategy 3)
        # --+---+-------------------------------------------------------------------
        # 1 | 0 | Trust the match with the b-jet, and set the other 2 subjets as the
        #   | 1 | l-jets (Strategy 4)
        #   | 2 | 
        # --+---+-------------------------------------------------------------------
        # 2 | 0 | Select the subjet that was matched to a b-jet with the highest
        #   | 1 | b-likelihood as the b, and set the remaining subjets as the l-jets
        # 3 | 0 | (Strategy 5)
        # --+---+-------------------------------------------------------------------

        # For every type of event, store the number of mismatches, the applied
        # strategy and the 'event_type_number' (see documentation, TODO):
        # ( nr_of_mismatches, strategy, event_type_number )

        #Do nothing; trust the btagFlags based on prefix
        if Match_subjet_bjet == 0 and Match_subjet_ljet == 0:
            ( nr_of_mismatches, strategy, event_type_number ) = ( 0, 1, 1 )
        elif Match_subjet_bjet == 0 and Match_subjet_ljet == 3:
            ( nr_of_mismatches, strategy, event_type_number ) = ( 1, 1, 2 )

        elif Match_subjet_bjet == 0 and Match_subjet_ljet == 1:

            # Get the subjet that was matched with a light jet
            tl_l_subjet = [ tl for tl in tl_subjets if hasattr(tl,'ljet_match') ][0]
            # Get the other two subjets
            tl_o_subjets = [ tl for tl in tl_subjets if not hasattr( tl,
                                                                     'ljet_match') ]
            ( nr_of_mismatches, strategy, event_type_number ) = ( 0, 1, 3 )

            if tl_l_subjet.prefix == 'sjNonW':
                # This means the 'sjNonW' subjet was matched with a light jet.
                # Trust the matching. From the remaining subjets, tag the one with
                # highest pt as the b
                ( nr_of_mismatches, strategy, event_type_number ) = ( 0, 2, 4 )
                setattr( tl_l_subjet, 'btagFlag', 0.0 )
                if tl_o_subjets[0].Pt() >  tl_o_subjets[1].Pt():
                    setattr( tl_o_subjets[0], 'btagFlag', 1.0 )
                    setattr( tl_o_subjets[1], 'btagFlag', 0.0 )
                else:
                    setattr( tl_o_subjets[0], 'btagFlag', 0.0 )
                    setattr( tl_o_subjets[1], 'btagFlag', 1.0 )
            # else: Trust btagFlags based on prefix

        elif Match_subjet_bjet == 0 and Match_subjet_ljet == 2:
            for tl in tl_subjets:
                if hasattr( tl, 'ljet_match' ): setattr( tl, 'btagFlag', 0.0 )
                else:                           setattr( tl, 'btagFlag', 1.0 )
            ( nr_of_mismatches, strategy, event_type_number ) = ( 0, 3, 5 )
            
        elif Match_subjet_bjet == 1 and Match_subjet_ljet in [0,1,2]:
            for tl in tl_subjets:
                if hasattr( tl, 'bjet_match' ): setattr( tl, 'btagFlag', 1.0 )
                else:                           setattr( tl, 'btagFlag', 0.0 )
            # Set appropiate ETN for easy accessing later
            if Match_subjet_ljet==0: ETN = 6
            if Match_subjet_ljet==1: ETN = 7
            if Match_subjet_ljet==2: ETN = 8
            ( nr_of_mismatches, strategy, event_type_number ) = ( 0, 4, ETN )

        elif Match_subjet_bjet in [2,3] and Match_subjet_ljet in [0,1]:

            tl_b_subjets = [ tl for tl in tl_subjets if hasattr(tl, 'bjet_match') ]
            tl_b_subjet = sorted(
                tl_b_subjets,
                key=lambda x: -getattr( x.bjet_match.origin_jet, self.bTagAlgo ))[0]

            for tl in tl_subjets:
                if tl == tl_b_subjet:           setattr( tl, 'btagFlag', 1.0 )
                else:                           setattr( tl, 'btagFlag', 0.0 )

            # Set appropiate ETN for easy accessing later
            if Match_subjet_bjet==2 and Match_subjet_ljet==0:
                NOM = 1
                ETN = 9
            if Match_subjet_bjet==2 and Match_subjet_ljet==1:
                NOM = 1
                ETN = 10
            if Match_subjet_bjet==3 and Match_subjet_ljet==0:
                NOM = 2
                ETN = 11
            ( nr_of_mismatches, strategy, event_type_number ) = ( NOM, 5, ETN )

        else:
            print 'Unassigned category! Create a strategy for this case.'
            print 'subjet-b matches: {0}, subjet-l matches = {1}'.format(
                Match_subjet_bjet, Match_subjet_ljet )
            return 0

        """
        # Check up printing - the httCandidate
        print '=====================================\nCheck print'
        for tl_subjet in tl_subjets:
            print '\nSubjet {0}'.format( tl_subjet.prefix )
            print '    btagFlag: {0}'.format( tl_subjet.btagFlag )
            if hasattr( tl_subjet, 'bjet_match' ):
                btag = getattr( tl_subjet.bjet_match.origin_jet, self.bTagAlgo ) 
                print '    original btag of the jet: {0}'.format(btag)
                self.Print_particle_lists(
                    ( [tl_subjet], 'TL', 'subjet'),
                    ( [tl_subjet.bjet_match], 'TL', 'Matching bjet')
                    )
            if hasattr( tl_subjet, 'ljet_match' ):

                self.Print_particle_lists(
                    ( [tl_subjet], 'TL', 'subjet'),
                    ( [tl_subjet.ljet_match], 'TL', 'Matching ljet')
                    )
        print '=====================================\n'

        # Check up printing - the input jets
        print '=====================================\n'
        print 'Input jets:'
        self.Print_particle_lists(
            ( event.btagged_jets, 'Class', 'event.btagged_jets'),
            ( event.wquark_candidate_jets, 'Class',
                'event.wquark_candidate_jets'),
            )
        print '=====================================\n'
        """

        # Remove all original instances and replace them with their matched subjets
        # If an original bjet is matched with a subjet that has btagFlag=0.0, the
        # original bjet is removed from the btagged_jets, and the subjet is appended
        # to event.wquark_candidate_jets.
        # (Analogous for original ljet matched with a subjet that has btagFlag=1.0)
        for tl_subjet in tl_subjets:
            if hasattr( tl_subjet, 'ljet_match' ):
                # Get the original ljet
                orig_ljet = tl_subjet.ljet_match.origin_jet
                # Remove it from the ljet list in the event
                event.wquark_candidate_jets.pop(
                    event.wquark_candidate_jets.index( orig_ljet ) )
            if hasattr( tl_subjet, 'bjet_match' ):
                # Get the original bjet
                orig_bjet = tl_subjet.bjet_match.origin_jet
                # Remove it from the bjet list in the event
                event.btagged_jets.pop( event.btagged_jets.index( orig_bjet ) )

            if tl_subjet.btagFlag == 1.0:
                # Append the subjet to btagged_jets list in the event
                event.btagged_jets.append( tl_subjet )
            elif tl_subjet.btagFlag == 0.0:
                # Append the subjet to wquark_candidate_jets list in the event
                event.wquark_candidate_jets.append( tl_subjet )

        """
        # Check up printing - the output jets
        print '=====================================\n'
        print 'Output jets:'
        self.Print_particle_lists(
            ( event.btagged_jets, 'Class', 'event.btagged_jets'),
            ( event.wquark_candidate_jets, 'Class',
                'event.wquark_candidate_jets'),
            )
        print '=====================================\n'
        """


        ########################################
        # Categorization
        ########################################

        if len( event.btagged_jets ) > 4:
            # The last added element is the subjet - keep that one, and the next
            # 3 btagged_jets
            event.btagged_jets = event.btagged_jets[-1:] + event.btagged_jets[:3]
        elif len( event.btagged_jets ) < 4:
            # Unpassable to MEM, needs precisely 4 btagged_jets
            event.cat_btag = 'NOCAT'

        if len( event.wquark_candidate_jets ) < 1:
            # Unpassable to MEM, can currently only deal with 1 missing light jet
            event.cat = 'NOCAT'

        if len( event.btagged_jets ) == 4:
            event.cat_btag = 'H'


        if len( event.btagged_jets ) == 4 and ( event.cat_btag == 'NOCAT' or event.cat == 'NOCAT' ):
            print 'Interesting event:'
            print 'This event should pass; nr of btagged_jets = 4, but cat_btag = NOCAT or cat = NOCAT'

            print '\nEvent suitability:'
            print 'n bjets  = {0}'.format( len( event.btagged_jets ) )
            print 'n ljets  = {0}'.format( len( event.wquark_candidate_jets ) )
            print 'cat_btag = {0}'.format( event.cat_btag )
            print 'cat      = {0}'.format( event.cat )


        ########################################
        # Write to event
        ########################################

        setattr( event, 'Matching_subjet_bjet', Match_subjet_bjet )
        setattr( event, 'Matching_subjet_ljet', Match_subjet_ljet )

        setattr( event, 'Matching_nr_of_mismatches', nr_of_mismatches )
        setattr( event, 'Matching_strategy', strategy )
        setattr( event, 'Matching_event_type_number', event_type_number )
        






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
    def Set_DelRwithLepton( self, event, tops ):

        # Get the lepton
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
    #--------------------------------------#
    

    def Get_Subjets( self, top ):

        tl_subjets = []
        prefixes = [ 'sjW1', 'sjW2', 'sjNonW' ]
        for (i_subjet, prefix) in enumerate( prefixes ):
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM(
                getattr( top, prefix + 'pt' ),
                getattr( top, prefix + 'eta' ),
                getattr( top, prefix + 'phi' ),
                getattr( top, prefix + 'mass' ) )
            setattr( x, 'prefix', prefix )

            # Set the btagFlag, based on prefix (only used if there is no other tag)
            if prefix == 'sjNonW':
                setattr( x, 'btagFlag', 1.0 )
            else:
                setattr( x, 'btagFlag', 0.0 )

            # Set pt, eta, phi and mass also as attributes
            # Needed for compatibility with mem
            setattr( x, 'pt', x.Pt() )
            setattr( x, 'eta', x.Eta() )
            setattr( x, 'phi', x.Phi() )
            setattr( x, 'mass', x.M() )

            tl_subjets.append( x )

        #Adding subjet transfer functions
        for subjet in tl_subjets:
            jet_eta_bin = 0
            if abs(subjet.Eta())>1.0:
                jet_eta_bin = 1

            #If True, TF [0] - reco, x - gen
            #If False, TF [0] - gen, x - reco
            eval_gen = False
            setattr( subjet, 'tf_b' ,
                self.conf.tf_sj_matrix['b'][jet_eta_bin].Make_Formula(eval_gen) )
            setattr( subjet, 'tf_l' ,
                self.conf.tf_sj_matrix['l'][jet_eta_bin].Make_Formula(eval_gen) )
            setattr( subjet, 'tf_b_lost' ,
                self.conf.tf_sj_matrix['b'][jet_eta_bin].Make_CDF() )
            setattr( subjet, 'tf_l_lost' ,
                self.conf.tf_sj_matrix['l'][jet_eta_bin].Make_CDF() )

            #Set jet pt threshold for CDF
            subjet.tf_b_lost.SetParameter(0, self.conf.jets["pt"])
            subjet.tf_l_lost.SetParameter(0, self.conf.jets["pt"])

        return tl_subjets



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
    #  - Output looks like:
    #    ( hadronic_BQuark, lightQuark1, lightQuark2, leptonic_BQuark )
    #    , where the all entries are TLorentzVector objects.
    def Get_tl_genquarks(self, event ):

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

        return ( tl_GenBQuarks[0],  # Hadronic
                 tl_GenWZQuarks[0], # Light 1
                 tl_GenWZQuarks[1], # Light 2
                 tl_GenBQuarks[1] ) # Leptonic
    #--------------------------------------#

    # Simple algorithm that matches the smallest delta R for two lists of TL vectors
    def Link_smallest_delR( self, tl_quarks, tl_jets ):

        n_jets = len(tl_jets)
        n_quarks = len(tl_quarks)

        Rmat = [[ (tl_quarks[i].DeltaR( tl_jets[j] )) \
            for j in range(n_jets)] for i in range(n_quarks) ]

        Rmin = 9999.0
        
        for (r, row) in enumerate(Rmat):
            for (c, ele) in enumerate(row):
                if ele < Rmin and ele < self.R_cut:
                    Rmin = ele
                    r_min = r
                    c_min = c

        if Rmin == 9999.0: return ( 'No link', 0, 0)

        return (r_min, c_min, Rmin)


    def Match_two_tl_lists( self, tls1_orig, label1, tls2_orig, label2 ):

        # If just a single TLorentzVector object was passed, convert it to a list
        if isinstance( tls1_orig, ROOT.TLorentzVector ):
            tls1_orig = [ tls1_orig ]
        if isinstance( tls2_orig, ROOT.TLorentzVector ):
            tls1_orig = [ tls2_orig ]

        # Create copies of the list, since entries will be popped
        tls1 = copy.deepcopy( tls1_orig )
        tls2 = copy.deepcopy( tls2_orig )

        # Save the original tl as an attribute
        for ( tl1, tl1_orig ) in zip( tls1, tls1_orig ):
            setattr( tl1, 'tl_origin', tl1_orig )
        for ( tl2, tl2_orig ) in zip( tls2, tls2_orig ):
            setattr( tl2, 'tl_origin', tl2_orig )

        # Attempt matching until the shortest list is depleted, or until there are
        # no more matches with delR < delR_cut
        n_matches = min( len(tls1), len(tls2) )

        for i_match in range(n_matches):

            # Attempt a match
            (i1, i2, delR) = self.Link_smallest_delR( tls1, tls2 )

            # Return the attempt number if no more matches could be made
            if i1 == 'No link':
                return i_match

            # Pop the matched tls from the lists for the next iteration
            matched_tl1 = tls1.pop(i1)
            matched_tl2 = tls2.pop(i2)

            # Record the match in the original tls
            setattr( matched_tl1.tl_origin,
                     '{0}_match'.format( label2 ),
                     matched_tl2.tl_origin )
            setattr( matched_tl2.tl_origin,
                     '{0}_match'.format( label1 ),
                     matched_tl1.tl_origin )

            # Record the delR value in the original tls
            setattr( matched_tl1.tl_origin,
                     '{0}_match_delR'.format( label2 ),
                     delR )
            setattr( matched_tl2.tl_origin,
                     '{0}_match_delR'.format( label1 ),
                     delR )

        return n_matches


#==========================END OF SUBJET ANALYZER==========================#

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

        Use_subjets = False

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

            print 'Number of btagged_jets:          {0}'.format(
                len( event.btagged_jets ) )
            print 'Number of wquark_candidate_jets: {0}'.format(
                len( event.wquark_candidate_jets ) )

            #print 'MEM currently disabled!'
            #return True
            

        else:
            #print 'Not calculating Matrix Element'
            #print 'event.cat_btag should be "H", but is "{0}"'.format(
            #    event.cat_btag )
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
