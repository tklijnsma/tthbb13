#!/usr/bin/env python
"""
Helpers to calculate mutual information.

Inspired by 1408.3122
"""
 
########################################
# Imports and setup ROOT with style
########################################

import math
import glob
import copy
import sys

import ROOT

from TTH.Plotting.python.Helpers.PrepareRootStyle import myStyle

ROOT.gStyle.SetPadLeftMargin(0.2)
ROOT.gStyle.SetPadRightMargin(0.1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.25)

ROOT.gROOT.SetStyle("myStyle")
ROOT.gROOT.ForceStyle()

ROOT.gStyle.SetPadLeftMargin(0.26)
ROOT.gStyle.SetPadRightMargin(0.1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.29)

ROOT.gROOT.ForceStyle()


# Black-Body-Spectrum
# Reasons why it is nicer:
# http://root.cern.ch/drupal/content/rainbow-color-map
ROOT.gStyle.SetPalette(53)
# make the colors smoother
ROOT.gStyle.SetNumberContours(100)

ROOT.gStyle.SetPaintTextFormat("3.2f");

ROOT.TH1.SetDefaultSumw2()

ROOT.gErrorIgnoreLevel = 1

# Global variable to count drawn histograms for unique naming
i_draw = 0 

########################################
# Import private support code
########################################

# initializer: simple creation of bag-of-object classes
from Initializer import initializer

import TTH.Plotting.python.Helpers.OutputDirectoryHelper as OutputDirectoryHelper


########################################
# Define and crate output directory
########################################

output_dir = "OutputMutualInformation/"

# Create directory and subdirectories for all filetypes
OutputDirectoryHelper.CreateOutputDirs( output_dir )


########################################
# Helper Classes
########################################

class mi():
   """ Helper class to store all information for a set of mutual information plots."""
   @initializer
   def __init__(self, 
                name,                
                sample_name_signal,
                sample_name_background,
                li_vars,
                fiducial_cut = "(1)"):
      """ Constructor. Arguments:
      name                    : (string) name of the mutual information set
      sample_name_signal      : (string) signal sample to process
      sample_name_background  : (string) background sample to process
      vars                    : list of variable (from VariableHelpers) objects 
      fiducial_cut            : (string) fiducial cut (numerator and denominator)
      """
      pass
# end of class mi


########################################
# CalculateEntropy
########################################

def CalculateEntropy(h, verbose = False):
   """ Calculate and return the Shannon Entropy (base 2) of a histogram."""

   total = h.Integral()
   checksum = 0
   entropy = 0.

   # 1d Histograms
   if type(h) in [ROOT.TH1F, ROOT.TH1D]:
      for i_bin in range(1, h.GetXaxis().GetNbins()+1):

         p = 1. * h.GetBinContent(i_bin)/total
         if p>0:
            checksum += p
            entropy -= p*math.log(p,2) 

   # 2d Histograms
   elif type(h) in [ROOT.TH2F, ROOT.TH2D]:
      
      for i_bin_x in range(1, h.GetXaxis().GetNbins()+1):
         for i_bin_y in range(1, h.GetYaxis().GetNbins()+1):
            
            p = 1. * h.GetBinContent(i_bin_x, i_bin_y)/total
            if p>0:
               checksum += p
               entropy -= p*math.log(p,2) 

   # 3d Histograms
   elif type(h) in [ROOT.TH3F, ROOT.TH3D]:
      
      for i_bin_x in range(1, h.GetXaxis().GetNbins()+1):
         for i_bin_y in range(1, h.GetYaxis().GetNbins()+1):
            for i_bin_z in range(1, h.GetZaxis().GetNbins()+1):
            
               p = 1. * h.GetBinContent(i_bin_x, i_bin_y, i_bin_z)/total
               if p>0:
                  checksum += p
                  entropy -= p*math.log(p,2) 

   # Invalid Histogram type
   else:
      print "Ivalid type {0} for histogram in CalculateEntropy".format(type(h))
      print "Exiting"
   # End of handling different histogram dimensions

   if verbose:
      print "checksum = ", checksum

   return entropy

# End of CalculateEntropy      


########################################
# MakeHistogram
########################################

def MakeHistogram(tree, variables, n_bins, cuts):
   """ Produce a histgrom for given variables, number of bins per axis
   and cuts from a ROOT tree"""

   global i_draw

   tmp_name = "htmp{0}".format(i_draw)
   i_draw += 1

   # Make a TH1
   if len(variables)==1:
      draw_string = "{0}>>{1}({2},{3},{4})".format(variables[0].name, 
                                                   tmp_name,
                                                   n_bins,
                                                   variables[0].range_min,
                                                   variables[0].range_max)

   # Make a TH2
   elif len(variables)==2:
      draw_string = "{0}:{1}>>{2}({3},{4},{5},{6},{7},{8})".format(variables[1].name,  # y
                                                                   variables[0].name,  # x
                                                                   tmp_name,
                                                                   n_bins,
                                                                   variables[0].range_min,
                                                                   variables[0].range_max,
                                                                   n_bins,
                                                                   variables[1].range_min,
                                                                   variables[1].range_max)

   # Make a TH3
   elif len(variables)==3:
      draw_string = "{0}:{1}:{2}>>{3}({4},{5},{6},{7},{8},{9},{10},{11},{12})".format(variables[2].name,  # z
                                                                                      variables[1].name,  # y
                                                                                      variables[0].name,  # x
                                                                                      tmp_name,
                                                                                      n_bins,
                                                                                      variables[0].range_min,
                                                                                      variables[0].range_max,
                                                                                      n_bins,
                                                                                      variables[1].range_min,
                                                                                      variables[1].range_max,
                                                                                      n_bins,
                                                                                      variables[2].range_min,
                                                                                      variables[2].range_max)

   # Otherwise fail
   else:
      print "Invalid number of variables: ", len(variables)
      print "Exiting.."
      sys.exit()

   tree.Draw(draw_string, cuts)
   h = getattr(ROOT, tmp_name).Clone()
   h.SetDirectory(0)

   return h
# End of MakeHistogram


########################################
# Prepare Canvas
########################################

c = ROOT.TCanvas("","",800,800)


########################################
# Make the Plots
########################################

def MakePlots(mis, files, input_treename = 'tree'):

   # Some configuration
   f = 0.5 # fraction of signal events in combined samplke
           # Normalize samples so this value come out
           # Has the advantage that H(T)=1
           
   # Binning
   n_bins_one_var_one_sample = 400
   n_bins_one_var_two_sample = 400
   n_bins_two_var_one_sample = 50  # this is per axis, total bins -> n^2
   n_bins_two_var_two_sample = 50  # this is per axis, total bins -> n^2
   
   # Loop over the list of plots
   li_h2d = []

   for mi in mis:

      # open files
      infile_sig = ROOT.TFile(files[mi.sample_name_signal])
      infile_bkg = ROOT.TFile(files[mi.sample_name_background])

      # Get the Trees
      input_tree_sig = infile_sig.Get(input_treename)
      input_tree_bkg = infile_bkg.Get(input_treename)

      # Count events
      passed_fiducial_sig = input_tree_sig.Draw("(1)", mi.fiducial_cut)
      passed_fiducial_bkg = input_tree_bkg.Draw("(1)", mi.fiducial_cut)

      # Calculate overall additional background weight so that:
      # signal / signal+background = f
      # (set extra weight for signal to 1 as the overall normalization does not matter!)
      extra_weight_sig = 1.
      extra_weight_bkg = (passed_fiducial_sig - f * passed_fiducial_sig)/(f * passed_fiducial_bkg)

      mi_result  = {}

      # Loop over pairs of variables:
      for ivar1, var1 in enumerate(mi.li_vars):
         
         mi_result[var1.name] = {}

         for ivar2, var2 in enumerate(mi.li_vars):

            # Diagonal
            if (ivar1 == ivar2):

               # For the diagonal we want to calculate:
               # I(T;A) = H(sig,bkg)[A] - f * H(sig)[A] - (1-f) * H(bkg)[A]

               # Total cut is
               # fiducial + range(var1) + extra cut(var1)
               cut = "("
               cut += "({0})".format(mi.fiducial_cut)
               cut += "&&({0}>={1})".format(var1.name, var1.range_min)
               cut += "&&({0}<={1})".format(var1.name, var1.range_max)
               cut += "&&({0})".format(var1.extra_cut)
               cut += ")"

               cut_and_weight_sig = "{0}*{1}".format(cut, extra_weight_sig)
               cut_and_weight_bkg = "{0}*{1}".format(cut, extra_weight_bkg)

               h_sig = MakeHistogram(input_tree_sig, [var1], n_bins_one_var_one_sample, cut_and_weight_sig)
               h_bkg = MakeHistogram(input_tree_bkg, [var1], n_bins_one_var_one_sample, cut_and_weight_bkg)
               h_sigbkg_sig = MakeHistogram(input_tree_sig, [var1], n_bins_one_var_two_sample, cut_and_weight_sig)
               h_sigbkg_bkg = MakeHistogram(input_tree_bkg, [var1], n_bins_one_var_two_sample, cut_and_weight_bkg)

            # Below Diagonal
            elif ivar2 < ivar1:

               # For the off-diagonal we want to calculate:
               # I(T;A,B) = H(sig,bkg)[A,B] - f * H(sig)[A,B] - (1-f) * H(bkg)[A,B]

               # Total cut is
               # fiducial + range(var1) + extra cut(var1) + range(var2) + extra_cut(var2)
               cut = "("
               cut += "({0})".format(mi.fiducial_cut)
               cut += "&&({0}>={1})".format(var1.name, var1.range_min)
               cut += "&&({0}<={1})".format(var1.name, var1.range_max)
               cut += "&&({0})".format(var1.extra_cut)
               cut += "&&({0}>={1})".format(var2.name, var2.range_min)
               cut += "&&({0}<={1})".format(var2.name, var2.range_max)
               cut += "&&({0})".format(var2.extra_cut)
               cut += ")"

               cut_and_weight_sig = "{0}*{1}".format(cut, extra_weight_sig)
               cut_and_weight_bkg = "{0}*{1}".format(cut, extra_weight_bkg)

               h_sig = MakeHistogram(input_tree_sig, [var1, var2], n_bins_two_var_one_sample, cut_and_weight_sig)
               h_bkg = MakeHistogram(input_tree_bkg, [var1, var2], n_bins_two_var_one_sample, cut_and_weight_bkg)
               h_sigbkg_sig = MakeHistogram(input_tree_sig, [var1, var2], n_bins_two_var_two_sample, cut_and_weight_sig)
               h_sigbkg_bkg = MakeHistogram(input_tree_bkg, [var1, var2], n_bins_two_var_two_sample, cut_and_weight_bkg)

            # Above Diagonal
            else:
               continue
            # End of Diagonal vs Off-Diagonal difference
            
            # Combine signal/background
            h_sigbkg = h_sigbkg_sig
            h_sigbkg.Add(h_sigbkg_bkg)

            # Calculate entropies
            entropy_sig = CalculateEntropy(h_sig)
            entropy_bkg = CalculateEntropy(h_bkg)
            entropy_sigbkg = CalculateEntropy(h_sigbkg)

            # Calculate Mutual Information
            I = entropy_sigbkg  - f * entropy_sig  - (1-f) * entropy_bkg
            
            print "H_sig = {0} \t H_bkg = {1} \t H_sigbkg = {2} \t I = {3}".format(entropy_sig,
                                                                                   entropy_bkg,
                                                                                   entropy_sigbkg,
                                                                                   I)

            mi_result[var1.name][var2.name] = I
         # End var2 loop
      # End var1 loop

      ROOT.gStyle.SetPalette(1)   
      
      h_mi = ROOT.TH2F("","", 
                       len(mi.li_vars), -0.5, -0.5 + len(mi.li_vars),
                       len(mi.li_vars), -0.5, -0.5 + len(mi.li_vars))      

      for ivar1, var1 in enumerate(mi.li_vars):         

         h_mi.GetXaxis().SetBinLabel(ivar1+1, var1.pretty_name)
         h_mi.GetYaxis().SetBinLabel(ivar1+1, var1.pretty_name)

         for ivar2, var2 in enumerate(mi.li_vars):

            if ivar2 > ivar1:
               continue

            h_mi.SetBinContent(ivar1+1, ivar2+1, mi_result[var1.name][var2.name])


      h_mi.LabelsOption("v","X")
      h_mi.GetXaxis().SetLabelSize(0.035)
      h_mi.GetYaxis().SetLabelSize(0.035)
      h_mi.GetZaxis().SetLabelSize(0.03)
      
      draw_opts = "COLZ TEXT"

      h_mi.Draw(draw_opts)
      OutputDirectoryHelper.ManyPrint(c, output_dir, "{0}_mi".format(mi.name))

   # End mi loop
# End MakePlots
            
            

