Proteomics MS/MS Spectrum Annotator

A Python command-line tool for visualising and annotating peptide fragmentation spectra from mzXML mass spectrometry files. Given a scan number and peptide sequence, the program computes theoretical b-ion and y-ion m/z values, matches them to observed peaks in the spectrum, and produces an annotated plot to assess whether the peptide is a good match to the spectrum.

Overview
Tandem mass spectrometry (MS/MS) generates fragmentation spectra that can be used to identify peptides. This tool automates the process of:

Parsing an mzXML file to extract peak data for a given scan
Computing theoretical m/z values for b-ions and y-ions from a peptide sequence
Matching theoretical values to observed peaks within a defined mass tolerance
Plotting the annotated spectrum so the user can visually evaluate the peptide-spectrum match


Requirements

Python 3.x
matplotlib

Install dependencies:
bashpip install matplotlib
Standard library modules used: sys, base64, array, xml.etree.ElementTree

Usage
bashpython annotator.py <mzXML_file> <scan_number> <peptide_sequence>

Arguments
ArgumentDescriptionmzXML_filePath to the mzXML input filescan_numberScan number to extract from the mzXML filepeptide_sequencePeptide sequence in single-letter amino acid code (e.g. PEPTIDE)

Example
bashpython MS_Visualization_Script.py sample.mzXML 1042 PEPTIDE

How It Works
1. Parsing the mzXML file
The program uses xml.etree.ElementTree to parse the mzXML file and locate the scan matching the provided scan number. Peak data (m/z values and intensities) are extracted by decoding the base64-encoded binary data using base64.b64decode and stored using array.

2. Computing theoretical ion m/z values — calculate_mz()
Given the peptide sequence, the program computes theoretical m/z values for:

b-ions — N-terminal fragments
y-ions — C-terminal fragments

3. Matching peaks — match_peaks()
Observed peaks are matched to theoretical ion m/z values using an absolute mass tolerance:
abs(observed m/z − theoretical m/z) < 0.05 Da
An intensity threshold of 5% of the most intense peak is applied before matching to reduce noise and focus annotation on significant peaks.

4. Plotting the spectrum — plot_spectrum()
The program generates an annotated spectrum plot using matplotlib, labelling matched peaks with their corresponding b-ion or y-ion assignments. The output plot allows the user to visually assess the quality of the peptide-spectrum match.

Parameters & Design Notes
ParameterValueNotesMass tolerance0.05 Da (absolute)May need adjustment for higher-resolution datasetsIntensity threshold5% of max intensityReduces noise; may exclude weak but meaningful peaks

Project Context
This project was completed as part of the M.S. in Bioinformatics programme at Georgetown University. It was developed to demonstrate applied skills in proteomics data processing, Python scripting, and scientific visualisation.

Author
Alma Ogunsina
M.S. Bioinformatics, Georgetown University
linkedin.com/in/almaogunsina
