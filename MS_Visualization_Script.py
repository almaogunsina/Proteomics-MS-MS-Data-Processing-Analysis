#Import modules
import sys
from base64 import b64decode
from array import array
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt


if len(sys.argv) < 4:
    print("Please provide the mzxml file, scan number, and peptide sequence")
    sys.exit(1)

mzxml_file = sys.argv[1]
scan_number = sys.argv[2]
peptide_seq = sys.argv[3]
peptide_seq = peptide_seq.upper()


ns = '{http://sashimi.sourceforge.net/schema/}'

mzs = None
ints = None
try:
    #Access the peak information for the scan number
    for event,ele in ET.iterparse(mzxml_file):
        if ele.tag == ns + "scan" and ele.attrib.get("num") == scan_number:
            peaks_element = ele.find(ns + "peaks")
            if peaks_element is not None:
                peaks_data = peaks_element.text

                #peakselt is the XML element corresponding to the peaks list
                peaks = array('f',b64decode(peaks_element.text))
                if sys.byteorder != 'big':
                    peaks.byteswap()

                mzs = peaks[::2]
                ints = peaks[1::2]

            ele.clear()

    #Checks if scan number exists (if there is data in mzs and ints)
    if mzs is None or ints is None:
        print("Scan number was not found in the mzxml file.")
        sys.exit(1)

except IOError:
    print("The mzxml file was not found.")
    sys.exit(1)


#Dictionary of molecular weights
mw = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
      'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
      'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
      'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }


#Function to calculate the molecular weight for each b-ion/y-ion
def calculate_mz(ion):
    molecular_weight = 0.00
    for aa in ion: #Iterate through each amino acid in the ion
        try:
            molecular_weight += mw[aa] #Find the MW for the aa in dictionary, then add the MW to the molecular_weight
        except KeyError:
            print("Amino acid in the peptide sequence is not recognized.")
            sys.exit(1)
    return molecular_weight


#Function to match m/z values to peaks
def match_peaks(mzs, ints, b_ions, y_ions):
    matched_b_ions = []
    matched_y_ions = []
    unmatched_peaks = []

    #Intensity threshold
    max_intensity = max(ints)
    intensity_threshold = max_intensity * 0.05

    #Iterate through each peak by index...
    for i in range(len(mzs)):
        mz = mzs[i] #Extract the m/z
        intensity = ints[i] #Extract the intensity
        label = None #Initialize the label

        #If the intensity is below the threshold, it should not be matched
        if intensity <= intensity_threshold:
            unmatched_peaks.append({'mz': mz, 'intensity': intensity})
            continue

        #Check if b_ion m/z values match any of the peaks
        best_b_match = None
        smallest_b_diff = 0.05  #Initialize to tolerance
        for b_ion, mz_value in b_ions.items():
            diff = abs(mz - mz_value)
            if diff < smallest_b_diff:
                smallest_b_diff = diff
                best_b_match = b_ion

        #Add best b-ion match
        if best_b_match is not None:
            matched_b_ions.append({'label': best_b_match, 'mz': mz, 'intensity': intensity})
            label = best_b_match

        #Check if y_ion m/z values match any of the peaks
        if label is None:
            best_y_match = None
            smallest_y_diff = 0.05  #Initialize to tolerance
            for y_ion, mz_value in y_ions.items():
                diff = abs(mz - mz_value)
                if diff < smallest_y_diff:
                    smallest_y_diff = diff
                    best_y_match = y_ion

            #Add best y-ion match
            if best_y_match is not None:
                matched_y_ions.append({'label': best_y_match, 'mz': mz, 'intensity': intensity})
                label = best_y_match

        #If there is no matched ion to peak, add unmatched peak
        if label is None:
            unmatched_peaks.append({'mz': mz, 'intensity': intensity})

    return matched_b_ions, matched_y_ions, unmatched_peaks


#Function to plot spectrum
def plot_spectrum(unmatched_peaks, matched_b_ions, matched_y_ions):
    plt.figure(figsize=(12,6))

    #Plot unmatched peaks
    unmatched_mzs = [peak['mz'] for peak in unmatched_peaks]
    unmatched_ints = [peak['intensity'] for peak in unmatched_peaks]
    plt.stem(unmatched_mzs, unmatched_ints, linefmt='k-', markerfmt='', basefmt='k-', label='Unmatched Peaks')

    #Plot matched B-ions
    if matched_b_ions:
        b_mzs = [peak['mz'] for peak in matched_b_ions]
        b_ints = [peak['intensity'] for peak in matched_b_ions]
        plt.stem(b_mzs, b_ints, linefmt='b-', markerfmt='', basefmt='k-', label='Matched B-ions')

    #Plot matched Y-ions
    if matched_y_ions:
        y_mzs = [peak['mz'] for peak in matched_y_ions]
        y_ints = [peak['intensity'] for peak in matched_y_ions]
        plt.stem(y_mzs, y_ints, linefmt='r-', markerfmt='', basefmt='k-', label='Matched Y-ions')

    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.title(f'Mass Spectrum for Scan Number {scan_number}, Peptide {peptide_seq}')
    plt.legend()

    #Annotate matched peaks with their labels (only if label is not None)
    for peak in matched_b_ions:
        plt.annotate(peak['label'], (peak['mz'], peak['intensity']), textcoords="offset points", xytext=(0,10), ha='center')
    for peak in matched_y_ions:
        plt.annotate(peak['label'], (peak['mz'], peak['intensity']), textcoords="offset points", xytext=(0,10), ha='center')

    #Adjust y-axis
    plt.ylim(0, max(ints) * 1.1)

    #Show the plot
    plt.tight_layout()
    plt.show()


#Determine the length of the peptide sequence
peptide_length = len(peptide_seq)

#Create dictionary of the B-ions
b_ions = {}
for seq in range(1, peptide_length + 1):
    #Extract each b-ion from the N-terminus
    b_ion = peptide_seq[:seq]

    #Calculate the MW for each b-ion, and store the key-value pair in the dictionary
    b_ions["b"+ str(seq)] = calculate_mz(b_ion) + 1




#Create dictionary of the y-ions
y_ions = {}
for seq in range(peptide_length, 0, -1):
    #Extract each y-ion from the C-terminus
    y_ion = peptide_seq[seq-1:]

    #Calculate the MW for each y-ion, and store the key-value pair in the dictionary
    y_ions["y"+ str(peptide_length-seq+1)] = calculate_mz(y_ion) + 19



matched_b_ions, matched_y_ions, unmatched_peaks = match_peaks(mzs, ints, b_ions, y_ions)


#Print matched b-ions
if matched_b_ions:
    print("\nMatched B-ions:")
    for peak in matched_b_ions:
        print(peak['label'], peak['mz'], peak['intensity'])
else:
    print("\nNo matched b-ions.")


#Print matched y-ions
if matched_y_ions:
    print("\nMatched Y-ions:")
    for peak in matched_y_ions:
        print(peak['label'], peak['mz'], peak['intensity'])
else:
    print("\nNo matched y-ions.")


#Print unmatched peaks
if unmatched_peaks:
    print("\nUnmatched ions:")
    for peak in unmatched_peaks:
        print(peak['mz'], peak['intensity'])
else:
    print("\nNo unmatched peaks.")


plot_spectrum(unmatched_peaks, matched_b_ions, matched_y_ions)

#Five test scan numbers and peptide sequences
#python finalproject.py 17mix_test2.mzxml 1298 TYDSYLGDDYVR
#python finalproject.py 17mix_test2.mzxml 1301 TYDSYLGDDYVR
#python finalproject.py 17mix_test2.mzxml 1505 LITAIGDVVNHDPVVGDR
#python finalproject.py 17mix_test2.mzxml 845  AAQKPDVLTTGGGNPVGDK
#python finalproject.py 17mix_test2.mzxml 1506 LITAIGDVVNHDPVVGDR
