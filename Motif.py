##############################################################################################
##### ARG PARSE ####################################################
##############################################################################################

import argparse
def get_arguments():
    parser = argparse.ArgumentParser(description="Motif-Mark: This python script is to visualize motifs on sequences. \
    It's Python 3 compatible, inputs FASTA file and a motifs file, takes multiple sequences, multiple motifs, handles ambiguous motifs and outputs in an svg file.")
    parser.add_argument('-f', '--fasta_file', help='insert absolute file path, Assumes sequences are reveresed complimented', required=True, type=str)
    parser.add_argument('-m', '--motif_file', help='insert absolute file path', required=True, type=str)

    return parser.parse_args()


args = get_arguments()
fasta_file = args.fasta_file
motif_file = args.motif_file

##############################################################################################
##### FUNCTION FOR IUPCA AMBIGUOUS MOTIFS ####################################################
##############################################################################################

from Bio import Seq #biopython
from itertools import product #IUPAC 

def get_ambiguous_motifs(motif):
    all_motifs = []
    motif = motif.upper()
    #NOTE: U cant be read, so reaplace U with T
    motif =  motif.replace('U', 'T')
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    for i in product(*[ d[j] for j in motif ]):
        all_motifs.append(''.join(i))

    return(all_motifs)


##############################################################################################
################### FILES ####################################################################
##############################################################################################

#fasta_file = "/Users/daisychilin/shell/Bi610/winter_term/genomics/Figure_1.fasta"
#motif_file = '/Users/daisychilin/shell/Bi610/winter_term/genomics/Fig_1_motifs.txt'

##############################################################################################
################ MOTIF SCRIPT ################################################################
##############################################################################################

import re
import sys
import cairo
import math 

original_motif_list =[] 
ambiguous_motif_list = []
sequence_dict = {}
motif_colors_dict = {'Red':(255,0,0), 'Lime':(0,255,0), 'Blue':(0,0,205), 'Yellow':(255,255,0), 'Aqua':(0,255,255), \
               'Magenta':(255,0,255), 'Brown':(165,42,42), 'Green':(0,128,0), 'Purple':(128,0,128), 'Teal':(0,128,128)}
colors=[*motif_colors]

with open (fasta_file, "r") as fasta, \
    open(motif_file) as motif:
    #open motif txt file and find all ambigous motifs and put them into a list
    #find motifs in the sequence
    for line in motif:
        line = line.strip()
        original_motif_list.append(line)
        ambiguous_motifs = get_ambiguous_motifs(line)
        ambiguous_motif_list.append(ambiguous_motifs)
        
            #open fasta and put data into a dictionary.
    #header as key and sequence as value {header:sequence}
    for line in fasta:
        line = line.strip()
        if line.startswith(">"):
            header = line
            sequence_dict[header] = " "
        else:
            sequence_dict[header] += line

for key, value in sequence_dict.items():
    
    ########## DRAWING INTRONS (LINE) ##########
    width, height = (len(value)+200), 800
    surface = cairo.SVGSurface("motif.svg",width, height)
    context = cairo.Context(surface)
    #Need to tell cairo where to put the brush, the color and width, and the shape you want it to draw
    context.set_line_width(3)
    context.move_to(0,100)
    #lenght of line is the lenght of the sequence
    context.line_to(len(value),100)
    context.set_source_rgb(0,0,0)
    context.stroke()
    #legend
    context.rectangle(400,200,40,5)
    context.set_source_rgb(0, 0, 0)
    context.fill()
    context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(20)
    context.move_to(450, 205)
    context.set_source_rgb(0, 0, 0)
    context.show_text("Intron")
    
    ########## DRAWING EXONS (RECTANGLE) ##########
    #find exons in the sequence
    #NOTE: exons are capitalized 
    exon = re.findall(r'[A-Z]+', value)
    #get the lenght of exon, this will used for the length of the rectangle
    for i in exon:
        exon_length = len(i)
    #get the start position of this exon to inniate the coordinate rectangle
    for m in re.finditer(r'[A-Z]+', value):
        exon_start_coordinate = (m.start(0))
    #now draw a rectangle
    context.rectangle(exon_start_coordinate,50,exon_length,80)
    context.set_source_rgb(0,0,0)
    context.fill()
    context.stroke()
    # exon legend
    context.rectangle(415,225,10,40)
    context.set_source_rgb(0, 0, 0)
    context.fill()
    context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(20)
    context.move_to(445, 250)
    context.set_source_rgb(0, 0, 0)
    context.show_text("Exon")
    
    ########## DRAWING MOTIFS (TICK MARKS) ##########
    #find motifs in the sequence, grab the length and start position
    #your motifs are all capitalized now, so capitilize the seq now
    ctr = 0
    seq = value.upper()
    for ambiguous_motifs in ambiguous_motif_list:
        #color
        ctr_color=colors[ctr]
        motif_color=motif_colors_dict[ctr_color]
        red = motif_color[0]
        blue = motif_color[1]
        green = motif_color[2]
        context.set_source_rgb(red, green,blue)
        #legend
        context.rectangle(30, (200 + (ctr*50)), 5, 20)
        context.fill()
        context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(20)
        context.move_to(40, (215 + (ctr*50)))
        context.show_text(original_motif_list[ctr])
        
        ctr += 1
        
        for motifs in ambiguous_motifs:
        #Find motifs in the sequence 
            motif = re.findall(motifs, seq)
            #lenght
            for i in motif:
                motif_length = len(i)
            #start coordinate
            for m in re.finditer(motifs, seq):
                motif_start_coordinate = m.start(0)
                context.rectangle(motif_start_coordinate,50,motif_length,80)        
                context.fill()     
                
    ##label your output images    
    context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(20)
    context.move_to(30, 30)
    context.set_source_rgb(0, 0, 0)

    context.show_text(key)

    surface.write_to_png(key + ".png")
    
