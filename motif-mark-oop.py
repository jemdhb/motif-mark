#!/usr/bin/env python 
GLOBAL_MOTIF_OBJ=None
import cairo
import argparse as arg

#set line and exon color to be black
LINE_COLOR=[0, 0, 0]
EXON_COLOR=(0,0,0)
RGB=255

#motif colors
MOTIF_COLORS=[(17, 214, 129),#teal
      (98, 93, 194),#light purple
      (230, 163, 18),#brown/orange
      (224, 92, 209),#orchid
      (212, 104, 68),#salmon
      (104, 107, 103), #green 
      (217, 186, 13),#yellow
      (134, 27, 196), #purple         
      ]

#one per record
class Gene:
   #is our seq rc
   rc=False
   #name from header
   name=""
   #length of sequence
   gene_len=0
   #set of introns and exons
   gf_set=()
   #whole gene sequence
   seq=""

   def __init__(self, myseq_tup):
      """constructor for our gene class

      Args:
          myseq_tup (tuple): fasta record info
      """

      #determine name from header
      self.name =self.parse_header(myseq_tup[0])

      #get length of our sequence
      self.gene_len=len(myseq_tup[1])

      #get sequence content
      self.seq=myseq_tup[1]
      
      #split sequence into introns and exons
      self.gf_set =self.parse_seq()


   def parse_header(self,header_line):
      
      """Clean up our headers for pycairo visualization

      Args:
          header_line (str): _description_

      Returns:
          str: our clean gene name
      """

      #remove rc 
      if "reverse" in header_line:
         gene_name=header_line[1:header_line.index("reverse")-1]
     
      else:
         #get gene name no carrot
         gene_name=header_line[1:]
     
      return gene_name
   
   def parse_seq(self):
      """split our gene into introns and exons as gene_feat objects

      Returns:
          my_feats: a list of the introns and exons within the gene
      """

      isLower=True
      start_index=0
      my_feats=[]

      #for each char at position in the sequence
      for index, char in enumerate(self.seq):

         #case of prev index
         curr_case=isLower
         #determine case of curr index
         if ord(char)<96:
            isLower=False
         else:
            isLower=True

         #if our curr case does not match the prev
         if curr_case != isLower and start_index!=index:
            #we have a feature to build
            new_gene_feat=gene_feat(curr_case,start_index,index)

            my_feats.append(new_gene_feat)
            #shift to next feature 
            start_index=index

      #grab tail end exon
      my_feats.append(gene_feat(curr_case,start_index,self.gene_len))

      return my_feats

class gene_feat:
   #det how to draw
   is_intron=True
   #det where to draw
   start=0
   end=0
   def __init__(self, is_intron, start,end):
      """gene feature constructor

      Args:
          is_intron (bool): is this feature an intron?
          start (int): where is does this feature start in the gene
          end (int): where does this feature end in the gene
      """
      self.is_intron=is_intron
      self.start=start
      self.end=end


class motifs:
   #how many ways can our motifs be created
   wobble=[]
   #how long is our motif
   length=0
   #default sequence color
   color=(0,0,0)

   location=0

   def __init__(self, seq, wobble):
      self.wobble=wobble
      self.length=len(seq)
      self.seq=seq
   

      
# read in two lines at a time. make motifs at the start

#for drawing it depends on if its one drawing per gene or one drawing per file
def get_wobble(seq):
      """determine all motif possibilites accounting for wobble sequences

      Args:
          seq (str): our motif sequence of interest

      Returns:
         all_seqs: a set of motif strings
      """
      wobble={'y':['c','t'], "w":["a","t"], "s":["c","g"], "m":["a","c"],
              "k":["g","t"], "r":["a","g"], "b":["c","g","t"], "d":["a","g","t"],
              "h":["h","c","t"], "v":["a","c","g"], "n":["a","t","c","g"]}
      #so not building so many strings
      all_seqs=[[]]
      #for each char in seq
      for char in seq.lower().replace("u","t"):
         if char not in wobble:
            #just add char to all possibilities 
            all_seqs=[item+[char] for item in all_seqs]
         #if its a wobble sequence, have to determine all of the possibilities 
         else:
            #grab possibilities from dict
            poss=wobble[char]
            poss_index=0
            num_poss=len(poss)
            
            one_poss=len(all_seqs)
            #determine how many possibilities this wobble will produce
            all_seqs=all_seqs*num_poss
            
            #for all possibities
            for i in range(len(all_seqs)):
               #determine what wobble to append
               if i>=one_poss:
                  one_poss+=one_poss
                  poss_index+=1

               all_seqs[i]=all_seqs[i]+[poss[poss_index]]
               
      #turn our list of lists into a list of strings
      all_seqs=set(["".join(item) for item in all_seqs])
      return all_seqs

def det_all_wobbles(infile):
   """determine all of the possible wobble permutations

   Args:
       infile (str): path to our input file

   Returns:
       motif_wobble_list, valid_lens: list of all possible motif strings, 
       list of all possible motif lengths
   """
   # for one per file could get all the characters to be able to do the scale thing
   motif_wobble_list=open(infile,"r").read().split()
   #build motifs object
   valid_lens=set()
   for index,item in enumerate(motif_wobble_list):
      m=get_wobble(item)
      motif_wobble_list[index]=m
      valid_lens.add(len(item))

   return motif_wobble_list,valid_lens

def read_fasta_for_genes(infile):

   """run through fasta file, builing gene and gene_feat objects as we go

   Args:
       infile (str): path to our input fasta file

   Returns:
       record_of_interest, longest_gene: list of all gene objects, length of 
       longest gene (for drawing purposes)
   """

   records=open(infile,"r")
   line_num=1
   #store gene objects
   record_of_interest=[]
   #curr header
   header_line=""
   #curr seq
   seq_line=""

   longest_gene=-1

   for line in records:
      #build seq info
      if line[0]!=">":
         seq_line+=line.strip()
      #otherwise header
      else:
         #if first header skip
         if line_num==1:
            header_line=line.strip()
            continue
         #otherwise ready to build gene
         curr_gene=Gene((header_line,seq_line))

         curr_gene_len=curr_gene.gene_len

         if curr_gene_len > longest_gene:
            longest_gene=curr_gene_len

         record_of_interest.append(curr_gene)
         #update to new header
         header_line=line.strip()
         #reset 
         seq_line=""
      line_num+=1

   #grab last record
   record_of_interest.append(Gene((header_line,seq_line)))

   if record_of_interest[-1].gene_len > longest_gene:
            longest_gene=record_of_interest[-1].gene_len

   records.close()
   return record_of_interest, longest_gene


def create_motif_at_location(g,index, all_motifs):
   """create a motif object NOT drawing for every relevant index in our genes

   Args:
       g (Gene): our current gene/record of interest
       index (int): our current location in our gene of interest
       all_motifs (list): a list of all our current motif objects
   """
   #for each motif
   for jindex,wobble in enumerate(MOTIF_WOBBLE_LIST):
            #grab relevant color
            color=MOTIF_COLORS[jindex]
            #det motif length 
            moi_len=len(min(wobble))
            #grab current section of gene
            my_motif=g.seq[index:index+moi_len].lower()
            #if this sequence slice is one of our possible motifs 
            if my_motif in wobble and len(my_motif)==moi_len:
                  #build a motif and add it to our list
                  m=motifs(min(wobble), wobble)
                  m.location=index
                  m.color=color
                  all_motifs.append(m)

def draw_empty_panel(longest_feature, record_of_interest):
   """draw an empty white grid

   Args:
       longest_feature (int): length of our longest feature 

   Returns:
       _type_: _description_
   """
   x_dim=longest_feature+50
   legend_size=len(open(MOTIFFILE).read().split())*25
   y_dim=300*(len(record_of_interest))+legend_size
   y_panel_dim=-200
   #get x from function
   surface=cairo.ImageSurface (cairo.FORMAT_ARGB32, x_dim, y_dim)  
   # position for the text 
   context = cairo.Context(surface)
   context.set_source_rgba(1,1,1, 1)
   context.rectangle(0, 0, x_dim, y_dim) 
   context.fill()

   return context, surface, y_panel_dim

def draw_gene(context, y_panel_dim, g):
   """draw our gene object

   Args:
       context (pycairo): image we are writing too
       y_panel_dim (int): where we are in the image
       g (Gene): what gene we are on
   """
   #text colors
   font_size=35
   intron_height=5
   #distance from name to gene
   header_y_padding=100

   #header info
   context.set_source_rgba(LINE_COLOR[0]/RGB, LINE_COLOR[1]/RGB,
                           LINE_COLOR[2]/RGB, 1)
   context.move_to(10,25+y_panel_dim)  
   context.set_font_size(font_size) 
   context.select_font_face( 
      "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
   context.show_text(f"{g.name}") 
   context.stroke() 

   #gene info
   context.move_to(0,header_y_padding+y_panel_dim)  
   context.rectangle(0, header_y_padding+y_panel_dim, g.gene_len,
                     intron_height)
   context.fill()

    

def draw_motifs(all_motifs, context, y_panel_dim):
   """draw all of our motif objects onto our surface

   Args:
       all_motifs (list): list of all motif objects to be drawn
       context (pycairo): what to write our pics to
       y_panel_dim (int): where to place this motif on our panel (y axis wise)
   """
   motif_height=50
 
   for mot in all_motifs:
      #determine color from motif object
      context.set_source_rgba(mot.color[0]/RGB, mot.color[1]/RGB,
                              mot.color[2]/RGB, 0.50)
      #draw motif based on object and panel info
      context.rectangle(mot.location, 75+y_panel_dim, mot.length, motif_height) 
      context.fill() 
      context.stroke()

def draw_gene_features(g,context, y_panel_dim):
   """place our gene features (introns and exons) onto our pycairo image

   Args:
       g (Gene): Current gene we are on 
       context (pycairo): surface we are writing to
       y_panel_dim (int): where we are in the image (height-wise)
   """
   #center on gene line
   exon_y_centering=92
   exon_height=20
   #every gene feature
   for gf in g.gf_set:
      #already drawn at gene level
      if gf.is_intron:
         pass
      else:
         #draw exon square
         context.set_source_rgba(EXON_COLOR[0]/RGB, EXON_COLOR[1]/RGB,
                                 EXON_COLOR[2]/RGB, 0.85)
         context.rectangle(gf.start, exon_y_centering+y_panel_dim,
                           gf.end-gf.start, exon_height) 
         context.fill() 

def draw_legend(context, y_panel_dim):
   """Create a legend for our pycairo object

   Args:
       context (pycairo object): where to write our images
       y_panel_dim (int): where are we in the image (height-wise)
   """
   #legend sizing
   legend_key_size=15
   legend_key_padding=10
   legend_y_padding=35
   text_x_padding=30
   #have some padding from last gene drawn
   before_legend_position=y_panel_dim+300
   
   #Legend label
   context.move_to(legend_key_padding,before_legend_position-legend_key_padding*2)
   context.set_source_rgba(LINE_COLOR[0]/RGB, LINE_COLOR[1]/RGB, LINE_COLOR[2]/RGB, 1)
   context.show_text(f"Legend") 


   #for each possible motif
   for index,mot in enumerate(open("Fig_1_motifs.txt","r").read().split()):
         
         #determine assigned color
         context.set_source_rgba(MOTIF_COLORS[index][0]/RGB, MOTIF_COLORS[index][1]/RGB,
                                 MOTIF_COLORS[index][2]/RGB, 0.50)
         #draw legend key (little square)
         context.rectangle(legend_key_padding, before_legend_position,
                           legend_key_size,
                           legend_key_size) 

         #place motif to the right of the square
         context.move_to(text_x_padding, before_legend_position+legend_key_size)
         context.set_source_rgba(MOTIF_COLORS[index][0]/RGB, MOTIF_COLORS[index][1]/RGB,
                                 MOTIF_COLORS[index][2]/RGB, 0.75)
         context.show_text(f"{mot}") 
         #write this legend row
         context.fill() 
         context.stroke()
         #move down to next motif
         before_legend_position+=legend_y_padding


def draw_everything(record_of_interest, longest_feature):
   """draw where motifs occurr in all of our genes

   Args:
       record_of_interest (list): list of all of our relevant gene objects
       longest_feature (int): length of the longest gene (for drawing purposes)
   """
   all_motifs=[]
   #draw white panel to put our info on
   context, surface, y_panel_dim =draw_empty_panel(longest_feature, record_of_interest)
   
   #each gene object
   for g in record_of_interest:
      #store all motifs found
      all_motifs=[]
      #move lower in the panel
      y_panel_dim+=250
      #parse through every char searching for our motifs
      for index in range(g.gene_len):
          create_motif_at_location(g,index, all_motifs)

      #draw header and the span of the whole gene
      draw_gene(context, y_panel_dim, g)
      
      #move lower to draw our introns and exons for this gene
      context.move_to(0,100+y_panel_dim)  
      draw_gene_features(g,context, y_panel_dim)
      #draw our motifs on top of all this gene info
      draw_motifs(all_motifs, context, y_panel_dim)


   draw_legend(context, y_panel_dim)
   
   surface.write_to_png (f"{OUTFILE}") # Output to PNG
         

def get_args():
    parser = arg.ArgumentParser(description="")
    parser.add_argument("-f", "--file",
                     help="Path to input sam file.",
                     required=True, type=str)
    parser.add_argument("-m", "--motifile",
                     help="Path to output sam file with our PCR duplicates removed.",
                     required=True, type=str)
    return parser.parse_args()

args = get_args()
INFILE=args.file
MOTIFFILE=args.motifile

#create unique outfile name
OUTFILE=INFILE[:INFILE.rfind(".")]+".png"
#determine our wobble sequences
MOTIF_WOBBLE_LIST, VALID_LENS=det_all_wobbles(MOTIFFILE)
#create gene objects and determine the longest gene
RECORD_OF_INTEREST, LONGEST_GENE =read_fasta_for_genes(INFILE)
#draw everything :)))))))
draw_everything(RECORD_OF_INTEREST, LONGEST_GENE)

