#!/usr/bin/python

# See below for name and description
# Copyright (C) 2020 Richard J. Edwards <dr.r.edwards@icloud.com>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       synbad
Description:  Synteny-based scaffolding assessment and adjustment
Version:      0.8.0
Last Edit:    18/04/21
GitHub:       https://github.com/slimsuite/synbad
Copyright (C) 2020  Richard J. Edwards - See source code for GNU License Notice

Function:
    SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
    between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.

    Synbad will use or create:

    1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
    spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
    correspond to misassemblies.

    2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.

    Pairwise hits between the genomes are filtered according to the `minlocid=PERC` and `minloclen=INT` criteria, which
    by default limits hits to be at least 1kb and 50% identity. Note that this is applied after GABLAM has run, so it
    should be possible to re-run with more relaxed constraints and re-use the GABLAM tables.

    Next, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
    sorted and `QryGap` and `SbjGap` fields added. Any local alignments with flanking hits are then flagged in these
    new fields with 5', 3' or Both.

    The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
    corresponding local hits on the Qry and Sbj genomes. If there is also an inversion, `SynSpan` will be negative.
    If the local hits are against two different sequences from the other genome, the two sequence names will be
    entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
    small gaps), `SynSpan` or `GapSpan` will have a value of zero.

    Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values, and subsequent local hit
    analysis to fix/update ratings:

    * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
    * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
    * `Div` = `Divergence` = Translocation or Breakpoint gap with flanking collinear hits.
    * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
    * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.
    * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
    * `Inv` = `Inversion` = Flanking hits are on alternative strands.
    * `InvBrk` = `Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.
    * `InvDupFix` = `Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)
    * `InvFix` = `Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)
    * `Long` = `Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).
    * `Null` = No mapping between genomes for that gap.
    * `Span` = `Spanned` = Any gaps without Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
    * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).
    * `Tran` = `Translocation` = `SynSpan` indicates matches are on different scaffolds.
    * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.

Dependencies:
    SynBad needs Minimap2 installed. For `gapass` gap mode, Flye also needs to be installed. To generate
    documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    For full documentation of the SynBad workflow, run with `dochtml=T` and read the `*.docs.html` file generated.

Commandline:
    ### ~ Main SynBad run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
    genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
    basefile=X      : Prefix for output files [synbad]
    gablam=X        : Optional prefix for GABLAM search [defaults to $BASEFILE.map]
    gapmode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
    maxsynskip=INT  : Maximum number of local alignments to skip for SynTrans classification [4]
    maxsynspan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [25000]
    synreadspan=INT : Minimum number of reads spanning a gap to change the rating to "Spanned" [5]
    spannedflank=INT: Required flanking distance for synreadspan "Spanned" rating [0]
    maxoverlap=INT  : Maximum overlap (bp) of adjacent local hits to allow compression [500]
    chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
    chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
    ### ~ Fragmentation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
    fragtypes=LIST  : List of SynBad ratings to trigger fragmentation [Brk,Inv,InvBrk,Frag,Tran]
    minreadspan=INT : Min number of Span0 reads in gaps table to prevent fragmentation [1]
    ### ~ Additional input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    bam1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
    paf1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    bam2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
    paf2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]
    reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    mapflanks1=FILE : Flanks fasta file from previous SynBad run for mapping genome 1 flank identifiers []
    mapflanks2=FILE : Flanks fasta file from previous SynBad run for mapping genome 2 flank identifiers []
    fullmap=T/F     : Whether to abort if not all flanks can be mapped [True]
    ### ~ HiC Gap Flank options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hicbam1=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
    hicbam2=FILE    : Optional BAM file of HiC reads mapped onto assembly 2 [$BASEFILE2.HiC.bam]
    gapflanks=INT   : Size of gap flank regions to output for HiC pairing analysis (0=off) [10000]
    pureflanks=T/F  : Whether to restrict gap flanks to pure contig sequence (True) or include good gaps (False) [True]
    hicscore=X      : HiC scoring mode (score/wtscore) [wtscore]
    hicmode=X       : Pairwise HiC assessment scoring strategy (synbad/pure/rand/full) [synbad]
    hicdir1=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.qryflanks/]
    hicdir2=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.hitflanks/]
    ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    newacc1=X       : Scaffold name prefix for updated Genome 1 output [None]
    newacc2=X       : Scaffold name prefix for updated Genome 2 output [None]
    bestpair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
    update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]
    force=T/F       : Whether to force regeneration of SynBad results tables [False]
    dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_db, rje_lrbridge, rje_rmd, rje_seqlist, rje_sequence
import diploidocus, gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial working version without fragment=T implementation.
    # 0.1.1 - Added minlocid=PERC and minloclen=INT.
    # 0.1.2 - Added additional translocation skipping for SynTrans rating.
    # 0.1.3 - Modified the SynBad classification text.
    # 0.1.4 - Modified code to be able to run without long read mapping. Added dochtml output.
    # 0.2.0 - Added fragment=T output.
    # 0.3.0 - Added chromosome scaffold Translocation restriction.
    # 0.4.0 - Added an Duplication rating in place of Breakpoint for overlapping flanking hits; added top sequence pairs.
    # 0.5.0 - Added HiFi read type. Changed default to gapmode=gapspan.
    # 0.5.1 - Added some extra bug-fixes for running Diploidocus checkpos.
    # 0.6.0 - Added additional gap-spanning classes and qry-hit pair compression. Added dev reworking on main compression.
    # 0.6.1 - Tidy up of code and transition of reworked code to main run (no longer dev=T only).
    # 0.7.0 - Added code for fixing Inversions and well-supported translocations. Updated defaults.
    # 0.8.0 - Added hicbam=FILE inputs and initial HiC assessment. Fixed inversion bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add chr1 and chr2 chromosome prefix identifers (PAFScaff prefixes) to restrict Translocations to chromosomes.
    # [ ] : Add GapGFF function from rje_genomics with additional Gap rating information. (Or BED?)
    # [X] : Add fullcollapse=T/F option to collapse entire (pairwise) matching regions between gaps.
    # [ ] : Add tidying of gap ratings and suggested modifications to scaffolds.
    # [ ] : Make sure that Frag gap type is being handled appropriately.
    # [ ] : Add runmode with different options.
    # [Y] : Add hicbam file loading and flank/join assessements.
    # [ ] : Add modified fasta output following edits.
    # [ ] : Build LocusFixer chassis and use for re-assembling duplication gap regions.
    # [Y] : Execute inversions and output updated fasta file.
    # [ ] : Empirically establish the best values for maxoverlap, maxsynskip and maxsynspan.
    # [ ] : Option to add DepthCharge gaps where there might be misassemblies.
    # [ ] : Option to add KAT kmer analysis of regions to help assess duplicity. (Can use with depth or alone.)
    # [ ] : Standardise the fields (CamelCase)
    # [?] : Switch early from Qry/Hit to SeqName/Hit and change output file names? (Or just do this at end.)
    # [ ] : Add final output where Qry and Hit are replaced for Genome1/2 output.
    # [ ] : Rationalise the final output to the really useful stuff.
    # [?] : Add minimum HiC read count and/or HiC score filter?
    # [ ] : Add hicbest table with the top pairs and mark as single or reciprocal.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SynBad', '0.8.0', 'April 2021', '2020')
    description = 'Synteny-based scaffolding assessment and adjustment'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
#i# Descriptions of SynBad gap types.
gapdesc = {'Aln':'`Aligned` = Gap is found in the middle of a local alignment to the Hit',
           'Syn':'`Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).',
           'Ins':'`Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.',
           'Brk':'`Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.',
           'Dup':'`Duplication` = Overlapping flanking hits on the same strand.',
           'Inv':'`Inversion` = Flanking hits are on alternative strands.',
           'Tran':'`Translocation` = `SynSpan` indicates matches are on different scaffolds.',
           'Frag':'`Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.',
           'Term':'`Terminal` = Gap is between a local alignment and the end of the query sequence.',
           'Span':'`Spanned` = Any gaps without Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.',
           'HiC':'`HiC-best` = Any translocation gaps that have reciprocal best-rated HiC support.',
           'DupHiC':'`HiC-best Duplication` = Any duplication gaps that have reciprocal best-rated HiC support.',
           'Null':'No mapping between genomes for that gap.',
           'Div':'`Divergence` = Translocation or Breakpoint gap with flanking collinear hits.',
           'InvBrk':'`Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.',
           'InvFix':'`Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)',
           'InvDupFix':'`Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)',
           'Long':'`Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).'}
#########################################################################################################################
#i# Setup lists of gap types for use in different situations
#i# For now, InvFix are not in goodgaps to avoid nested inversions etc. Recommend re-running on fixed fasta.
goodgaps = ['Syn', 'Aln', 'Span', 'Div', 'Ins', 'Long','HiC','DupHiC']
skipgaps = ['Syn', 'Aln', 'Span', 'Long', 'Div', 'Dup', 'Ins','DupHiC','HiC']  #X# 'InvFix', 'InvDupFix',
fraggaps = ['Brk', 'Inv', 'InvBrk', 'Frag', 'Tran']
hicgaps = ['Brk', 'InvBrk', 'Frag', 'Tran', 'Null']
#########################################################################################################################
## SYNBAD Database Table Details
#########################################################################################################################
#i# Genertic database field formatting dictionary for use with all tables except where noted
dbformats = {}
for field in ['seqlen','start','end','Span0','AlnNum','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd']: dbformats[field] = 'int'
for field in ['gaplen','MaxFlank5','MaxFlank3','Span0','Span100','Span1000','Span5000','GapSpan','HitStart','HitEnd','Non','Alt','Trans','Inv','Dup']: dbformats[field] = 'int'
for field in ['score','ctglen','id1','id2','pairs']: dbformats[field] = 'int'
for field in ['wtscore','score']: dbformats[field] = 'num'
#########################################################################################################################
#i# Database table fields and keys
dbfields = {}

dbkeys = {}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class SynBad(rje_obj.RJE_Object):
    '''
    Class. Author: Rich Edwards (2015).

    Str:str
    - BAM1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
    - BAM2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
    - Chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
    - Chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
    - GABLAM=X        : Optional prefix for GABLAM search [defaults to $BASEFILE.map]
    - GapMode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    - Genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
    - Genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
    - HiCBAM1=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
    - HiCBAM2=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
    - MapFlanks1=FILE : Flanks fasta file from previous SynBad run for mapping genome 1 flank identifiers []
    - MapFlanks2=FILE : Flanks fasta file from previous SynBad run for mapping genome 2 flank identifiers []
    - HiCScore=X      : HiC scoring mode (score/wtscore) [wtscore]
    - HiCMode=X       : Pairwise HiC assessment scoring strategy (synbad/pure/rand/full) [synbad]
    - HiCDir1=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.qryflanks/]
    - HiCDir2=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.hitflanks/]
    - NewAcc1=X       : Scaffold name prefix for updated Genome 1 output [None]
    - NewAcc2=X       : Scaffold name prefix for updated Genome 2 output [None]
    - PAF1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    - PAF2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]

    Bool:boolean
    - BestPair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - Fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
    - FullMap=T/F     : Whether to abort if not all flanks can be mapped [True]
    - PureFlanks=T/F  : Whether to restrict gap flanks to pure contig sequence (True) or include good gaps (False) [True]
    - Update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]

    Int:integer
    - GapFlanks=INT   : Size of gap flank regions to output for HiC pairing analysis (0=off) [10000]
    - MaxOverlap=INT  : Maximum overlap (bp) of adjacent local hits to allow compression [500]
    - MaxSynSkip=INT  : Maximum number of local alignments to skip for SynTrans classification [4]
    - MaxSynSpan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [25000]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
    - MinReadSpan=INT     : Min number of Span0 reads in gaps table to prevent fragmentation [1]
    - SpannedFlank=INT : Required flanking distance for synreadspan "Spanned" rating [0]
    - SynReadSpan=INT : Minimum number of reads spanning a gap to change the rating to "Spanned" [5]

    Num:float
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]

    File:file handles with matching str filenames

    List:list
    - FragTypes=LIST  : List of SynBad ratings to trigger fragmentation [Brk,Inv,InvBrk,Frag,Tran]
    - Reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - Reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    - ReadType2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    - qry = Assembly Map list for Genome 1
    - hit = Assembly Map list for Genome 2

    Dict:dictionary
    - FlankMap = {qry:{GapFlank:OldFlank},hit:{GapFlank:OldFlank}}

    Obj:RJE_Objects
    - DB = Database object
    - Genome1 = Genome 1 SeqList Object
    - Genome2 = Genome 2 SeqList Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BAM1','BAM2','Chr1','Chr2','GABLAM','GapMode','Genome1','Genome2','HiCBAM1','HiCBAM2',
                        'MapFlanks1','MapFlanks2','HiCScore','HiCMode','HiCDir1','HiCDir2','NewAcc1','NewAcc2','PAF1','PAF2']
        self.boollist = ['BestPair','DocHTML','Fragment','FullMap','Update']
        self.intlist = ['GapFlanks','MaxOverlap','MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan','SpannedFlank','SynReadSpan']
        self.numlist = ['MinLocID']
        self.filelist = []
        self.listlist = ['FragTypes','Reads1','Reads2','ReadType1','ReadType2','qry','hit']
        self.dictlist = ['FlankMap']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'GapMode':'gapspan','HiCScore':'wtscore','HiCMode':'synbad'})
        self.setBool({'BestPair':False,'Fragment':False,'PureFlanks':True,'FullMap':True,'Update':False})
        self.setInt({'GapFlanks':10000,'MaxOverlap':500,'MaxSynSkip':4,'MaxSynSpan':25000,'MinLocLen':1000,'MinReadSpan':1,'SpannedFlank':0,'SynReadSpan':5})
        self.setNum({'MinLocID':50.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.list['FragTypes'] = fraggaps
        self.dict['FlankMap'] = {'qry':{}, 'hit':{}}
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ###
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ###
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['Chr1','Chr2','GABLAM','GapMode','HiCScore','HiCMode','NewAcc1','NewAcc2'])   # Normal strings
                self._cmdReadList(cmd,'path',['HiCDir1','HiCDir2'])  # String representing directory path
                self._cmdReadList(cmd,'file',['BAM1','BAM2','Genome1','Genome2','HiCBAM1','HiCBAM2','MapFlanks1','MapFlanks2','PAF1','PAF2'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BestPair','DocHTML','Fragment','FullMap','PureFlanks','Update'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['GapFlanks','MaxOverlap','MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan','SpannedFlank','SynReadSpan'])   # Integers
                self._cmdReadList(cmd,'perc',['MinLocID']) # Percentage
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['FragTypes','ReadType1','ReadType2'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Reads1','Reads2']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # SynBad:  Synteny-based scaffolding adjustment

        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.

        SynBad will use or create:

        1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
        spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
        correspond to misassemblies.

        2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.

        Pairwise hits between the genomes are filtered according to the `minlocid=PERC` and `minloclen=INT` criteria, which
        by default limits hits to be at least 1kb and 50% identity. Note that this is applied after GABLAM has run, so it
        should be possible to re-run with more relaxed constraints and re-use the GABLAM tables.

        Next, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
        sorted and `QryGap` and `SbjGap` fields added. Any local alignments with flanking hits are then flagged in these
        new fields with the number of flanking 5', 3' gaps.

        The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
        corresponding local hits on the Qry and Sbj genomes. If there is also an inversion, `SynSpan` will be negative.
        If the local hits are against two different sequences from the other genome, the two sequence names will be
        entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
        small gaps), `SynSpan` or `GapSpan` will have a value of zero.

        Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:

        * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
        * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).
        * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
        * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
        * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
        * `Inv` = `Inversion` = Flanking hits are on alternative strands.
        * `Tran` = `Translocation` = `SynSpan` indicates matches are on different scaffolds.
        * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.
        * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.
        * `Span` = `Spanned` = Any gaps without Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
        * `Null` = No mapping between genomes for that gap.

        If `chr1=X` and/or `chr2=X` chromosome scaffold prefixes are provided then `Translocation` will be restricted to
        matches between two different chromosome scaffolds. Matches including one or more non-chromosome scaffolds will
        be classed as `Fragmentation`.

        Following gap fixes for inversions and translocations/breakpoints (below), following

        After the initial mapping and GABLAM unique hit reduction, SynBad will further compress the Qry-Hit pairs by
        combining collinear hits between Qry-Hit pairs into larger blocks. Merged hits must be within `maxsynspan=INT`
        bp, overlap by no more than `maxoverlap=INT` bp, and have no intervening assembly gaps. An initial pass compares
        adjacent hits in position order. This is followed by iterative compressions where up to `maxsynskip=INT`
        intervening hits will be incorporated, providing there are no assembly gaps between the merged hits. This will
        start with the longest local hit, scanning downstream and then upstream for hits to merge, with increasing numbers
        of intervening hits up to `maxsynskip=INT`. Once no more merging is found for that hit, the next longest hit will
        be considered. If any hits are merged during this process, it will be repeated until convergence is reached.

        During merging, `Length` and `Identity` statistics will be summed, and additional statistics will be added:

        * `Non` = The unaligned bases in the query between the GABLAM local hits.
        * `Alt` = The number of bases (`Length`) of hits to different hit scaffolds
        * `Trans` = The number of bases (`Length`) of hits to the same scaffold but outside the merged hits.
        * `Inv` = The number of bases (`Length`) of intermediate hits that are in the right place but inverted strand.
        * `Dup` = The number of overlapping bases in the query scaffolds between adjacent hits.

        NOTE: This will be done for both genomes. For Genome 2, the "hits" are the `Qry` scaffolds, and all assessments
        of collinearity will be performed on the `Hit` scaffolds.

        Next, SynBad will assess inversions and apparent translocations for flanking collinear hits and instigate a
        number of updated ratings or suggested fixes. The following new SynBad gap ratings may be added:

        * `Div` = `Divergence` = Translocation or Breakpoint gap with flanking collinear hits.
        * `InvBrk` = `Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.
        * `InvFix` = `Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)
        * `InvDupFix` = `Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)
        * `Long` = `Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).

        If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
        `minreadspan=INT` reads span the gap.

        A future release of Synbad will optionally re-arrange the two assemblies, incorporating gapass assemblies
        where possible. [GapSpanner](https://github.com/slimsuite/gapspanner) can also be used to gap-fill spanned gaps
        prior to running SynBad.

        ## Dependencies

        SynBad needs Minimap2 installed. For `gapass` gap mode, Flye also needs to be installed. To generate
        documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For full documentation of the SynBad workflow, run with `dochtml=T` and read the `*.docs.html` file generated.


        ## Commandline options

        ```
        ### ~ Main SynBad run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
        genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
        basefile=X      : Prefix for output files [synbad]
        gablam=X        : Optional prefix for GABLAM search [defaults to $BASEFILE.map]
        gapmode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
        minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
        minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
        maxsynskip=INT  : Maximum number of local alignments to skip for SynTrans classification [4]
        maxsynspan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [25000]
        synreadspan=INT : Minimum number of reads spanning a gap to change the rating to "Spanned" [5]
        spannedflank=INT: Required flanking distance for synreadspan "Spanned" rating [0]
        maxoverlap=INT  : Maximum overlap (bp) of adjacent local hits to allow compression [500]
        chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
        chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
        ### ~ Fragmentation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
        fragtypes=LIST  : List of SynBad ratings to trigger fragmentation [Brk,Inv,InvBrk,Frag,Tran]
        minreadspan=INT : Min number of Span0 reads in gaps table to prevent fragmentation [1]
        ### ~ Additional input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        bam1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
        paf1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
        reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        bam2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
        paf2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]
        reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        mapflanks1=FILE : Flanks fasta file from previous SynBad run for mapping genome 1 flank identifiers []
        mapflanks2=FILE : Flanks fasta file from previous SynBad run for mapping genome 2 flank identifiers []
        fullmap=T/F     : Whether to abort if not all flanks can be mapped [True]
        ### ~ HiC Gap Flank options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        hicbam1=FILE    : Optional BAM file of HiC reads mapped onto assembly 1 [$BASEFILE1.HiC.bam]
        hicbam2=FILE    : Optional BAM file of HiC reads mapped onto assembly 2 [$BASEFILE2.HiC.bam]
        gapflanks=INT   : Size of gap flank regions to output for HiC pairing analysis (0=off) [10000]
        pureflanks=T/F  : Whether to restrict gap flanks to pure contig sequence (True) or include good gaps (False) [True]
        hicscore=X      : HiC scoring mode (score/wtscore) [wtscore]
        hicmode=X       : Pairwise HiC assessment scoring strategy (synbad/pure/rand/full) [synbad]
        hicdir1=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.qryflanks/]
        hicdir2=PATH    : Path to HiC read ID lists for genome 1 [$BASEFILE.hitflanks/]
        ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        newacc1=X       : Scaffold name prefix for updated Genome 1 output [None]
        newacc2=X       : Scaffold name prefix for updated Genome 2 output [None]
        bestpair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
        update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]
        force=T/F       : Whether to force regeneration of SynBad results tables [False]
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        ## SynBad workflow

        _Details of the SynBad workflow will be added in a future release._

        **NOTE:** SynBad is still under development. Features and details will continue to be updated.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return self.docHTML()
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.synBad()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if self.getStrLC('GapMode') not in ['gapspan','gapass','gapfill']:
                self.printLog('#MODE','GapMode "{0}" not recognised. Setting gapmode=gapspan.'.format(self.getStrLC('GapMode')))
                self.setStr({'GapMode':'gapspan'})
            if self.getStrLC('GapMode') not in ['gapspan','gapass']:
                self.printLog('#MODE','GapMode "{0}" no longer supported. Please run Diploidocus separately or use gapmode=gapspan.'.format(self.getStrLC('GapMode')))
                if rje.yesNo('Set gapmode=gapspan?'):
                    self.setStr({'GapMode':'gapspan'})
                    self.printLog('#MODE','Setting gapmode=gapspan.')
                else: return False
            ### ~ [2] Check for files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('Genome1')): raise IOError('Genome1 file "{0}" not found!'.format(self.getStr('Genome1')))
            if not rje.exists(self.getStr('Genome2')): raise IOError('Genome2 file "{0}" not found!'.format(self.getStr('Genome2')))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def docHTML(self):  ### Generate the Diploidocus Rmd and HTML documents.                                        # v0.1.0
        '''Generate the Diploidocus Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            info = self.log.obj['Info']
            prog = '%s V%s' % (info.program,info.version)
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
            rtxt += string.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2020 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown SynBad documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### Loading and formatting data methods                                                                     #
#########################################################################################################################
    def dbTable(self,tclass,ttype='',expect=True):     ### Returns tclass.ttype table, else None/Error
        '''
        Returns tclass.ttype table, else None/Error.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tname = '{0}.{1}'.format(tclass.lower(),ttype.lower())
            if not ttype: tname = tclass.lower()
            devtab = '{0}{1}'.format(tclass.lower(),ttype.lower())
            table = self.db(tname)
            if table: return table
            elif expect:
                if self.db(devtab):
                    self.warnLog('Using dev table name: "{0}"'.format(devtab))
                    return self.db(devtab)
                raise ValueError('Cannot find table: "{0}"'.format(tname))
            return None
        except:
            raise
#########################################################################################################################
    def seqObjSetup(self):  ### Loads the two genomes into sequence list objects
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:  ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Genome1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'SeqList1' not in self.obj or not self.obj['SeqList1']:
                base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
                cmd1 = ['seqin={0}'.format(self.getStr('Genome1')),
                        'runmode={0}'.format(self.getStrLC('GapMode')),
                        'basefile={0}'.format(base1)]
                seqcmd = self.cmd_list + cmd1 + ['summarise=F', 'gapstats=F', 'raw=F', 'dna=T']
                self.obj['SeqList1'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
            ## ~ [1b] Genome2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'SeqList2' not in self.obj or not self.obj['SeqList2']:
                base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
                cmd2 = ['seqin={0}'.format(self.getStr('Genome2')),
                        'runmode={0}'.format(self.getStrLC('GapMode')),
                        'basefile={0}'.format(base2)]
                seqcmd = self.cmd_list + cmd2 + ['summarise=F', 'gapstats=F', 'raw=F', 'dna=T']
                self.obj['SeqList2'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
            return
        except:
            self.errorLog('%s.seqObjSetup error' % self.prog())
#########################################################################################################################
    ### <4> ### Main SynBad Methods                                                                                     #
#########################################################################################################################
    #i# SynBad processing workflow:
    #i# [1] Load genomes = self.diploidocusGapRun()
    #i# [2] Extract gaps = self.diploidocusGapRun()
    #i# [3] GABLAM Search of Qry versus Hit = self.runGABLAM()
    #i# [4] Use the genome mapping to rate all the gaps = self.synBadMapping(cdb1,cdb2)
    #i# [5] Generate (or map) flanks and map on the HiC data if available = self.synBadHiCMap(): return False
        #!# -> Add the update of the HiC score support to update like the 'Div' updates
    #i# [6]
    # if not self.synBadCompress(): return False
    # if not self.tidyTables(): return False
    # if not self.saveTables(): return False
    #     def makeAssemplyMap(self): return
    # if not self.corrections(): return False
    # return self.fragment()
    # def synBadSummarise(self): return
#########################################################################################################################
    def synBad(self):   ### Main SynBad run method
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#
            #!# Add possibility to run without a second genome - pure read-spanning and HiC verification

            ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GAPSPANNER GAP ANALYSIS', line='=')
            self.infoLog('SynBad uses GapSpanner to extract a table of assembly gaps from each input assembly.')
            self.infoLog('If long-read sequencing data is provided, gaps will also be assessed for spanning read support.')
            (cdb1,cdb2) = self.diploidocusGapRun()
            if not cdb1 or not cdb2: return False

            ### ~ [2] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GENOME VS GENOME GABLAM SEARCH',line='=')
            self.infoLog('GABLAM is used to search Genome1 against Genome2. Local alignments are filtered by length and identity.')
            self.infoLog('For each genome, local hits are reduced to unique (non-overlapping) alignments with the other genome.')
            self.infoLog('This unique hit reduction is performed based on number of identical aligned based.')
            self.infoLog('Partially overlapping local hits are trimmed.')
            if not self.runGABLAM(): return False

            ### ~ [3] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.synBadMapping(cdb1,cdb2): return False

            #!# Will want to split this method up, but potentially also move the mapping etc. higher up and classification later
            if not self.synBadHiCMap(): return False
            self.headLog('ASSEMBLY MAPS',line='=')
            for qh in ['qry','hit']:
                #!# Add checking and reloading from OLDmakeAssemplyMap()
                if not self.makeAssemplyMap(qh): return False
                if not self.saveAssemblyMaps(qh,mapname='map'): return False

            if not self.synBadCompress(): return False
            if not self.tidyTables(): return False
            if not self.saveTables(backup=False): return False
            if not self.synBadSummarise(): return False


            #!# Move corrections and fragmentation outside this method and make self-sufficient (loading tables if missing)

            ### ~ [4] Update assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.mapCorrections(): return False

            ### ~ [5] Fragment Gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GAP FRAGMENTATION',line='=')
            # If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
            # `minreadspan=INT` reads span the gap.
            if self.getBool('Fragment'):
                return self.fragment()
            else:
                self.printLog('#FRAG','No fragmentation (fragment=F).')

            # A future release of Synbad will optionally re-arrange the two assemblies, incorporating gapass assemblies
            # where possible.

            return True
        except: self.errorLog('%s.synBad error' % self.prog()); return False
#########################################################################################################################
    ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    def diploidocusGapRun(self):    ### Runs Diploidocus and returns the gap tables
        '''
        Runs Diploidocus and returns the gap tables
        :return: (cdb1,cdb2)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            dcmd = ['bam=', 'paf=', 'reads=', 'readtype=ont']
            wanted = ['qry.gap.tdt','hit.gap.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad *.gap.tdt files found: skipping GapSpanner extraction (force=F)')
                cdb1 = db.addTable('%s.qry.gap.tdt' % db.baseFile(), mainkeys=['seqname', 'start', 'end'], name='qry.gap', ignore=[], expect=True)
                cdb1.dataFormat(dbformats)
                cdb2 = db.addTable('%s.hit.gap.tdt' % db.baseFile(), mainkeys=['seqname', 'start', 'end'], name='qry.gap', ignore=[], expect=True)
                cdb2.dataFormat(dbformats)
                return (cdb1,cdb2)

            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#
            #!# Update this to cycle qry and hit as for other methods

            ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # SynBad will use or create:
            #
            # 1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
            # spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
            # correspond to misassemblies.
            #
            # ==> tigersnake.v2.7.pafscaff.checkpos.tdt <==
            # #       seqname start   end     seqlen  gaplen  gap     MaxFlank5       MaxFlank3       Span0   Span100 Span1000        Span5000
            # 209     NSCUCHR1.01     57638   57737   341225390       100     NSCUCHR1.01.57638-57737 57637   341167653       1       1       1       0
            # 2       NSCUCHR1.01     103648  104147  341225390       500     NSCUCHR1.01.103648-104147       103647  341121243       14      14      12      5

            ## ~ [1a] Genome1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
            cmd1 = ['seqin={0}'.format(self.getStr('Genome1')),
                    'runmode={0}'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base1)]
            rundip = False
            if self.getStrLC('BAM1'): cmd1.append('bam={0}'.format(self.getStr('BAM1')))
            if self.getStrLC('PAF1'): cmd1.append('paf={0}'.format(self.getStr('PAF1'))); rundip = True
            if self.list['Reads1']: cmd1.append('reads={0}'.format(','.join(self.list['Reads1']))); rundip = True
            if self.list['ReadType1']: cmd1.append('readtype={0}'.format(','.join(self.list['ReadType1'])))
            dip1 = diploidocus.Diploidocus(self.log, self.cmd_list + dcmd + cmd1)
            cdb1 = None
            if rundip:
                dip1.run()
                cdb1 = dip1.db('checkpos')
                if cdb1:
                    cdb1.baseFile(basefile)
                    cdb1.setStr({'Name': 'qry.gap'})
                    self.db().list['Tables'].append(cdb1)
                else:
                    cdb1 = db.addTable('%s.checkpos.tdt' % base1, mainkeys=['seqname', 'start', 'end'], name='qry.gap', ignore=[], expect=True)
            if not cdb1:
                if rundip:
                    self.errorLog('Diploidocus checkpos run appears to have failed', printerror=False, quitchoice=self.i() >= 0)
                else:
                    self.printLog('#READS', 'No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base1))
                if not rje.checkForFiles(filelist=['.gaps.tdt'], basename=base1, log=self.log):
                    seqcmd = self.cmd_list + cmd1 + ['summarise=T', 'gapstats=T', 'raw=F', 'dna=T']
                    self.obj['SeqList1'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
                cdb1 = db.addTable('%s.gaps.tdt' % base1, mainkeys=['seqname', 'start', 'end'], name='qry.gap', ignore=[], expect=True)
                cdb1.addField('Span0', evalue=0)
            cdb1.dataFormat({'seqlen': 'int', 'start': 'int', 'end': 'int', 'Span0': 'int'})
            cdb1.addField('GapSpan', evalue='.')
            cdb1.addField('SynSpan', evalue='.')
            cdb1.addField('SynBad', evalue='Null')
            cdb1.index('seqname')
            if '#' in cdb1.fields(): cdb1.dropField('#')

            ## ~ [1b] Genome2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
            cmd2 = ['seqin={0}'.format(self.getStr('Genome2')),
                    'runmode={0}'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base2)]
            rundip = False
            if self.getStrLC('BAM2'): cmd2.append('bam={0}'.format(self.getStr('BAM2')))
            if self.getStrLC('PAF2'): cmd2.append('paf={0}'.format(self.getStr('PAF2'))); rundip = True
            if self.list['Reads2']: cmd2.append('reads={0}'.format(','.join(self.list['Reads2']))); rundip = True
            if self.list['ReadType2']: cmd2.append('readtype={0}'.format(','.join(self.list['ReadType2'])))
            dip2 = diploidocus.Diploidocus(self.log, self.cmd_list + dcmd + cmd2)
            cdb2 = None
            if rundip:
                dip2.run()
                cdb2 = dip2.db('checkpos')
                if cdb2:
                    cdb2.baseFile(basefile)
                    cdb2.setStr({'Name': 'hit.gap'})
                    self.db().list['Tables'].append(cdb2)
                else:
                    cdb2 = db.addTable('%s.checkpos.tdt' % base2, mainkeys=['seqname', 'start', 'end'], name='hit.gap', ignore=[], expect=True)
            if not cdb2:
                if rundip:
                    self.errorLog('Diploidocus checkpos run appears to have failed', printerror=False, quitchoice=self.i() >= 0)
                else:
                    self.printLog('#READS', 'No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base2))
                if not rje.checkForFiles(filelist=['.gaps.tdt'], basename=base2, log=self.log):
                    seqcmd = self.cmd_list + cmd2 + ['summarise=T', 'gapstats=T', 'raw=F', 'dna=T']
                    self.obj['SeqList2'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
                cdb2 = db.addTable('%s.gaps.tdt' % base2, mainkeys=['seqname', 'start', 'end'], name='hit.gap', ignore=[], expect=True)
                cdb2.addField('Span0', evalue=0)
            cdb2.dataFormat({'seqlen': 'int', 'start': 'int', 'end': 'int', 'Span0': 'int'})
            cdb2.addField('GapSpan', evalue='.')
            cdb2.addField('SynSpan', evalue='.')
            cdb2.addField('SynBad', evalue='Null')
            cdb2.index('seqname')
            if '#' in cdb2.fields(): cdb2.dropField('#')

            return(cdb1,cdb2)
        except: self.errorLog('%s.diploidocusGapRun error' % self.prog()); raise
#########################################################################################################################
    ### ~ [2] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    def runGABLAM(self):   ### Runs GABLAM of two assemblies
        '''
        Runs GABLAM of two assemblies and loads results into database tables.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            wanted = ['qry.tdt','hit.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad qry.tdt and hit.tdt files found: skipping GABLAM (force=F)')
                return True
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#

            ### ~ [1] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gabbase = basefile + '.map'
            if self.getStrLC('GABLAM'): gabbase = self.getStr('GABLAM')
            self.printLog('#GABLAM','GABLAM output basefile: {0}'.format(gabbase))
            # 2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.
            # ==> tiger.v.najna.XXXunique.tdt <==
            # Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd
            # In each case, Qry is Genome1 and Hit is Genome2
            quniq = '{0}.qryunique.tdt'.format(gabbase)
            huniq = '{0}.hitunique.tdt'.format(gabbase)
            ## ~ [2a] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.force() or not rje.exists(quniq) or not rje.exists(huniq):
                gabcmd = ['seqin={0}'.format(self.getStr('Genome1')),'searchdb={0}'.format(self.getStr('Genome2')),'mapper=minimap','minlocid=0','minloclen=0','basefile={0}'.format(gabbase)]
                gabobj = gablam.GABLAM(self.log,self.cmd_list+gabcmd)
                gabobj.gablam()
            ## ~ [2b] Load tables and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                ufile = '{0}.{1}unique.tdt'.format(gabbase,qh.lower())
                udb = db.addTable(ufile,mainkeys=['Qry','Hit','AlnNum'],name=qh.lower(),ignore=[],expect=True)
                #udb = db.addTable(ufile,mainkeys=[qh,'{0}Start'.format(qh),'{0}End'.format(qh)],name=qh.lower(),ignore=[],expect=True)
                udb.dataFormat({'AlnNum':'int','Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
                lenx = 0; idx = 0
                for entry in udb.entries():
                    if entry['Length'] < self.getInt('MinLocLen'): udb.dropEntry(entry); lenx += 1
                    elif (100.0 * entry['Identity'] / entry['Length']) < self.getNum('MinLocID'): udb.dropEntry(entry); idx += 1
                self.printLog('#MINCUT','Dropped %s entries < %s bp and %s < %.1f%% identity' % (rje.iStr(lenx),rje.iStr(self.getInt('MinLocLen')),rje.iStr(idx),self.getNum('MinLocID')))
                udb.addField('HitStart',evalue=0)
                udb.addField('HitEnd',evalue=0)
                udb.addField('Strand',evalue='+')
                for entry in udb.entries():
                    if entry['SbjStart'] > entry['SbjEnd']:
                        entry['HitStart'] = entry['SbjEnd']
                        entry['HitEnd'] = entry['SbjStart']
                        entry['Strand'] = '-'
                    else:
                        entry['HitStart'] = entry['SbjStart']
                        entry['HitEnd'] = entry['SbjEnd']
                udb.newKey([qh,'{0}Start'.format(qh),'{0}End'.format(qh)])
                udb.setFields(['Qry','QryStart','QryEnd','Hit','HitStart','HitEnd','Strand','AlnNum','Length','Identity'])
                udb.addField('Non',evalue=0)
                udb.addField('QryGap',evalue='')
                udb.addField('HitGap',evalue='')
                udb.addField('Qry5',evalue='')
                udb.addField('Qry3',evalue='')
                udb.addField('Hit5',evalue='')
                udb.addField('Hit3',evalue='')

            return True
        except: self.errorLog('%s.runGABLAM error' % self.prog()); return False
#########################################################################################################################
    def synBadMapping(self,cdb1,cdb2):   ### Main SynBad gap mapping and compression
        '''
        Main SynBad gap mapping and compression.
        >> cdb1: Table of genome 1 gaps
        >> cdb2: Table of genome 2 gaps
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD GAP MAPPING',line='=')
            base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
            base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
            self.printLog('#NOTE','Genome1 ({0}) is always Qry and Genome2 ({1}) is always Hit'.format(base1,base2))
            db = self.db()
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            wanted = ['qry.tdt','hit.tdt','qry.gap.tdt','hit.gap.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad qry.tdt, hit.tdt and *.gap.tdt files found: skipping SynBad mapping (force=F)')
                return True

            ### ~ [1] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # First, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
            # sorted and `QryGap` and `HitGap` fields added. Any local alignments with flanking hits are then flagged in these
            # new fields with the number of flanking 5', 3' gaps indicated by < and > characters.
            #
            # The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
            # corresponding local hits on the Qry and Hit genomes. If there is also an inversion, `SynSpan` will be negative.
            # If the local hits are against two different sequences from the other genome, the two sequence names will be
            # entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
            # small gaps), `SynSpan` or `GapSpan` will have a value of zero.
            ## ~ [1a] Work through query and hit tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug(self.db().tableNames())
            for qh in ('Qry','Hit'):
                qry = qh
                hit = {'Qry':'Hit','Hit':'Qry'}[qry]
                hitchr = {'Qry':chr1,'Hit':chr2}[hit]
                # gapdb = self.db('{0}gap'.format(qh.lower()))
                gapdb = self.dbTable(qh,'gap')
                spanfield = 'Span{0:d}'.format(self.getInt('SpannedFlank'))
                spancheck = spanfield in gapdb.fields()
                if self.getInt('SynReadSpan') > 0 and not spancheck:
                    self.warnLog('synreadspan={0} but "{0}" not found in fields. Check spannedflank=INT.'.format(self.getInt('SynReadSpan'),spanfield))
                spancheck = spancheck and self.getInt('SynReadSpan') > 0
                if spancheck: gapdb.dataFormat({spanfield:'int'})
                #self.debug(gapdb)
                altdb = {cdb1:cdb2,cdb2:cdb1}[gapdb]
                # locdb = self.db(qh.lower())
                locdb = self.dbTable(qh)
                locdb.index(qry)
                locdb.index(hit)
                # First, deal with QryGaps
                for seqname in gapdb.indexKeys('seqname'):
                    if seqname not in locdb.index(qry):
                        self.printLog('#ORPHAN','No alignments found for {0}'.format(seqname))
                        continue
                    self.progLog('\r#SYNBAD','SynBad mapping: {0}'.format(seqname))
                    gap = gapdb.index('seqname')[seqname][0:]
                    loc = [(seqname,0,0)] + locdb.index(qry)[seqname] + [(seqname,-1,-1)]
                    while gap and len(loc) > 1 :
                        # Cycle to get the first gap in the list between two local alignment entries
                        while gap and gap[0][1] <= loc[0][2]: gap.pop(0)
                        while len(loc) > 1 and loc[1][2] < gap[0][1]: loc.pop(0)
                        if not gap or not len(loc) > 1: continue
                        # Should now have loc[0][2] < gap[0][1] & loc[1][2] > gap[0][1]
                        # Also want to have gap[0][2] < loc[1][1] for a gap properly flanked by two local alignments
                        entry5 = None
                        entry3 = None
                        if loc[0][1] > 0:
                            entry5 = locdb.data(loc[0])
                            entry5['{0}Gap'.format(qry)] += '>'
                            if not entry5['{0}3'.format(qry)]: entry5['{0}3'.format(qry)] = []
                            entry5['{0}3'.format(qry)].append(gap[0])
                        if loc[1][1] > 0:
                            entry3 = locdb.data(loc[1])
                            entry3['{0}Gap'.format(qry)] += '<'
                            if not entry3['{0}5'.format(qry)]: entry3['{0}5'.format(qry)] = []
                            entry3['{0}5'.format(qry)].append(gap[0])
                        thisgap = gap.pop(0)
                        gentry = gapdb.data(thisgap)
                        if entry3 and thisgap[2] >= loc[1][1]: # Gap overlapping local alignment
                            gentry['GapSpan'] = 0
                            gentry['SynSpan'] = '{0}:0'.format(entry3[hit])
                            gentry['SynBad'] = 'Aligned'
                            continue
                        # At this point, loc[0][2] < thisgap[1] & thisgap[2] < loc[1][1]
                        if loc[1][1] == -1: gentry['GapSpan'] = gentry['seqlen'] - loc[0][2]
                        else: gentry['GapSpan'] = loc[1][1] - loc[0][2]
                        if not entry5 or not entry3:
                            gentry['SynBad'] = 'Terminal'
                            if spancheck and gentry[spanfield] >= self.getInt('SynReadSpan'):
                                gentry['SynBad'] = 'Spanned'
                            continue
                        # Assess synteny and classify
                        hitstart = '{0}Start'.format(hit)
                        hitend = '{0}End'.format(hit)
                        gap5 = {'+':entry5[hitend],'-':entry5[hitstart]}[entry5['Strand']]
                        gap3 = {'-':entry3[hitend],'+':entry3[hitstart]}[entry3['Strand']]
                        synspan = max(gap5,gap3) - min(gap5,gap3)
                        #i# Look for picking up synteny by skipping fragments from other sequences
                        #i# Where a local alignment switches the hit sequence, SynBad will look downstream for another
                        #i# aligment to the same sequence. A maximum of `maxsynskip=INT` alignments will be skipped, up to
                        #i# a maximum distance of `maxsynspan=INT` bp.
                        maxskip = self.getInt('MaxSynSkip')     # Maximum number of local alignments to skip for SynTrans classification
                        entryskip = None
                        egap = -1
                        if entry3[hit] != entry5[hit]:
                            for i in range(maxskip):
                                if len(loc) > (2+i) and loc[2+i][1] > 0:
                                    entryskip = locdb.data(loc[2+i])
                                    if not entryskip:
                                        self.warnLog('Problem finding local hit entry for {0}'.format(str(loc[i+2])))
                                        continue
                                    goodskip = True
                                    for field in [hit,'Strand']:
                                        if entryskip[field] != entry5[field]: goodskip = False
                                    if entry5['Strand'] == '+' and entryskip[hitstart] < entry5[hitend]: goodskip = False
                                    elif entry5['Strand'] == '-' and entryskip[hitend] > entry5[hitstart]: goodskip = False
                                    elif entry5['Strand'] == '+': egap = entryskip[hitstart] - entry5[hitend]
                                    elif entry5['Strand'] == '-': egap = entryskip[hitend] - entry5[hitstart]
                                    if egap > self.getInt('MaxSynSpan'): goodskip = False
                                    if loc[2+i][1] - loc[0][2] > self.getInt('MaxSynSpan'): goodskip = False
                                    if goodskip: break
                                    entryskip = None
                        # Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:
                        #
                        # * `Aligned` = Gap is found in the middle of a local alignment to the Hit
                        # * `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
                        # * `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
                        # * `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
                        # * `Duplication` = Overlapping flanking hits on the same strand.
                        # * `Inversion` = Flanking hits are on alternative strands.
                        # * `Translocation` = `SynSpan` indicates matches are on different scaffolds.
                        # * `Terminal` = Gap is between a local alignment and the end of the query sequence.
                        # * `Null` = No mapping between genomes for that gap.
                        if entryskip:
                            gentry['SynSpan'] = '{0}:{1}:{2}'.format(entry5[hit],egap,entryskip[hit])
                            gentry['SynBad'] = 'Insertion'
                        elif entry3[hit] != entry5[hit]:
                            gentry['SynSpan'] = '{0}::{1}'.format(entry5[hit],entry3[hit])
                            if hitchr and (not entry5[hit].startswith(hitchr) or not entry3[hit].startswith(hitchr)):
                                gentry['SynBad'] = 'Fragmentation'
                            else:
                                gentry['SynBad'] = 'Translocation'
                        elif entry3['Strand'] != entry5['Strand']:
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Inversion'
                        elif (entry3['Strand'] == '+' and gap3 < gap5) or (entry3['Strand'] == '-' and gap3 > gap5):
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],-synspan)
                            gentry['SynBad'] = 'Duplication'
                        elif (synspan - gentry['GapSpan']) > self.getInt('MaxSynSpan'):
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Breakpoint'
                        else:
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Syntenic'
                        if spancheck and gentry[spanfield] >= self.getInt('SynReadSpan') and gentry['SynBad'] not in ['Aligned','Syntenic']:
                            gentry['SynBad'] = 'Spanned'

                # Shorten classes
                devswap = {'Aligned':'Aln','Syntenic':'Syn','Insertion':'Ins','Breakpoint':'Brk','Duplication':'Dup',
                           'Inversion':'Inv','Translocation':'Tran','Fragmentation':'Frag','Terminal':'Term','Spanned':'Span','Null':'Null'}
                #!# Make this an option? if self.dev():
                for entry in gapdb.entries():
                    if entry['SynBad'] in devswap: entry['SynBad'] = devswap[entry['SynBad']]

                # Then, deal with Hit Gaps
                for seqname in altdb.indexKeys('seqname'):
                    if seqname not in locdb.index(hit):
                        continue
                    self.progLog('\r#SYNBAD','SynBad mapping: {0}'.format(seqname))
                    gap = altdb.index('seqname')[seqname][0:]
                    loc = [(seqname,0,0)]
                    for entry in locdb.indexEntries(hit,seqname):
                        loc.append((entry[hit],entry['{0}Start'.format(hit)],entry['{0}End'.format(hit)],entry))
                    loc.sort()
                    loc.append((seqname,-1,-1))
                    while gap and len(loc) > 1 :
                        # Cycle to get the first gap in the list between two local alignment entries
                        while gap and gap[0][1] <= loc[0][2]: gap.pop(0)
                        while len(loc) > 1 and loc[1][2] < gap[0][1]: loc.pop(0)
                        if not gap or not len(loc) > 1: continue
                        # Should now have loc[0][2] < gap[0][1] & loc[1][2] > gap[0][1]
                        # Also want to have gap[0][2] < loc[1][1] for a gap properly flanked by two local alignments
                        thisgap = gap.pop(0)
                        if loc[0][1] > 0:
                            loc[0][3]['{0}Gap'.format(hit)] += '>'
                            if not loc[0][3]['{0}3'.format(hit)]: loc[0][3]['{0}3'.format(hit)] = []
                            loc[0][3]['{0}3'.format(hit)].append(thisgap)
                        if loc[1][1] > 0:
                            loc[1][3]['{0}Gap'.format(hit)] += '<'
                            if not loc[1][3]['{0}5'.format(hit)]: loc[1][3]['{0}5'.format(hit)] = []
                            loc[1][3]['{0}5'.format(hit)].append(thisgap)
                self.progLog('\r#SYNBAD','{0} SynBad mapping complete.               '.format(qh))
                self.printLog('\r#SYNBAD','{0} SynBad mapping complete.'.format(qh))

            ### ~ [2] Update tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qrygapdb = self.dbTable('qry','gap')
            hitgapdb = self.dbTable('hit','gap')
            qrydb = self.db('qry')
            hitdb = self.db('hit')
            ex = 0.0; etot = qrydb.entryNum() + hitdb.entryNum()
            for entry in qrydb.entries() + hitdb.entries():
                self.progLog('\r#UPDATE','Updating flanking gap SynBad ratings: %.2f%%' % (ex/etot)); ex += 100
                for field in ('Qry5','Qry3'):
                    if entry[field]:
                        gaps = []
                        for gap in entry[field]: gaps.append(qrygapdb.data(gap)['SynBad'])
                        entry[field] = ';'.join(rje.sortUnique(gaps))
                    else: entry[field] = '-'
                for field in ('Hit5','Hit3'):
                    if entry[field]:
                        gaps = []
                        for gap in entry[field]: gaps.append(hitgapdb.data(gap)['SynBad'])
                        entry[field] = ';'.join(rje.sortUnique(gaps))
                    else: entry[field] = '-'
                if not entry['QryGap']: entry['QryGap'] = '.'
                if not entry['HitGap']: entry['HitGap'] = '.'
            self.printLog('\r#UPDATE','Updating flanking gap SynBad ratings complete.')

            return True
        except: self.errorLog('%s.synBadMapping error' % self.prog()); return False
#########################################################################################################################
    def synBadCompressTable(self,qh,locdb,maxskip=0,quick=False):   ### Main SynBad gap mapping compression for a local hits table.
        '''
        Main SynBad gap mapping compression for a local hits table.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            if 'Alt' not in locdb.fields():
                locdb.addField('Alt',evalue=0,after='Non')
                locdb.addField('Trans',evalue=0,after='Alt')
                locdb.addField('Inv',evalue=0,after='Trans')
                locdb.addField('Dup',evalue=0,after='Inv')
            tname = locdb.name()
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: quick = True
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = []

            ### ~ [1] Cycle through entries, identifying and merging flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            merged = True
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            while entries:
                self.progLog('\r#CMPRSS','Compress %s table by %s: %.2f%%' % (tname,qh,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if not merged and not lookback:
                        lookback = True
                    elif not merged:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                elif not merged:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break
                elif eskip > 0:
                    # if merging was performed for higher eskip -> will need to reset eskip and cycle again
                    eskip = 0
                    lookback = False
                merged = False

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\nSkip={0}; Lookback={1}'.format(eskip,lookback))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.': continue #i# Don't merge gaps
                #i# - (2) Lack of internal qry gaps
                if '>' in merge1[qgap] or '<' in merge2[qgap]:
                    if debug: self.bugPrint('QryGap...')
                    continue
                skipgap = False
                for entry in toskip:
                    if entry[qgap] != '.': skipgap = True
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                #i# - (3) Lack of hit gaps
                #?# Currently eliminating any internal entries with gaps. May want to refine this.
                for entry in [merge1,merge2] + toskip:
                    if entry[hgap] != '.': skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                if qdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                if hdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                    continue

                ## ~ [1c] Classify and merge intermediates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Next, work through, classify and merge each intermediate entry and merge2
                if debug: self.bugPrint('-> merge ->')
                merged = True
                mstart = merge1[hend] - self.getInt('MaxOverlap')
                mend = merge2[hstart] + self.getInt('MaxOverlap')
                if bwd:
                    mstart = merge2[hend] - self.getInt('MaxOverlap')
                    mend = merge1[hstart] + self.getInt('MaxOverlap')
                prev = merge1
                merge1['Dup'] += max(0,-hdist)
                for entry in toskip + [merge2]:
                    ex += 100.0
                    if entry == merge2: etype = 'Merge'
                    #i# (i) Different sequence (AltChrom)
                    elif entry[alt] != merge1[alt]: etype = 'Alt'
                    #i# (i) Same seq, collinear
                    elif entry[hstart] > mstart and entry[hend] < mend:
                        if entry['Strand'] == merge1['Strand']: etype = 'Col'
                    #i# (ii) Same seq, inside, inverted (Inv)
                        else: etype = 'Inv'
                    #i# (iii) Same seq, outside (Trans)
                    else: etype = 'Trans'
                    #i# Merge into merge1
                    for field in ['Length','Identity','Alt','Trans','Inv','Dup']:
                        merge1[field] += entry[field]
                    inslen = entry[hend] - entry[hstart] + 1
                    if etype in ['Alt','Trans','Inv']: merge1[etype] += inslen
                    merge1['Non'] += max(0,entry[qstart]-prev[qend]-0)
                    entries.remove(entry)
                    if entry in sorted: sorted.remove(entry)
                    prev = entry

                ## ~ [1d] Update and remove the merged entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                merge1[qend] = merge2[qend]
                if bwd: merge1[hstart] = merge2[hstart]
                else: merge1[hend] = merge2[hend]
                merge1['{0}3'.format(qh)] = merge2['{0}3'.format(qh)]
                if merge1[qgap] != '.' and merge2[qgap] != '.':
                    merge1[qgap] = merge1[qgap] + merge2[qgap]
                elif merge2[qgap] != '.':
                    merge1[qgap] = merge2[qgap]
                locdb.dropEntryList(toskip + [merge2],log=False)
                if debug: self.debug('%s' % (locdb.entrySummary(merge1,collapse=True)))
            self.printLog('\r#CMPRSS','Compressed %s table by %s: %s alignments -> %s' % (tname,qh,rje.iStr(etot),rje.iStr(locdb.entryNum())))
            locdb.remakeKeys()

            return (etot - locdb.entryNum())
        except: self.errorLog('%s.synBadCompressTable error' % self.prog())
#########################################################################################################################
    def synBadCompress(self):   ### Main SynBad gap mapping compression
        '''
        Main SynBad gap mapping compression.
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            self.headLog('SYNBAD HIT COMPRESSION',line='=')
            self.printLog('#CMPRSS','Maximum distance between local hits for compression: %s' % rje_seqlist.dnaLen(self.getInt('MaxSynSpan')))
            self.printLog('#CMPRSS','Maximum overlap between local hits for compression: %s' % rje_seqlist.dnaLen(self.getInt('MaxOverlap')))
            self.printLog('#CMPRSS','Maximum number of local hits to bypass for compression: %s' % self.getInt('MaxSynSkip'))
            dformats = {}
            for field in ['start','end','gaplen','seqlen','gaplen','MaxFlank5','MaxFlank3','Span0','Span100','Span1000',
                          'Span5000','GapSpan','QryStart','QryEnd','HitStart','HitEnd','Length','Identity','Non',
                          'Alt','Trans','Inv','Dup']:
                dformats[field] = 'int'
            ## ~ [0a] Check for tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #X# wanted = ['qryfull.tdt','qry.tdt','qrypairs.tdt','qrygap.tdt','hitfull.tdt','hit.tdt','hitpairs.tdt','hitgap.tdt']
            #i# Removing pairs from the check: not used for anything!
            wanted = ['qry.full.tdt','qry.tdt','qry.gap.tdt','hit.full.tdt','hit.tdt','hit.gap.tdt']
            if rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                if self.force():
                    self.printLog('#FORCE','SynBad output files found, but regenerating (force=T).')
                else:
                    for tdtfile in wanted:
                        tname = tdtfile[:-4]
                        if 'full' in tname: continue
                        if self.getBool('Update'):
                            if self.db(tname): db.deleteTable(tname)
                            if 'pairs' in tname: continue
                        if 'gap' in tname:
                            tkeys = ['seqname','start','end']
                        elif tname.startswith('qry'):
                            tkeys = ['Qry','QryStart','QryEnd']
                            if 'pairs' in tname: tkeys.insert(1,'Hit')
                        else:
                            tkeys = ['Hit','HitStart','HitEnd']
                            if 'pairs' in tname: tkeys.insert(1,'Qry')
                        table = db.addTable(mainkeys=tkeys,name=tname,ignore=[],expect=True,replace=True)
                        table.dataFormat(dformats)
                    for qh in ('Qry','Hit'):
                        self.updateHitGHaps(qh)
                    if self.getBool('Update'):
                        self.printLog('#UPDATE','SynBad output files found and loaded, but updating (force=F update=T).')
                    else:
                        self.printLog('#SKIP','SynBad output files found and loaded (force=F update=F).')
                        return True
            else:
                self.printLog('#CHECK','Not all SynBad output files found: generating')
            ## ~ [0b] Save full tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                tname = qh.lower()
                locdb = self.db(tname)
                if not locdb:
                    tkeys = ['Hit','HitStart','HitEnd']
                    if tname.startswith('qry'): tkeys = ['Qry','QryStart','QryEnd']
                    locdb = db.addTable(mainkeys=tkeys,name=tname,ignore=[],expect=True,replace=True)
                    locdb.dataFormat(dformats)
                #!# Don't output if qry and qryfull both found?
                locdb.info['Name'] = '{0}.full'.format(qh.lower())
                locdb.saveToFile()
                locdb.info['Name'] = qh.lower()
            ## ~ [0c] Add gaps to tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                alt = {'Qry':'Hit','Hit':'Qry'}[qh]
                gapdb = self.dbTable(qh,'gap')
                gapdb.dataFormat({'start':'int','end':'int','gaplen':'int'})
                for gentry in gapdb.entries():
                    if gentry['SynBad'] in ['Aln','Aligned']: continue
                    #  seqname start end seqlen gaplen gap MaxFlank5  MaxFlank3  Span0  Span100  Span1000  Span5000  GapSpan  SynSpan SynBad
                    # Qry           QryStart   QryEnd     Hit           HitStart   HitEnd     Strand  Length   Identity  Non  Alt    Trans   Inv    QryGap  HitGap  Qry5       Qry3      Hit5       Hit3
                    addentry = {qh:gentry['seqname'],'{0}Start'.format(qh):gentry['start'],'{0}End'.format(qh):gentry['end'],
                                alt:gentry['SynSpan'],'{0}Start'.format(alt):1,'{0}End'.format(alt):gentry['gaplen'],
                                '{0}Gap'.format(qh):gentry['SynBad'],'{0}Gap'.format(alt):'.','Strand':'.',
                                'Non':gentry['gaplen'],'Alt':0,'Trans':0,'Inv':0,'Dup':0,
                                'Length':gentry['gaplen'],'Identity':0,'Qry5':'-','Qry3':'-','Hit5':'-','Hit3':'-'}
                    locdb.addEntry(addentry,overwrite=False)
                self.printLog('#GAPS','Added {0} {1} gaps to {2} table -> {3} entries'.format(rje.iStr(gapdb.entryNum()),qh,qh.lower(),rje.iStr(locdb.entryNum())))

            ### ~ [1] Quick merge of adjacent flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ('Qry','Hit'):
                locdb = self.dbTable(qh)
                skipx = 0
                self.headLog('SynBad {0} quick compression: merge adjacent local hits'.format(qh),line='~')
                self.synBadCompressTable(qh,locdb,skipx,quick=True)

            ### ~ [2] Remove internal odd entries and merge flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('MaxSynSkip') > 0:
                for qh in ('Qry','Hit'):
                    alt = {'Qry':'Hit','Hit':'Qry'}[qh]
                    locdb = self.dbTable(qh)
                    cyc = 0
                    prex = 0
                    while locdb.entryNum() != prex:
                        prex = locdb.entryNum(); cyc += 1
                        # for skipx in range(self.getInt('MaxSynSkip'),-1,-1):
                        #     self.headLog('SynBad {0} compression cycle {1} (SynSkip={2})'.format(qh,cyc,skipx),line='~')
                        #     self.synBadCompressTable(qh,locdb,skipx)
                        self.headLog('SynBad {0} compression cycle {1} (MaxSynSkip={2})'.format(qh,cyc,self.getInt('MaxSynSkip')),line='~')
                        self.synBadCompressTable(qh,locdb,self.getInt('MaxSynSkip'))
                    continue
            else: self.printLog('#MERGE','Adjacent merging only (maxsynskip=0).')

            return True
        except: self.errorLog('%s.synBadCompress error' % self.prog()); return False
#########################################################################################################################
    ### <4> ### SynBad Tidy Methods                                                                                     #
#########################################################################################################################
    def tidyTables(self):   ### SynBad tidying of gaps with clear fixes
        '''
        This method works through a table and identifies regions between pairs of gaps that might possibly be able to be
        fixed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if not self.dev():
            #    self.printLog('#TIDY','Table tidy functions dev only (dev=T)')
            #    return True
            self.headLog('SYNBAD GAP & SCAFFOLD TIDY PROCESSING',line='=')

            #i# Notes for future updates ...
            #?# Add a gaphide=LIST mode to "hide" gaps with a certain final classification.
            #?# Add a gapshow=FASFILE mode to resinstate hidden gaps in a given file.
            #?# Add a gapignore=LIST option to ignore certain gaps during compression? (Or just encourage hide -> re-run -> show?)
            #i# NOTE: Don't actually want to hide these gaps, I think. It is possible that maxsynskip settings could
            #i# result in a Syn gap that could still contain another mis-placed contig.

            ### ~ [1] Real translocations -> Divergent/Long Synteny ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Any Tran/Brk gaps that link regions of two scaffolds where scaffold 1 is unambiguously linked up/downstream
            #i# to scaffold 2, should be reclassified as "Div" (Divergent). Make sure the Hit does not have any gaps here too?
            #i# Probably also want to have some kind of proximity filter. Could this be the same as the merging code,
            #i# except re-labelling rather than merging due to the presence of the gap(s)?
            #i# Any Brk gaps that are flanked by collinear hits > maxsynpan bp apart will be re-classified as 'Long'
            for qh in ('Qry','Hit'):
                locdb = self.dbTable(qh)
                gapdb = self.dbTable(qh,'gap')
                self.db().addEmptyTable('{0}.corrections'.format(qh.lower()),['seqname','start','end','flank1','flank2','edit','details','synbad'],['seqname','start','end','edit'])
                # cordb = self.dbTable(qh,'corrections')
                skipx = self.getInt('MaxSynSkip') + 2
                self.headLog('SynBad {0} translocation/breakpoint assessment'.format(qh),line='~')
                self.synBadDivergence(qh,locdb,gapdb,skipx,quick=True)

            ### ~ [2] HiC support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Updates based on best HiC support (BestFlank matching GapFlank)
            for qh in ('Qry','Hit'):
                locdb = self.db(qh.lower())
                gapdb = self.dbTable(qh,'gap')
                hx = 0; dx = 0; px = 0
                self.headLog('SynBad {0} HiC pair assessment'.format(qh),line='~')
                #!# Should there be a toggle to allow reclassification of partial HiC support? Would this indicate a duplication problem?!
                for gentry in gapdb.entries():
                    if gentry['BestFlank3'] == gentry['GapFlank3'] and gentry['BestFlank5'] == gentry['GapFlank5']:
                        if gentry['SynBad'] in hicgaps:
                            gentry['SynBad'] = 'HiC'; hx += 1
                            locdb.data((gentry['seqname'],gentry['start'],gentry['end']))['QryGap'] = 'HiC'
                        elif gentry['SynBad'] in ['Dup']:
                            gentry['SynBad'] = 'DupHiC'; dx += 1
                            locdb.data((gentry['seqname'],gentry['start'],gentry['end']))['QryGap'] = 'DupHiC'
                    elif gentry['SynBad'] not in hicgaps + ['Dup']:
                        px += 1
                    elif gentry['BestFlank3'] == gentry['GapFlank3'] or gentry['BestFlank5'] == gentry['GapFlank5']:
                        self.printLog('#HIC','{0} Gap {1}::{2} has partial best-HiC support'.format(gentry['SynBad'],gentry['GapFlank5'],gentry['GapFlank3']))
                self.printLog('#HIC','{0} {1} gaps with best HiC-scoring flank pairs reclassified to "HiC"'.format(hx,qh))
                self.printLog('#HIC','{0} {1} gaps with best HiC-scoring flank pairs reclassified to "DupHiC"'.format(dx,qh))
                self.printLog('#HIC','{0} {1} gaps with one-directional best HiC-scoring flank pair'.format(px,qh))

            ### ~ [3] Inversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Inversions consist of a region flanked by two Inv QryGaps. In the simplest case, there are no intervening
            #i# gaps, but the following gaps can also be skipped: Syn, Span, Dup. (Aln should not appear!)
            for qh in ('Qry','Hit'):
                locdb = self.db(qh.lower())
                gapdb = self.dbTable(qh,'gap')
                #skipx = self.getInt('MaxSynSkip') + 2
                self.headLog('SynBad {0} inversion assessment'.format(qh),line='~')
                self.synBadInversions(qh,locdb,gapdb,0,quick=True)
                self.synBadInvertedDuplications(qh,locdb,gapdb,0,quick=True)

            ### ~ [4] Update hitgaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ('Qry','Hit'):
                self.updateHitGHaps(qh)

            #># Generate a new fragments table that skips through Syn, Aln, Span, InvFix, InvDupFix, Long and Div gaps. Sort by size?
            # and look for out of place neighbours. Build a collinear backbone in size order? Base this on the ends.
            # Identity all that don't fit. Sub classify as cis and trans. Look at whether cis can be moved into any gaps.
            # Need to identify cis and trans gaps too - include the Div gaps as cis, ignoring the non-cis region.
            # This latter bit is basically looking at the pairs tables but including the gaps in every pair?
            # Need to rename the Div gaps so they have the cis scaffolds in their names.
            # Can then add cis gaps to the pairs table and work off that.
            for qh in ('Qry','Hit'):
                self.headLog('SynBad {0} synteny blocks and cis/trans assignment'.format(qh),line='~')
                self.synBadSyntenyBlocks(qh)


            ### ~ [!] Duplications ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# These need the addition of read depth data. Here, the duplicated regions should be identified and then
            #i# analysed for the relative probability of being <=0.5N versus >=1N. If the region between two gaps,
            #i# including a Dup gap, is <=0.5N then it should be removed: will need to be careful not to remove two!
            #i# If a region is spanned by two Dup gaps then it is the best candidate: focus on these first!
            #?# NOTE: I think this needs to be a different tool. SynBad should output plots that can highlight issues.

            ### ~ [!] Extract debris ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# If the BestFlank3 is downstream and the GapFlank3 has no HiC support, and the BestFlank5 gap also has
            #i# no HiC support, then remove the entire section between GapFlank5 and BestFlank5 and move to "debris", i.e.
            #i# extract into its own sequence.

            ### ~ [!] Relocate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Look for a chunk between flanks where the ends (1) do not have BestFlank support, and (2) are the
            #i# two BestFlank entries for the same gap -> lift out that chunk and insert into the gap. (May need to invert.)

            ### ~ [!] Translocate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Look for a pair of gaps without BestFlank support
            # reciprocal swap where the BestFlank5 of Gap 1 is the GapFlank5 of Gap2 and vice versa.
            #i# Then check that the BestFlank3 for both is also not the current GapFlank3. (Would this indicate a duplication?)
            #i# If the BestFlank is not GapFlank5, identify the BestFlank5 for the current GapFlank3 of the gap that
            #i# has the focal BestFlank5 as its GapFlank5. Then scan downstream for


            #i# ~ [!] Breakpoints ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# When a Brk gap is encountered, look for something that would fit in the gap. If something is found,
            #i# and has Tran/Brk flanks, replace them with 'Mis' ("Misassembly").


            ### ~ [!] Translocations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Find regions that are flanked by two Trans gaps with the same Hit. Then find another Trans gap that is
            #i# adjacent and close enough to the flanked region to move.

            #!# Will want to differentiate edits between transpositions and single break re-join swaps.

            #!# Identify a gap that has a different BestJoin to the current GapFlank
            # -> Identify the current contig that is on that side
            # -> Identify the better contig to be on that side => Check the score for reciprocal improvement
            # -> Q. Should we check the other end of each contig? I think not - let's hope the single rearrangements solve it!


            #i# ~ [!] Dodgy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Identify dodgy regions that have non-syntenic gaps and duplications etc. without a clear fix.
            #!# After other fixes are applies, identify the Dup regions and output to a file. (Maybe with the gaps
            #!# themselves?) Then run rje_lrbridge.regionSpan(self,regfile). If no edits made, can just use the BAM
            #!# file already made of mapped reads. Otherwise, will need to re-map read data to corrected output.
            #!# -> Use the assembly map to identify the regions?

            return True
        except: self.errorLog('%s.tidyTables error' % self.prog()); return False
#########################################################################################################################
    def updateHitGHaps(self,qh): ### Uses alt gapdb gap ratings to update qh locdb.
        '''
        Uses alt gapdb gap ratings to update qh locdb.
        >> qh:str = Qry or Hit
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Updating SynBad {0} hit gaps'.format(qh),line='~')
            locdb = self.db(qh.lower())
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            gapdb = self.dbTable(qh,'gap')
            hitsort = []    # (hitseq,hitpos,start/end/gap,entry)
            for field in ['Hit3','Hit5','Qry3','Qry5']:
                if field not in locdb.fields(): locdb.addField(field,after='HitGap',evalue='.')
            q5 = '{0}5'.format(qh)
            q3 = '{0}3'.format(qh)
            qgap = '{0}Gap'.format(qh)
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            ### ~ [1] Sort QryGaps and build hitsort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            entries = locdb.entries(sorted=True)
            etot = len(entries) * 4 + gapdb.entryNum() * 3
            ex = 100.0
            entry = entries.pop(0)
            hitsort.append((entry[alt],entry[hstart],'start',entry))
            hitsort.append((entry[alt],entry[hend],'end',entry))
            while entries:
                self.progLog('\r#SYNGAP','Updating {0} {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                prev = entry
                entry = entries.pop(0)
                entry[hgap] = ''
                entry[h5] = entry[h3] = '-'
                if self.isGap(entry[qgap]): continue
                entry[qgap] = ''
                if self.isGap(prev[qgap]):
                    entry[qgap] += '<'; entry[q5] = prev[qgap]
                if entries and self.isGap(entries[0][qgap]):
                    entry[qgap] += '>'; entry[q3] = entries[0][qgap]
                if not entry[qgap]: entry[qgap] = '.'
                hitsort.append((entry[alt],entry[hstart],'start',entry))
                hitsort.append((entry[alt],entry[hend],'end',entry))
            ## ~ [1a] Add hit gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gentry in gapdb.entries():
                self.progLog('\r#UPDATE','Updating {0} flanking {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                hitsort.append((gentry['seqname'],gentry['start'],'gapstart',gentry['SynBad']))
                hitsort.append((gentry['seqname'],gentry['end'],'gapend',gentry['SynBad']))
            hitsort.sort()
            ### ~ [2] Update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(hitsort)):
                self.progLog('\r#UPDATE','Updating {0} flanking {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                #i# Find gap and look each side
                if hitsort[i][2] == 'gapstart':
                    if i == 0: continue
                    j = i-1
                    if hitsort[j][2].startswith('gap'): continue
                    jend = hitsort[j][1]
                    while j >= 0 and hitsort[j][1] == jend:
                        if hitsort[j][2] == 'end':
                            hitsort[j][3][hgap] += '>'
                            hitsort[j][3][h3] = hitsort[i][3]
                        j -= 1
                if hitsort[i][2] == 'gapend':
                    j = i+1
                    if j >= len(hitsort): continue
                    if hitsort[j][2].startswith('gap'): continue
                    jstart = hitsort[j][1]
                    while j < len(hitsort) and hitsort[j][1] == jstart:
                        if hitsort[j][2] == 'start':
                            hitsort[j][3][hgap] += '<'
                            hitsort[j][3][h5] = hitsort[i][3]
                        j += 1
            for entry in locdb.entries():
                self.progLog('\r#UPDATE','Updating {0} flanking {1}: {2:.1f}%%'.format(qh,hgap,ex/etot)); ex += 100.0
                if not entry[hgap]: entry[hgap] = '.'
                if entry[hgap] == '><': entry[hgap] = '<>'
            self.printLog('\r#UPDATE','Updating {0} flanking {1}, {2} and {3} complete.'.format(qh,hgap,h5,h3))

            return True
        except: self.errorLog('%s.updateHitGHaps error' % self.prog()); return False
#########################################################################################################################
    def isGap(self,gapstr): ### Returns whether gapstr is a gap entry.
        if type(gapstr) == int: return False
        if ('>' in gapstr or '<' in gapstr or gapstr == '.'): return False
        return True
#########################################################################################################################
    def synBadSyntenyBlocks(self,qh,quick=False,skiptypes=None):   ### SynBad collinear chunks and Cis/Trans assessment
        '''
        SynBad collinear chunks and Cis/Trans assessment.
        #># Generate a new fragments table that skips through Syn, Aln, Span, Long and Div gaps. Sort by size?
        # and look for out of place neighbours. Build a collinear backbone in size order? Base this on the ends.
        # Identity all that don't fit. Sub classify as cis and trans. Look at whether cis can be moved into any gaps.
        # Need to identify cis and trans gaps too - include the Div gaps as cis, ignoring the non-cis region.
        # This latter bit is basically looking at the pairs tables but including the gaps in every pair?
        # Need to rename the Div gaps so they have the cis scaffolds in their names.
        # Can then add cis gaps to the pairs table and work off that.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        >> skiptypes:list [] = List of gap types to skip when making chunks. (None=Default)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if skiptypes == None:
                skiptypes = skipgaps
            locdb = self.db(qh.lower())
            gapdb = self.dbTable(qh,'gap')
            #if 'SynType' not in gapdb.fields(): gapdb.addField('SynType')
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            ## ~ [0a] Generate synteny blocks table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            syntable = '{0}.blocks'.format(qh.lower())
            syndb = self.db().copyTable(locdb,syntable,replace=True,add=True)
            #syndb.addField('GapNum',evalue=0)
            syndb.addField('SynType',evalue="")
            syndb.dropFields(['Alt','Trans','Inv','Dup',hgap,'Qry5','Qry3','Hit5','Hit3'])
            entries = syndb.entries(sorted=True)

            ### ~ [1] Generate synteny blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while entries:
                #i# Start new block
                block = entries.pop(0)
                #i# Check for intervening gap to skip
                if self.isGap(block[qgap]):
                    block['SynType'] = 'Gap'
                    continue
                block[qgap] = 0
                #i# Find end of block
                blockend = block
                prev = block
                while entries and (entries[0][qh] == block[qh]) and ('>' in entries[0][qgap] or '<' in entries[0][qgap] or entries[0][qgap] == '.' or entries[0][qgap] in skiptypes):
                    if entries[0][qgap] in skiptypes: block[qgap] += 1
                    blockend = entries.pop(0)
                    #i# Will probably want to fix or drop these at some point
                    for field in ['Length','Identity','Non']:
                        block[field] += blockend[field]
                    block['Non'] += max(0,blockend[qstart]-prev[qend]-0)
                    prev = blockend
                    syndb.dropEntry(blockend)
                #i# Update entry
                if blockend:
                    block[qend] = blockend[qend]
                    bstart = block[hstart]
                    bend = blockend[hend]
                    if block['Strand'] == '-': bstart = block[hend]
                    if blockend['Strand'] == '-': bend = blockend[hstart]
                    block[hstart] = bstart
                    block[hend] = bend
                    block['Strand'] = block['Strand'] + blockend['Strand']
                    if block[alt] != blockend[alt]:
                        block['SynType'] = 'Trans'
                        block[alt] = '{0}::{1}'.format(block[alt],blockend[alt])
                    else: block['SynType'] = 'Cis'
            #syndb.dropField(qgap)
            syndb.indexReport('SynType')
            syndb.saveToFile()

            return True
        except: self.errorLog('%s.synBadChunkTable error' % self.prog()); return False
#########################################################################################################################
    def synBadDivergence(self,qh,locdb,gapdb,maxskip=0,quick=False,gaptypes=('Brk','Tran')):   ### Main SynBad gap mapping compression for a local hits table.
        '''
        Main SynBad gap mapping compression for a local hits table.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: quick = True
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = [('NSCUCHR1.01',95640,477499),('NSCUCHR1.01',1312891,1314913),
                       ('NSCUCHR1.01',1320580,1323528),
                       ('NSCUCHR1.01',1316193,1319918)]
            debugme = []

            ### ~ [1] Cycle through entries, identifying and assessing flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            dx = 0
            while entries:
                self.progLog('\r#DIV','Assessing %s table divergence: %.2f%%' % (tname,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if lookback:
                        lookback = True
                    else:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                else:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = False # self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\nSkip={0}; Lookback={1}'.format(eskip,lookback))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.': continue #i# Don't merge gaps
                #i# - (2) Currently only want one internal gap of chosen type
                skipgap = False
                gapnum = 0
                trangap = None
                for entry in toskip:
                    if '>' in entry[qgap] or '<' in entry[qgap] or entry[qgap] == '.': continue
                    if entry[qgap] not in gaptypes: skipgap = True
                    else: gapnum += 1; trangap = entry
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                elif gapnum != 1:
                    if debug: self.bugPrint('NoFocalGap...')
                    continue
                #i# - (3) Lack of hit gaps - certain good gaps are allowed
                for entry in [merge1,merge2] + toskip:
                    #X#if entry[hgap] != '.': skipgap = True
                    if entry[h5] not in goodgaps + ['-'] and entry[h3] not in goodgaps + ['-']: skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                maxspan = False
                if qdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                    maxspan = True
                if hdist > self.getInt('MaxSynSpan'):
                    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                    maxspan = True
                #i# maxspan=True when maximum length breached: used for converting Brk to Long
                if maxspan:
                    if trangap[qgap] != 'Brk' or eskip != 1: continue

                ## ~ [1c] Re-classify gap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Next, work through, classify and merge each intermediate entry and merge2
                try:
                    gkey = locdb.makeKey(trangap)
                    newgap = 'Div'
                    if maxspan: newgap = 'Long'
                    self.printLog('\r#DIV','Reclassifying {0} {1}-{2} {3} {4} -> {5}'.format(gkey[0],gkey[1],gkey[2],trangap[alt],trangap[qgap],newgap))
                    gentry = gapdb.data(gkey)
                    gentry['SynBad'] = newgap
                    trangap[qgap] = newgap
                    if debug:
                        self.bugPrint('%s' % (locdb.entrySummary(trangap,collapse=True)))
                        self.debug('%s' % (gapdb.entrySummary(gentry,collapse=True)))
                    dx += 1
                except:
                    self.errorLog('%s.synBadDivergence error' % self.prog())
                    self.warnLog('Failed to update gap class')
            self.printLog('\r#DIV','Assessed %s table divergence: %s gaps reclassified' % (tname,rje.iStr(dx)))

            return dx
        except: self.errorLog('%s.synBadDivergence error' % self.prog()); raise
#########################################################################################################################
    def synBadInversions(self,qh,locdb,gapdb,maxskip=0,quick=False):   ### SynBad inversion correction.
        '''
        SynBad inversion correction.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            tdb = self.dbTable(qh,'corrections')
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: maxskip = etot
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = [('NSCUCHR1.01',11501373,11553304),('NSCUCHR1.01',58321100,58755424)]

            ### ~ [1] Cycle through entries, identifying and assessing flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            dx = 0
            focus = None
            while entries:
                if focus and locdb.makeKey(focus) in debugme: self.deBug('?')
                self.progLog('\r#DIV','Assessing simple %s table inversions: %.2f%%' % (tname,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if lookback:
                        lookback = True
                    else:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                else:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\n\nSkip={0}; Lookback={1}; Focus={2}'.format(eskip,lookback,locdb.makeKey(focus)))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                # if debug:
                #     for entry in [merge1] + toskip:
                #         self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                #     self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (0) Inversions
                if len(toskip) < 3: continue
                if toskip[0][qgap] != 'Inv':
                    eskip = maxskip
                    continue
                if toskip[-1][qgap] != 'Inv': continue
                if toskip[1][alt] != merge1[alt]:
                    eskip = maxskip
                    continue
                if toskip[-2][alt] != merge1[alt]:
                    continue
                eskip = maxskip
                #if not debug: self.bugPrint('\n\nSkip={0}; Lookback={1}; Focus={2}'.format(eskip,lookback,locdb.makeKey(focus)))
                #debug = self.debugging()
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.bugPrint('%s' % (locdb.entrySummary(merge2,collapse=True)))
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.':
                    continue #i# Don't merge gaps
                #i# - (2) Will allow skipping of "good" gaps, including Dup gaps.
                skipgap = False
                trangap = None
                for entry in toskip[1:-1]:
                    #X#if '>' in entry[qgap] or '<' in entry[qgap] or entry[qgap] == '.': continue
                    if self.isGap(entry[qgap]) and entry[qgap] not in skipgaps:
                        skipgap = True
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                #i# - (3) Lack of hit gaps - certain good gaps are allowed
                for entry in [merge1,merge2] + toskip:
                    #X#if entry[hgap] != '.': skipgap = True
                    if entry[h5] not in goodgaps + ['-'] and entry[h3] not in goodgaps + ['-']: skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                #if qdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                #    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                #if hdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                #    continue
                #i# ~ (5) Positions of inverted region
                (flankstart, flankend) = (merge1[hend],merge2[hstart])
                (invstart, invend) = (toskip[-2][hstart],toskip[1][hend])
                if bwd:
                    (flankstart, flankend) = (merge2[hstart],merge1[hend])
                    (invstart, invend) = (toskip[1][hstart],toskip[-2][hend])
                #self.deBug('\n{0} -| {1} -> {2} |- {3} ?'.format(flankstart,invstart,invend,flankend))
                if invstart - flankstart + self.getInt('MaxOverlap') < 0:
                    if debug: self.bugPrint('ToSkip start too small')
                    continue
                if flankend - invend + self.getInt('MaxOverlap') < 0:
                    if debug: self.bugPrint('ToSkip end too big')
                    continue

                ## ~ [1c] Invert sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Next, work through, classify and merge each intermediate entry and merge2
                for invgap in (toskip[0],toskip[-1]):
                    gkey = locdb.makeKey(invgap)
                    self.printLog('\r#INV','Reclassifying {0} {1}-{2} {3} {4} -> InvFix'.format(gkey[0],gkey[1],gkey[2],invgap[alt],invgap[qgap]))
                    gentry = gapdb.data(gkey)
                    gentry['SynBad'] = 'InvFix'
                    invgap[qgap] = 'InvFix'
                    if debug:
                        self.bugPrint('%s' % (locdb.entrySummary(invgap,collapse=True)))
                        self.debug('%s' % (gapdb.entrySummary(gentry,collapse=True)))
                    dx += 1
                #i# Fix intervening sequences: invert and reverse order
                # i1 = merge1 -> +1 = Inv1
                # i2 = merge2 -> -1 = Inv2
                #!# No longer doing this -> might cause issues. Re-run on the edited fasta files instead.
                invstart = toskip[0][qend] + 1
                invend = toskip[-1][qstart] - 1
                #i# Make sure that this inversion is documented -> will need to check/execute later (or reverse)
                self.printLog('#FIX','Invert {0}:{1}-{2}'.format(merge1[qh],invstart,invend))
                # tdb.addEntry({'seqname':merge1[qh],'start':invstart,'end':invend,'edit':'invert','details':'Inversion in-place.'})
                flank1 = gapdb.data((toskip[0][qh],toskip[0][qstart],toskip[0][qend]))['GapFlank3']
                flank2 = gapdb.data((toskip[-1][qh],toskip[-1][qstart],toskip[-1][qend]))['GapFlank5']
                tdb.addEntry({'seqname':merge1[qh],'start':invstart,'end':invend,'flank1':flank1,'flank2':flank2,'edit':'invert','details':'Inversion in-place.'})
            locdb.remakeKeys()
            self.printLog('\r#INV','Assessed simple %s table inversions: %s gaps reclassified' % (tname,rje.iStr(dx)))
            if dx: self.warnLog('NOTE: Fixed inversions will not carry over into %s table. (Fix assemblies and re-run SynBad.)' % alt)
            return dx
        except: self.errorLog('%s.synBadInversions error' % self.prog()); raise
#########################################################################################################################
    def synBadInvertedDuplications(self,qh,locdb,gapdb,maxskip=0,quick=False):   ### Extended SynBad inversion correction.
        '''
        Extended SynBad inversion correction.
        >> qh:str = Whether the Query focus is the 'Qry' or 'Hit'
        >> locdb:Table = Local hits table to compress
        >> eksip:int [0] = Number of intervening entries to skip and remove/incorporate. (Exact)
        >> quick:bool [False] = Whether to use quick (position sorted) processing order or slow (length sorted)
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alt = {'Qry':'Hit','Hit':'Qry'}[qh]
            tname = locdb.name()
            tdb = self.dbTable(qh,'corrections')
            #i# qstart and qend will always be sorted in order
            qstart = '{0}Start'.format(qh)
            qend = '{0}End'.format(qh)
            qgap = '{0}Gap'.format(qh)
            #i# hstart and hend will always be sorted in reverse order for a pair of -ve strand entries
            hstart = '{0}Start'.format(alt)
            hend = '{0}End'.format(alt)
            hgap = '{0}Gap'.format(alt)
            h5 = '{0}5'.format(alt)
            h3 = '{0}3'.format(alt)
            #i# Setup sorted entries to process
            entries = locdb.entries(sorted=True)
            etot = len(entries); ex = 0.0
            sorted = entries[0:]    # List of entries to work through as the focal entry (merge1 or merge2)
            if maxskip == 0: maxskip = etot
            if not quick: sorted = locdb.sortedEntries('Length',reverse=True)

            debugme = [('NSCUCHR1.01',11501373,11553304),('NSCUCHR1.01',58321100,58755424)]

            ### ~ [1] Cycle through entries, identifying and assessing flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lookback = False    # Whether to look back in this round, i.e. the focal entry is merge2
            eskip = 0
            dx = 0
            while entries:
                self.progLog('\r#DIV','Assessing complex %s table inversions: %.2f%%' % (tname,ex/etot))
                #i# In quick mode, will either merge and keep focus, or not merge and move on to next entry
                #i# In slow mode, we only worry about merge after lookback has been activate, i.e. if lookback=True
                #i# If lookback is false and no merge, then switch lookback to True and keep the focus
                if not quick:
                    if lookback:
                        lookback = True
                    else:
                        lookback = False
                        if eskip < maxskip: eskip += 1
                        else:
                            eskip = 0
                            entries.remove(sorted.pop(0)); ex += 100.0
                            if not entries: break
                else:
                    if eskip < maxskip:
                        eskip += 1
                    else:
                        eskip = 0
                        entries.remove(sorted.pop(0)); ex += 100.0
                        if not entries: break

                ## ~ [1a] Setup entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                focus = sorted[0]
                debug = self.debugging() and locdb.makeKey(focus) in debugme
                if debug: self.bugPrint('\nSkip={0}; Lookback={1}'.format(eskip,lookback))
                if focus not in entries: continue
                i1 = entries.index(focus)
                i2 = i1 + eskip + 1
                if lookback: i2 = i1; i1 -= (eskip + 1)
                if i1 < 0 or i2 >= len(entries): continue
                merge1 = entries[i1]
                merge2 = entries[i2]
                toskip = entries[i1+1:i2]
                # if debug:
                #     for entry in [merge1] + toskip:
                #         self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                #     self.debug('%s' % (locdb.entrySummary(merge2,collapse=True)))

                ## ~ [1b] Identity whether flanks can be merged ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# - (0) Inversions
                if len(toskip) < 3: continue
                if toskip[0][qgap] != 'Inv':
                    eskip = maxskip
                    continue
                if toskip[-1][qgap] != 'Inv': continue
                if toskip[1][alt] != merge1[alt]:
                    eskip = maxskip
                    continue
                if toskip[-2][alt] != merge1[alt]:
                    continue
                eskip = maxskip
                debug = self.debugging()
                if debug:
                    for entry in [merge1] + toskip:
                        self.bugPrint('%s' % (locdb.entrySummary(entry,collapse=True)))
                    self.bugPrint('%s' % (locdb.entrySummary(merge2,collapse=True)))
                #i# - (1) Colinearity
                matched = True
                for field in ['Qry','Hit','Strand']:
                    matched = matched and merge1[field] == merge2[field]
                if not matched: continue
                if merge1['Strand'] == '.':
                    continue #i# Don't merge gaps
                #i# - (2) Will allow skipping of "good" gaps, including Dup gaps.
                skipgap = False
                trangap = None
                for entry in toskip[1:-1]:
                    #X#if '>' in entry[qgap] or '<' in entry[qgap] or entry[qgap] == '.': continue
                    if self.isGap(entry[qgap]) and entry[qgap] not in skipgaps:
                        skipgap = True
                if skipgap:
                    if debug: self.bugPrint('SkipGap...')
                    continue
                #i# - (3) Lack of hit gaps - certain good gaps are allowed
                for entry in [merge1,merge2] + toskip:
                    #X#if entry[hgap] != '.': skipgap = True
                    if entry[h5] not in goodgaps + ['-'] and entry[h3] not in goodgaps + ['-']: skipgap = True
                if skipgap:
                    if debug: self.bugPrint('HitGap...')
                    continue
                #i# - (4) Proximity
                bwd = merge1['Strand'] == '-'
                qdist = merge2[qstart] - merge1[qend]
                hdist = merge2[hstart] - merge1[hend]
                if bwd: hdist = merge1[hstart] - merge2[hend]
                if debug: self.bugPrint('QryDist={0}; HitDist={1}; (-ve:{2})'.format(qdist,hdist,bwd))
                if self.getInt('MaxOverlap') + qdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), qdist))
                    continue
                #if qdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), qdist))
                #    continue
                if self.getInt('MaxOverlap') + hdist < 0:
                    if debug: self.bugPrint('{0} + {1} < 0'.format(self.getInt('MaxOverlap'), hdist))
                    continue
                #if hdist > self.getInt('MaxSynSpan'):
                #    if debug: self.bugPrint('{0} < {1}'.format(self.getInt('MaxSynSpan'), hdist))
                #    continue
                #i# ~ (5) Positions of inverted region
                (flankstart, flankend) = (merge1[hend],merge2[hstart])
                (invstart, invend) = (toskip[-2][hstart],toskip[1][hend])
                if bwd:
                    (flankstart, flankend) = (merge2[hstart],merge1[hend])
                    (invstart, invend) = (toskip[1][hstart],toskip[-2][hend])
                #self.deBug('\n{0} -| {1} -> {2} |- {3} ?'.format(flankstart,invstart,invend,flankend))

                ## ~ [1c] Invert sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                invert = False
                #i# Next, work through, classify and merge each intermediate entry and merge2
                if invstart - flankstart + self.getInt('MaxOverlap') >= 0:
                    toskip[0][qgap] = 'InvFix'; invert = True
                elif invstart < flankstart < invend:
                    toskip[0][qgap] = 'InvDupFix'; invert = True
                else:
                    toskip[0][qgap] = 'InvBrk'
                if flankend - invend + self.getInt('MaxOverlap') >= 0:
                    toskip[-1][qgap] = 'InvFix'; invert = True
                elif invstart < flankend < invend:
                    toskip[-1][qgap] = 'InvDupFix'; invert = True
                else:
                    toskip[-1][qgap] = 'InvBrk'

                for invgap in (toskip[0],toskip[-1]):
                    gkey = locdb.makeKey(invgap)
                    gentry = gapdb.data(gkey)
                    gentry['SynBad'] = invgap[qgap]
                    self.printLog('\r#INV','Reclassifying {0} {1}-{2} {3} Inv -> {4}'.format(gkey[0],gkey[1],gkey[2],invgap[alt],invgap[qgap]))
                    if debug:
                        self.bugPrint('%s' % (locdb.entrySummary(invgap,collapse=True)))
                        self.debug('%s' % (gapdb.entrySummary(gentry,collapse=True)))
                    dx += 1


                #i# Fix intervening sequences: invert and reverse order
                # i1 = merge1 -> +1 = Inv1
                # i2 = merge2 -> -1 = Inv2
                if invert:
                    invstart = toskip[0][qend] + 1
                    invend = toskip[-1][qstart] - 1
                    #i# Make sure that this inversion is documented -> will need to check/execute later (or reverse)
                    self.printLog('#FIX','Invert {0}:{1}-{2}'.format(merge1[qh],invstart,invend))
                    # tdb.addEntry({'seqname':merge1[qh],'start':invstart,'end':invend,'edit':'invert','details':'Inversion in-place.'})
                    flank1 = gapdb.data((toskip[0][qh],toskip[0][qstart],toskip[0][qend]))['GapFlank3']
                    flank2 = gapdb.data((toskip[-1][qh],toskip[-1][qstart],toskip[-1][qend]))['GapFlank5']
                    tdb.addEntry({'seqname':merge1[qh],'start':invstart,'end':invend,'flank1':flank1,'flank2':flank2,'edit':'invert','details':'Inversion in-place.'})
                    # tdb.addEntry({'seqname':merge1[qh],'flank1':invstart,'flank2':invend,'edit':'invert','details':'Inversion in-place.'})
            locdb.remakeKeys()
            self.printLog('\r#INV','Assessed complex %s table inversions: %s gaps reclassified' % (tname,rje.iStr(dx)))
            if dx: self.warnLog('NOTE: Fixed inversions will not carry over into %s table. (Fix assemblies and re-run SynBad.)' % alt)

            return dx
        except: self.errorLog('%s.synBadInvertedDuplications error' % self.prog()); raise
#########################################################################################################################
    ### <5> ### SynBad Output Methods                                                                                   #
#########################################################################################################################
    def saveTables(self,backup=True):   ### Saves the output tables, once all processing is complete.
        '''
        Saves the output tables, once all processing is complete.
        '''
        try:### ~ [1] Save gap and local tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SAVE SYNBAD TABLES',line='=')
            for qh in ('Qry','Hit'):
                self.dbTable(qh,'gap').saveToFile(backup=backup)
                self.dbTable(qh,'corrections').saveToFile(backup=backup)
                locdb = self.db(qh.lower())
                locdb.dropField('AlnNum')
                outfields = locdb.fields()[0:]
                try:
                    outfields.remove('{0}5'.format(qh))
                    outfields.remove('{0}3'.format(qh))
                except:
                    pass    # Will need to update pairs tables!
                locdb.saveToFile(savefields=outfields,backup=backup)

            ### ~ [2] Paired scaffold output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## ~ [2e] Top-matched pairs only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bestpair = self.getBool('BestPair')
            if not bestpair: return True
            pairs = []
            for qh in ('Qry','Hit'):
                topdb = db.copyTable(self.db(qh.lower()),'top{0}'.format(qh.lower()),replace=True,add=True)
                topdb.keepFields(['Qry','Hit','Length']+topdb.keys())
                topdb.compress(['Qry','Hit'],default='sum')
                topdb.keepFields(['Qry','Hit','Length'])
                topdb.rankFieldByIndex(qh,'Length',newfield='Rank',rev=True,absolute=True,lowest=True,unique=False,warn=True,highest=False)
                topdb.dropEntriesDirect('Rank',[1],inverse=True,log=True,force=False)
                pairs += topdb.dataKeys()
            ## ~ [2b] Output pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                topdb = db.copyTable(self.db(qh.lower()),'{0}.pairs'.format(qh.lower()),replace=True,add=True)
                if bestpair:
                    ex = 0.0; etot = topdb.entryNum()
                    for ekey in topdb.datakeys()[0:]:
                        self.progLog('\r#PAIRS','Reducing %s table to top-aligned pairs: %.2f%%' % (qh,ex/etot)); ex += 100
                        entry = topdb.data(ekey)
                        if (entry['Qry'],entry['Hit']) not in pairs: topdb.dict['Data'].pop(ekey)
                    self.printLog('\r#PAIRS','Reduced %s table to %s alignments between top sequence pairs.' % (qh,rje.iStr(topdb.entryNum())))
                else:
                    if qh == 'Qry': topdb.newKey(['Qry','Hit','QryStart','QryEnd'])
                    else: topdb.newKey(['Hit','Qry','HitStart','HitEnd'])
                topdb.saveToFile(backup=backup)

            return True
        except: self.errorLog('%s.saveTables error' % self.prog()); return False
#########################################################################################################################
    def synBadSummarise(self):  ### Summarises ratings etc.
        '''
        Summarises ratings etc.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Add a summary of each class with descriptions and [Good/Fix/Dup/Bad] followed by percentage classes
            base = {'qry':rje.baseFile(self.getStr('Genome1'), strip_path=True),
                    'hit':rje.baseFile(self.getStr('Genome2'), strip_path=True)}
            ### ~ [1] Summary reports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gapx = {'Total':0,'Good':0,'Fix':0,'Dup':0,'Bad':0}
            for qh in ['qry','hit']:
                self.headLog('SYNBAD {0} {1}'.format(qh,base[qh]),line='-')
                gapdb = self.dbTable(qh,'gap')
                for gtype in gapdb.indexKeys('SynBad'):
                    if gtype not in gapdesc: gapdesc[gtype] = 'See docs for details'
                    gclass = 'Bad'
                    if gtype in goodgaps: gclass = 'Good'
                    elif 'Fix' in gtype: gclass = 'Fix'
                    elif 'Dup' in gtype: gclass = 'Dup'
                    gx = len(gapdb.index('SynBad')[gtype])
                    self.printLog('#SYNBAD','{0} {1} ({2}): {3} [{4}]'.format(gx,gtype,qh,gapdesc[gtype],gclass))
                    gapx['Total'] += gx
                    gapx[gclass] += gx

                for gclass in ['Good','Fix','Dup','Bad']:
                    self.printLog('#SYNBAD','{0:.1f}% {1} gaps in {2} SynBad class'.format(100.0*gapx[gclass]/gapx['Total'],base[qh],gclass))

            return True
        except:
            self.errorLog('%s.synBadSummarise error' % self.prog())
            return False
#########################################################################################################################
    def corrections(self):  ### Loads corrections table, edits sequences and outputs corrected assembly
        '''
        Loads corrections table, edits sequences and outputs corrected assembly.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seqObjSetup()
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            base = {'qry':rje.baseFile(self.getStr('Genome1'),strip_path=True),'hit':rje.baseFile(self.getStr('Genome2'),strip_path=True)}
            basefile = self.baseFile(strip_path=True)
            goodedits = ['invert']
            db = self.db()

            ### ~ [1] Check edits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                if self.dbTable(qh,'corrections'): db.deleteTable('corrections')
                # cdb = db.addTable(mainkeys=['seqname','start','end','edit'],name='corrections',replace=True,ignore=[],expect=True)
                # cdb.dataFormat({'start':'int','end':'int'})
                cdb = db.addTable(mainkeys=['seqname','start','end','edit'],name='{0}.corrections'.format(qh),replace=True,ignore=[],expect=True)
                prev = None
                for entry in cdb.entries(sorted=True):
                    if prev and prev['end'] >= entry['start'] and prev['seqname'] == entry['seqname']:
                        self.warnLog('Dropped correction {0} {1}-{2} due to overlapping edit region.'.format(entry['seqname'],entry['start'],entry['end']))
                        cdb.dropEntry(entry)
                        continue
                    prev = entry
                cdb.indexReport('seqname')
                cdb.indexReport('edit')
                for etype in cdb.index('edit'):
                    if etype not in goodedits:
                        self.warnLog('Edit type "{0}" not implemented.')

            #!# Update workflow to (1) correct Assembly Map, and (2) Regenerate Assembly from contigs and map

            #!# This workflow will need editing when going beyond inversions: need to be careful of clashes.
            ### ~ [2] Perform edits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                cdb = self.dbTable(qh,'corrections')
                seqlist = seqobj[qh]
                sx = 0.0; stot = seqlist.seqNum()
                fasout = '{0}.{1}.fasta'.format(base[qh],basefile)
                #!# Consider making a new set of identifiers
                seqdict = {}    # seqname:(fullname,sequence)
                edits = {}      # seqname:[edit entries in order of edit]
                self.progLog('\r#SEQ','Setting up {0} sequences'.format(qh))
                for seq in seqlist.seqs():
                    seqname = seqlist.shortName(seq)
                    seqdict[seqname] = seqlist.getSeq(seq)
                    edits[seqname] = []
                self.printLog('\r#SEQ','Setting up {0} sequences complete.'.format(qh))
                ## ~ [2a] Inversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#INV','Processing {0} inversions'.format(qh))
                editx = 0
                for ekey in cdb.indexDataKeys('edit','invert'):
                    seqname = ekey[0]
                    if seqname not in seqdict: continue  # Qry and Hit in same file
                    i = ekey[1]
                    j = ekey[2]
                    sequence = seqdict[seqname][1]
                    invreg = rje_sequence.reverseComplement(sequence[i-1:j])
                    sequence = sequence[:i] + invreg + sequence[j:]
                    seqdict[seqname] = (seqdict[seqname][0], sequence)
                    edits[seqname].append(cdb.data(ekey))
                    editx += 1
                self.printLog('\r#INV','Processed {0} inversions: {1} edits'.format(qh,editx))

            ### ~ [3] Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                rje.backup(self,fasout)
                FASOUT = open(fasout,'w'); sx = 0.0
                for seq in seqlist.seqs():
                    self.progLog('\r#FASOUT','Corrected {0} sequence output: {1:.1f}%%'.format(qh,sx/stot)); sx += 100.0
                    seqname = seqlist.shortName(seq)
                    (newname,sequence) = seqdict[seqname]
                    newname += ' ({0} SynBad edits)'.format(len(edits[seqname]))
                    FASOUT.write('>{0}\n{1}\n'.format(newname,sequence))
                FASOUT.close()
                self.printLog('\r#FASOUT','Corrected {0} sequence output to {1}'.format(qh,fasout))

            return True
        except:
            self.errorLog('%s.corrections error' % self.prog())
            return False
#########################################################################################################################
    def fragment(self):   ### SynBad fragmentation
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seqObjSetup()
            seqobj = {}
            seqobj['qryfrag'] = self.obj['SeqList1']
            seqobj['hitfrag'] = self.obj['SeqList2']
            base1 = rje.baseFile(self.getStr('Genome1'),strip_path=True)
            base2 = rje.baseFile(self.getStr('Genome2'),strip_path=True)
            db = self.db()
            basefile = self.baseFile()
            qrygap = self.db('qry.gap')
            if not qrygap:
                qrygap = db.addTable('%s.gaps.tdt' % base1,mainkeys=['seqname','start','end'],name='qry.gap',ignore=[],expect=True)
            qrygap.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            qryfrag = db.copyTable(qrygap,'qry.frag',replace=True,add=True)
            hitgap = self.db('hit.gap')
            if not hitgap:
                hitgap = db.addTable('%s.gaps.tdt' % base2,mainkeys=['seqname','start','end'],name='hit.gap',ignore=[],expect=True)
            hitgap.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            hitfrag = db.copyTable(hitgap,'hit.frag',replace=True,add=True)
            ### ~ [1] Filter gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
            # * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
            # * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
            # * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
            # * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
            # * `Inv` = `Inversion` = Flanking hits are on alternative strands.
            # * `Tran` = `Translocation` = `SynSpan` indicates matches are on different scaffolds.
            # * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different scaffolds, 1+ of which is not a chromosome scaffold.
            # * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.
            # * `Span` = `Spanned` = Any gaps with Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
            # * `Null` = No mapping between genomes for that gap.
            badgaptypes = self.list['FragTypes']
            self.headLog('SYNBAD GAP FILTER',line='=')
            for table in (qryfrag,hitfrag):
                table.addField('fragid',evalue='')
                for entry in table.entries():
                    if entry['Span0'] < self.getInt('MinReadSpan') and entry['SynBad'] in badgaptypes:
                        entry['fragid'] = 'Filter'
                table.dropEntriesDirect('fragid',['Filter'],inverse=True)
                table.newKey(['seqname','start','end','fragid'])
                table.addField('syn5',evalue='')
                table.addField('syn3',evalue='')
            ### ~ [2] Fragment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('FRAGMENT ON GAPS',line='=')
            for frag in ('qry.frag','hit.frag'):
                table = self.db(frag)
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
                seqlist = seqobj[frag]
                prev = {'seqname':None}
                fragname = ''; fx = 1
                for ekey in table.dataKeys():
                    entry = table.data(ekey)
                    #i# Add extra 5' fragment for each sequence
                    if entry['seqname'] != prev['seqname']:
                        prev = rje.combineDict({},entry)
                        prev['start'] = 1
                        prev['end'] = entry['start'] - 1
                        fragname = entry['seqname']; fx = 1
                        if rje.matchExp('\S+_\S+__(\S+)',fragname): fragname = rje.matchExp('\S+_\S+__(\S+)',fragname)[0]
                        prev['fragid'] = '{0}.{1}'.format(fragname,fx); fx += 1
                        prev['syn5'] = 'Terminal'
                        prev['syn3'] = 'Terminal'    #i# May be over-written
                        prev = table.addEntry(prev,remake=False)
                    #i# Update details of this entry
                    prev['syn3'] = entry['SynBad']
                    prev['end'] = entry['start'] - 1
                    entry['start'] = entry['end'] + 1
                    entry['fragid'] = '{0}.{1}'.format(fragname,fx); fx += 1
                    entry['syn5'] = entry['SynBad']
                    entry['syn3'] = 'Terminal'                          #i# May be over-written
                    entry['end'] = entry['seqlen']                      #i# May be overwritten
                    prev = entry
                    #self.bugPrint(table.entrySummary(entry,collapse=True))
                    if entry['seqname'] == 'NSCUSCAFF155': self.debug('?')
                #i# Add sequences without any gaps
                table.remakeKeys()
                seqnames = rje.sortKeys(table.index('seqname'))
                self.printLog('#FRAG','{0} fragments from {1} sequences with gaps'.format(rje.iStr(table.entryNum()),rje.iLen(seqnames)))
                sx = 0
                for seq in seqlist.seqs():
                    seqname = seqlist.shortName(seq)
                    if seqname not in seqnames:
                        fragname = seqname
                        if rje.matchExp('\S+_\S+__(\S+)',fragname): fragname = rje.matchExp('\S+_\S+__(\S+)',fragname)[0]
                        fragid = '{0}.1'.format(fragname)
                        table.addEntry({'fragid':fragid,'seqname':seqname,'start':1,'seqlen':seqlist.seqLen(seq),'end':seqlist.seqLen(seq),'syn5':'Terminal','syn3':'Terminal'})
                        sx += 1
                self.printLog('#SEQ','Added {0} full-length sequences without gaps'.format(rje.iStr(sx)))
                #i# Tidy up fragment table
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
                #table.newKey(['fragid'])
                table.setFields(['fragid','seqname','start','end','seqlen','syn5','syn3'])
                table.saveToFile()
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
            ### ~ [3] Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                seqdict = seqlist.seqNameDic()
                fragfas = '{0}.{1}.fasta'.format(basefile,frag)
                rje.backup(self,fragfas)
                FRAGFAS = open(fragfas,'w'); ex = 0
                for entry in table.entrySort():
                    self.progLog('\r#FASOUT','{0} fragments output to {1}'.format(rje.iStr(ex),fragfas)); ex += 1
                    #self.bugPrint(table.entrySummary(entry,collapse=True))
                    seq = seqdict[entry['seqname']]
                    (seqname,sequence) = seqlist.getSeq(seq)
                    FRAGFAS.write('>{0}\n'.format(entry['fragid']))
                    FRAGFAS.write(sequence[entry['start']-1:entry['end']]+'\n')
                FRAGFAS.close()
                self.printLog('\r#FASOUT','{0} fragments output to {1}'.format(rje.iStr(table.entryNum()),fragfas))
            return True
        except:
            self.errorLog('%s.fragment error' % self.prog())
            return False
#########################################################################################################################
    ### <6> ### SynBad HiC and Contig Flank Methods                                                                     #
#########################################################################################################################
    def synBadHiCMap(self):   ### Main SynBad gap mapping and compression
        '''
        Main SynBad gap mapping and compression.
        >> cdb1: Table of genome 1 gaps
        >> cdb2: Table of genome 2 gaps
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Might want to split this up a bit - separate HiC bit from the rest
            self.dict['FlankMap'] = {'qry':{}, 'hit':{}}
            if not self.getBool('PureFlanks'):
                self.setBool({'PureFlanks':True})
                self.warnLog('pureflanks=F not currently supported: setting pureflanks=T')
                #!# Will just need to extract full BED file for this mode
            db = self.db()
            self.seqObjSetup()
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            hicbam = {}
            hicbam['qry'] = self.getStr('HiCBAM1')
            hicbam['hit'] = self.getStr('HiCBAM2')
            bamdir = {}
            bamdir['qry'] = self.getStr('HiCDir1')
            bamdir['hit'] = self.getStr('HiCDir2')
            base = {'qry':rje.baseFile(self.getStr('Genome1'),strip_path=True),'hit':rje.baseFile(self.getStr('Genome2'),strip_path=True)}
            flanklen = self.getInt('GapFlanks')
            if flanklen < 1:
                self.printLog('#HIC','No BED file output for HiC BAM reduction (gapflanks={0})'.format(flanklen))
                return False
            else:
                self.printLog('#HIC','Gap Flanks for BED file output: {0}'.format(rje_seqlist.dnaLen(flanklen)))
            bamout = True
            if os.popen('samtools --version').read():
                self.printLog('#SYS',' '.join(os.popen('samtools --version').read().split()))
            else:
                self.printLog('Cannot open samtools: no HiCBAM analysis')
                bamout = False

            ## ~ [0a] Check and load outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gformats = {}
            for field in ['start','end','gaplen','seqlen','gaplen','MaxFlank5','MaxFlank3','Span0','Span100','Span1000',
                          'Span5000','GapSpan','QryStart','QryEnd','HitStart','HitEnd','Length','Identity','Non',
                          'Alt','Trans','Inv','Dup']:
                gformats[field] = 'int'
            dformats = {}
            for field in ['start','end','score','ctglen','id1','id2','pairs','score']: dformats[field] = 'int'
            for field in ['wtscore']: dformats[field] = 'num'
            pformats = rje.combineDict({},dformats)
            pformats['score'] = 'num'
            wanted = ['qry.tdt','hit.tdt']
            for qh in ['qry','hit']:
                wanted.append('{0}.flanks.bed'.format(qh))
                wanted.append('{0}.flanks.fasta'.format(qh))
                for tname in ['flanks','contigs','full','gap']:
                    wanted.append('{0}.{1}.tdt'.format(qh,tname))
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad Flank and HiC Map files found: loading (force=F)')
                for qh in ['qry','hit']:
                    gdb = db.addTable(name='{0}.gap'.format(qh),mainkeys=['seqname','start','end'],expect=True,ignore=[],replace=True)
                    gdb.dataFormat(gformats)
                    bed = db.addTable(name='{0}.flanks'.format(qh),mainkeys=['seqname','start','end'],expect=True,ignore=[],replace=True)
                    bed.dataFormat(dformats)
                    cdb = db.addTable(name='{0}.contigs'.format(qh),mainkeys=['seqname','start','end'],expect=True,ignore=[],replace=True)
                    cdb.dataFormat(dformats)
                    pdb = db.addTable(name='{0}.hicpairs'.format(qh),mainkeys=['region1','region2'],expect=False,ignore=[],replace=True)
                    if pdb: pdb.dataFormat(pformats)
                return True

            ### ~ [1] Extract gap flanks to BED and fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            flankpairs = {}     # Dictionary of contig flanks pairs for internal paired read removal
            pcount = {}         # HiCPairs counter used for reloading data
            for qh in ('qry','hit'):
                gdb = self.dbTable(qh,'gap')
                for field in ['GapFlank5','GapFlank3']:
                    if field not in gdb.fields(): gdb.addField(field,evalue='')
                bedfile = '{0}.{1}.flanks.bed'.format(db.baseFile(),qh)
                bed = db.addEmptyTable('{0}.flanks'.format(qh),['seqname','start','end','name','score','strand','pair'],['seqname','start','end'],log=True)
                cdb = db.addEmptyTable('{0}.contigs'.format(qh),['seqname','start','end','ctglen','name','flank5','flank3','synbad'],['seqname','start','end'],log=True)
                entry = None
                entries = gdb.entries(sorted=True)
                while entries:
                    prev = entry
                    entry = entries.pop(0)
                    #self.bugPrint('\nGap :: %s' % (gdb.entrySummary(entry,collapse=True)))
                    # 5' flank
                    i = entry['start'] - flanklen
                    if prev and prev['seqname'] == entry['seqname'] and self.getBool('PureFlanks'): i = max(prev['end'] + 1,i)
                    else: i = max(1,i)
                    j = entry['start'] - 1
                    score = 0
                    #if 'Span0' in gdb.fields(): score = entry['Span0']
                    bentry = bed.addEntry({'seqname':entry['seqname'],'start':i,'end':j,'name':'{0}.{1}-{2}'.format(entry['seqname'],i,j),'score':score,'strand':'+'},warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    entry['GapFlank5'] = '{0}.{1}-{2}'.format(entry['seqname'],i,j)
                    # 3' flank
                    i = entry['end'] + 1
                    j = min(entry['seqlen'],entry['end'] + flanklen)
                    if entries and entries[0]['seqname'] == entry['seqname']:
                        j = min(entries[0]['start']-1,j)
                    bentry = bed.addEntry({'seqname':entry['seqname'],'start':i,'end':j,'name':'{0}.{1}-{2}'.format(entry['seqname'],i,j),'score':score,'strand':'+'},warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    entry['GapFlank3'] = '{0}.{1}-{2}'.format(entry['seqname'],i,j)
                    #i# Extra end of sequence contig
                    if prev and entry['seqname'] != prev['seqname']:
                        i = max(prev['end']+1,prev['seqlen']-flanklen+1)
                        j = prev['seqlen']
                        bed.addEntry({'seqname':prev['seqname'],'start':i,'end':j,
                                      'name':'{0}.{1}-{2}'.format(prev['seqname'],i,j),'score':0,'strand':'+'},warn=False)
                        centry = {'seqname':prev['seqname'],'start':prev['end']+1,'end':prev['seqlen'],
                                  'flank5':prev['GapFlank3'],'flank3':'{0}.{1}-{2}'.format(prev['seqname'],i,prev['seqlen']),
                                  'synbad':'{0}-End'.format(prev['SynBad'])}
                        centry['name'] = '{0}.{1}-{2}'.format(centry['seqname'],centry['start'],centry['end'])
                        centry['ctglen'] = centry['end'] - centry['start'] + 1
                        cdb.addEntry(centry)
                        flankpairs[centry['flank5']] = centry['flank3']
                        flankpairs[centry['flank3']] = centry['flank5']
                        #self.bugPrint('-> Extra end contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                    #i# Start of sequence contig
                    if not prev or prev['seqname'] != entry['seqname']:
                        i = 1
                        j = min(entry['start'] - 1,flanklen)
                        bentry = bed.addEntry({'seqname':entry['seqname'],'start':i,'end':j,
                                      'name':'{0}.{1}-{2}'.format(entry['seqname'],i,j),'score':0,'strand':'+'},warn=False)
                        #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                        centry = {'seqname':entry['seqname'],'start':1,'end':entry['start'] - 1,
                                  'flank5':'{0}.{1}-{2}'.format(entry['seqname'],i,j),'flank3':entry['GapFlank5'],
                                  'synbad':'End-{0}'.format(entry['SynBad'])}
                        centry['name'] = '{0}.{1}-{2}'.format(centry['seqname'],centry['start'],centry['end'])
                        centry['ctglen'] = centry['end'] - centry['start'] + 1
                        cdb.addEntry(centry)
                        #self.bugPrint('-> Start contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                    #i# Middle of sequence contig
                    else:
                        centry = {'seqname':entry['seqname'],'start':prev['end']+1,'end':entry['start'] - 1,
                                  'flank5':prev['GapFlank3'],'flank3':entry['GapFlank5'],
                                  'synbad':'{0}-{1}'.format(prev['SynBad'],entry['SynBad'])}
                        centry['name'] = '{0}.{1}-{2}'.format(centry['seqname'],centry['start'],centry['end'])
                        centry['ctglen'] = centry['end'] - centry['start'] + 1
                        cdb.addEntry(centry)
                        #self.bugPrint('-> Mid contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                    flankpairs[centry['flank5']] = centry['flank3']
                    flankpairs[centry['flank3']] = centry['flank5']
                    #self.bugPrint('Gap = {0} => {1} flankpairs'.format(gdb.entrySummary(entry,collapse=True),len(flankpairs)))
                prev = entry
                if prev:
                    i = max(prev['end']+1,prev['seqlen']-flanklen+1)
                    j = prev['seqlen']
                    bentry = bed.addEntry({'seqname':prev['seqname'],'start':i,'end':j,
                                  'name':'{0}.{1}-{2}'.format(prev['seqname'],i,j),'score':0,'strand':'+'},warn=False)
                    #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                    centry = {'seqname':prev['seqname'],'start':prev['end']+1,'end':prev['seqlen'],
                              'flank5':prev['GapFlank3'],'flank3':'{0}.{1}-{2}'.format(prev['seqname'],i,prev['seqlen']),
                              'synbad':'{0}-End'.format(prev['SynBad'])}
                    centry['name'] = '{0}.{1}-{2}'.format(centry['seqname'],centry['start'],centry['end'])
                    centry['ctglen'] = centry['end'] - centry['start'] + 1
                    cdb.addEntry(centry)
                    #self.bugPrint('-> Final end contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                    flankpairs[centry['flank5']] = centry['flank3']
                    flankpairs[centry['flank3']] = centry['flank5']
                    #self.bugPrint('=> {0} flankpairs'.format(len(flankpairs)))

                #i# Add sequences without gaps to flanks and contigs
                seqdict = seqobj[qh].seqNameDic()
                sx = 0
                for seq in seqobj[qh].seqs():
                    seqname = seqobj[qh].shortName(seq)
                    if seqname not in gdb.index('seqname'):
                        seqlen = seqobj[qh].seqLen(seq)
                        bentry = {'seqname':seqname,'start':1,'end':min(flanklen,seqlen),'score':0,'strand':'+'}
                        bentry['name'] = '{0}.{1}-{2}'.format(bentry['seqname'],bentry['start'],bentry['end'])
                        bed.addEntry(bentry,warn=False)
                        #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                        centry = {'seqname':seqname,'start':1,'end':seqlen,'flank5':bentry['name'],'synbad':'End-End'}
                        bentry = {'seqname':seqname,'start':max(1,seqlen-flanklen),'end':seqlen,'score':0,'strand':'+'}
                        bentry['name'] = '{0}.{1}-{2}'.format(bentry['seqname'],bentry['start'],bentry['end'])
                        bed.addEntry(bentry,warn=False)
                        #self.bugPrint('Flank: %s' % (bed.entrySummary(bentry,collapse=True)))
                        centry['flank3'] = bentry['name']
                        centry['name'] = '{0}.{1}-{2}'.format(centry['seqname'],centry['start'],centry['end'])
                        centry['ctglen'] = centry['end'] - centry['start'] + 1
                        cdb.addEntry(centry); sx += 1
                        flankpairs[centry['flank5']] = centry['flank3']
                        flankpairs[centry['flank3']] = centry['flank5']
                        #self.bugPrint('-> Full-length contig: %s' % (cdb.entrySummary(centry,collapse=True)))
                self.printLog('#CONTIG','Added flanks for {0} sequences without gaps'.format(rje.iStr(sx)))

                #i# Save flanks to fasta file
                fasout = '{0}.{1}.flanks.fasta'.format(db.baseFile(),qh)
                if rje.exists(fasout) and not self.force():
                    self.printLog('#FLANK','{0} found (force=F)'.format(fasout))
                else:
                    FASOUT = open(fasout,'w'); sx = 0
                    for seqname in bed.index('seqname'):
                        sequence = seqobj[qh].getSeq(seqdict[seqname])[1]
                        #self.bugPrint('{0}: {1} bp'.format(seqname,len(sequence)))
                        for entry in bed.indexEntries('seqname',seqname):
                            self.progLog('\r#FLANK','{0} {1} gap flanks output to {2}'.format(rje.iStr(sx),qh,fasout)); sx += 1
                            flankseq = sequence[entry['start']-1:entry['end']]
                            if len(flankseq) != entry['end'] - entry['start'] + 1:
                                self.debug(bed.entrySummary(entry,collapse=True))
                                raise ValueError('Flanking sequence length mismatch for {0}'.format(entry['name']))
                            FASOUT.write('>{0}\n{1}\n'.format(entry['name'],flankseq))
                    FASOUT.close()
                    self.printLog('\r#FLANK','{0} {1} gap flanks output to {2}'.format(rje.iStr(bed.entryNum()),qh,fasout))
                bed.saveToFile(bedfile,delimit='\t',headers=False,savefields=['seqname','start','end','name'])

                #i# Save contigs to fasta file
                fasout = '{0}.{1}.contigs.fasta'.format(db.baseFile(),qh)
                if rje.exists(fasout) and not self.force():
                    self.printLog('#CTG','{0} found (force=F)'.format(fasout))
                else:
                    FASOUT = open(fasout,'w'); sx = 0
                    seqlist = seqobj[qh]
                    cx = 0.0; ctot = cdb.entryNum()
                    for seq in seqlist.seqs():
                        seqname = seqlist.shortName(seq)
                        if seqname in cdb.index('seqname'):
                            sequence = seqlist.getSeq(seq)[1]
                            for centry in cdb.indexEntries('seqname',seqname):
                                self.progLog('\r#CTG','Extracting contig sequences: {0:.1f}%'.format(cx/ctot)); cx += 100.0
                                flankseq = sequence[centry['start']-1:centry['end']]
                                if len(flankseq) != centry['end'] - centry['start'] + 1:
                                    self.debug(cdb.entrySummary(centry,collapse=True))
                                    raise ValueError('Contig sequence length mismatch for {0}'.format(centry['name']))
                                FASOUT.write('>{0}\n{1}\n'.format(centry['name'],flankseq)); sx += 1
                        else:
                            self.warnLog('Sequence {0} missing from contigs table'.format(seqname))
                    FASOUT.close()
                    self.printLog('\r#CTG','Extracted {0} sequences for {1} contigs'.format(rje.iStr(sx),rje.iStr(ctot)))
                    if sx != ctot: raise ValueError('Failure to generate full contigs fasta file!')

            ### ~ [X] Map old flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                mdb = self.mapFlanksToGenome(qh)
                if mdb:
                    for entry in mdb.entries():
                        self.dict['FlankMap'][qh][entry['GapFlank']] = entry['OldFlank']
            #!# Peform mdb = self.mapFlanksToGenome(qh) and use the mapping ot
                # (A) Update all the tables
                # (B) Add oldflank and oldpair to the bed table and update the name and pair fields
            #!# possibly move this all to function: mapNewFlanks(self): return

            ### ~ [2] Extract regions to BAM files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                extract = True  # Whether to extract reads from BAM file
                hicdir = True   # Whether HiC directory given to find existing files
                if not bamout: extract = False
                if hicbam[qh].lower() in ['','none']: hicbam[qh] = '{0}.HiC.bam'.format(base[qh])
                if not rje.exists(hicbam[qh]):
                    self.printLog('#HICBAM','HiC BAM file for {0} not found: {1}'.format(qh,hicbam[qh]))
                    extract = False
                if bamdir[qh].lower() in ['','none']: bamdir[qh] = '{0}.{1}flanks/'.format(db.baseFile(),qh)
                if not rje.exists(bamdir[qh]):
                    self.printLog('#HICDIR','HiC extraction directory for {0} not found: {1}'.format(qh,bamdir[qh]))
                    hicdir = False
                    if extract: rje.mkDir(self,bamdir[qh])
                if extract:
                    self.printLog('#HICBAM','HiC BAM file for {0} found: {1}'.format(qh,hicbam[qh]))
                elif hicdir and not self.force():
                    self.printLog('#HICDIR','HiC BAM file for {0} not found but HiC directory exists for reading previous results: {1}'.format(qh,bamdir[qh]))
                else:
                    continue
                self.setBool({'FullHiC':False})
                dformats = {}
                for field in ['start','end','score','ctglen','id1','id2','pairs','score']: dformats[field] = 'int'
                for field in ['wtscore']: dformats[field] = 'num'
                pformats = rje.combineDict({},dformats)
                pformats['score'] = 'num'
                pdb = None; pcount[qh] = 0
                if not self.force():
                    pdb = db.addTable(name='{0}.hicpairs'.format(qh),mainkeys=['region1','region2'],expect=False,ignore=[],replace=True)
                if pdb:
                    pdb.dataFormat(pformats)
                    pcount[qh] = pdb.entryNum()
                else:
                    pdb = db.addEmptyTable('{0}.hicpairs'.format(qh),['region1','region2','id1','id2','pairs','type','score','wtscore'],['region1','region2'])
                    uniqfile = '{0}{1}.uniq.id'.format(bamdir[qh],qh)
                    if rje.exists(uniqfile) and not self.force():
                        uniqid = string.split(open(uniqfile,'r').read())
                        self.printLog('\r#HICBAM','{0} {1} BAM flank read IDs read from {2} (force=F).'.format(rje.iLen(uniqid),qh,uniqfile))
                    elif extract:
                        bx = 0; fx = 0; tx = 0; ex = 0.0; etot = bed.entryNum()
                        readids = []
                        for entry in bed.entries():
                            self.progLog('\r#HICBAM','Extracting {0} BAM flank regions and reads: {1:.1f}%'.format(qh,ex/etot)); ex += 100.0
                            region = '{0}:{1}-{2}'.format(entry['seqname'],entry['start'],entry['end'])
                            regtmp = '%s%s.tmp' % (bamdir[qh],entry['name'])
                            sampipe = "samtools view %s %s | awk '{print $1;}' | sort | uniq -c | awk '$1 == 1' | awk '{print $2;}' | tee %s" % (hicbam[qh],region,regtmp)
                            regids = os.popen(sampipe).readlines()
                            readids += regids
                            entry['score'] = min(1000,len(regids))
                            if rje.exists(regtmp): tx += 1
                            else:
                                self.warnLog('Read ID tmp file {0} not created'.format(regtmp))
                                self.debug(sampipe)
                            if self.dev():  #?# Add as options?
                                regbam = '%s%s.bam' % (bamdir[qh],entry['name'])
                                fasout = '%s%s.reads.fasta' % (bamdir[qh],entry['name'])
                                samcmd = 'samtools view -h -F 0x104 {0} {1} | samtools view -h -b -o {2} - 2>&1'.format(hicbam[qh],region,regbam)
                                os.popen(samcmd).read()
                                if rje.exists(regbam): bx += 1
                                else: self.warnLog('BAM file {0} not created'.format(regbam))
                                bam2fas = 'samtools fasta {0} > {1} 2>&1'.format(regbam,fasout)
                                os.popen(bam2fas).read()
                                #!# Parse processed 3452 reads
                                if rje.exists(fasout): fx += 1
                                else: self.warnLog('Read fasta file {0} not created'.format(fasout))
                        self.printLog('\r#HICBAM','Extracted {0} read IDs from {1} BAM flank regions: {2} tmp files'.format(rje.iLen(readids),qh,rje.iStr(tx)))
                        if self.dev():
                            self.printLog('\r#HICBAM','Extracted {0} BAM flank regions and reads: {1} BAM and {2} fasta files'.format(qh,rje.iStr(bx),rje.iStr(fx)))
                        ## ~ Identify region-identical paired IDs ~ #
                        #!# NOTE: This would not work with pureflanks=F
                        tmpfile = '{0}{1}.tmp'.format(bamdir[qh],qh)
                        open(tmpfile,'w').writelines(readids)
                        uniqid = string.split(os.popen("sort %s | uniq -c | awk '$1 == 1' | awk '{print $2;}' | tee %s" % (tmpfile,uniqfile)).read())
                        self.printLog('\r#HICBAM','{0} of {1} {2} BAM flank read IDs found in one flank only.'.format(rje.iLen(uniqid),rje.iLen(readids),qh))
                    ## ~ Reduce to non-unique read IDs per region ~ ##
                    #!# This is very slow. Consider speeding up with diff and/or forking #!#
                    nx = 0; ix = 0; ex = 0.0; etot = bed.entryNum()
                    for entry in bed.entries():
                        self.progLog('\r#HICBAM','Reducing {0} BAM flank IDs to partial pairs: {1:.2f}%'.format(qh,ex/etot)); ex += 100.0
                        regfile = '%s%s.id' % (bamdir[qh],entry['name'])
                        if entry['name'] in self.dict['FlankMap'][qh]: regfile = '%s%s.id' % (bamdir[qh],self.dict['FlankMap'][qh][entry['name']])
                        if rje.exists(regfile) and not self.force():
                            regids = string.split(open(regfile,'r').read())
                        elif extract:
                            regtmp = '%s%s.tmp' % (bamdir[qh],entry['name'])
                            regids = rje.listDifference(string.split(open(regtmp,'r').read()),uniqid)
                            if flankpairs[entry['name']] != entry['name']:
                                regtmp2 = '%s%s.tmp' % (bamdir[qh],flankpairs[entry['name']])
                                regids = rje.listDifference(regids,string.split(open(regtmp2,'r').read()))
                            #regids = string.split(os.popen("diff %s %s | grep '^<' | awk '{print $2;}' | tee %s" % (regtmp,uniqfile,regfile)).read())
                            open(regfile,'w').write('\n'.join(regids))
                        else:
                            self.warnLog('No ID file for {0} and no BAM for extraction'.format(entry['name']))
                            regids = []
                        if not regids: nx += 1
                        ix += len(regids)
                        entry['score'] = min(1000,len(regids))
                        entry['pair'] = flankpairs[entry['name']]
                    self.printLog('\r#HICBAM','Reduced {0} BAM flank IDs to partial pairs -> {1} read IDs; {2} regions without reads'.format(qh,rje.iStr(ix),rje.iStr(nx)))
                    #i# Cleanup
                    if not self.dev() and not self.debugging():
                        for entry in bed.entries():
                            regtmp = '%s%s.tmp' % (bamdir[qh],entry['name'])
                            if rje.exists(regtmp): os.unlink(regtmp)

            ### ~ [3] Calculate pair overlaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                for field in ['HiCScore','BestFlank5','BestFlank3']: gdb.addField(field)
                hicscore = self.getStr('HiCScore')
                gx = 0.0; gtot = gdb.entryNum()
                for entry in gdb.entries():
                    self.progLog('\r#HICBAM','Calculating {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,gx/gtot)); gx += 100.0
                    regions = [entry['GapFlank5'],entry['GapFlank3']]
                    pentry = self.hicPair(regions,qh,pdb,entry['SynBad'])
                    entry['HiCScore'] = pentry[hicscore]
                    #?# Should this be moved to be part of the correction method?
                    if entry['HiCScore'] == 'NA':
                        entry['BestFlank5'] = '?'
                        entry['BestFlank3'] = '?'
                        continue
                    elif entry['HiCScore'] >= 0.5:
                        entry['BestFlank3'] = entry['GapFlank3']
                        entry['BestFlank5'] = entry['GapFlank5']
                        flankcheck = []
                    elif entry['SynBad'] in ['Syn', 'Aln', 'Long', 'InvFix', 'Dup', 'InvDupFix'] and self.getStrLC('HiCMode') != 'full':
                        #i# For the same-sequence gaps, just compare to others in same sequence
                        #i# Note that the reciprocal search should still identify a better match for out of place contigs
                        flankcheck = bed.indexEntries('seqname',entry['seqname'])
                    else:
                        flankcheck = bed.entries()
                    #i# 5' Join
                    bestscore = entry['HiCScore']
                    bestjoin = entry['GapFlank5']
                    for jentry in flankcheck:
                        regions = [entry['GapFlank3'],jentry['name']]
                        pentry = self.hicPair(regions,qh,pdb,'Check')
                        if pentry and pentry[hicscore] != 'NA' and pentry[hicscore] > bestscore:
                            bestscore = pentry[hicscore]
                            bestjoin = jentry['name']
                    entry['BestFlank5'] = bestjoin
                    #i# 3' Join
                    bestscore = entry['HiCScore']
                    bestjoin = entry['GapFlank3']
                    for jentry in flankcheck:
                        regions = [entry['GapFlank5'],jentry['name']]
                        pentry = self.hicPair(regions,qh,pdb,'Check')
                        if pentry and pentry[hicscore] != 'NA' and pentry[hicscore] > bestscore:
                            bestscore = pentry[hicscore]
                            bestjoin = jentry['name']
                    entry['BestFlank3'] = bestjoin
                self.printLog('\r#HICBAM','Calculating {1} {0} BAM flank ID pair overlaps complete'.format(qh,rje.iStr(pdb.entryNum())))

            ### ~ [4] Add random/full data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                #i# NOTE: Random pairs are very unlikely to share IDs. Need a good way to score match and then compare.
                #i# Need to take into account the number of IDs in each flank. Naively assume half match in each direction.
                #?# Should we have some positive control regions of the adjacent regions within a chunk
                #i# NOTE: Duplication gaps are going to cause issues due to repeated sequences?
                #!# Check the setting for HiCMode
                if self.getStrLC('HiCMode') == 'random':
                    regnames = bed.index('name').keys()
                    #!# Add some kind of safety check, e.g. if total combinations < 2x replicates -> do them all!
                    replicates = 10000
                    randx = 0
                    reglist1 = rje.randomList(regnames)
                    reglist2 = rje.randomList(regnames)
                    while randx < replicates:
                        self.progLog('\r#RANDOM','Generating random {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,100.0*randx/replicates))
                        if not reglist1:
                            reglist1 = rje.randomList(regnames)
                            reglist2 = rje.randomList(regnames)
                        regions = [reglist1.pop(0),reglist2.pop(0)]
                        self.hicPair(regions,qh,pdb,'Rand')
                        randx += 1
                    self.printLog('\r#RANDOM','Generated {1} random {0} BAM flank ID pair overlaps.'.format(qh,rje.iStr(replicates)))
                elif self.getStrLC('HiCMode') == 'full':
                    regnames = bed.index('name').keys()
                    replicates = len(regnames)
                    randx = 0
                    for region1 in regnames:
                        self.progLog('\r#RANDOM','Generating full {0} BAM flank ID pair overlaps: {1:.2f}%'.format(qh,100.0*randx/replicates))
                        for region2 in regnames:
                            regions = [region1,region2]
                            self.hicPair(regions,qh,pdb,'Rand')
                        randx += 1
                    self.printLog('\r#RANDOM','Generated {1} full BAM flank ID pair overlaps.'.format(qh,rje.iStr(replicates)))

            ### ~ [5] Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Without HiC data, this will only have partial data
            for qh in ['qry','hit']:
                for tname in ['flanks','contigs','hicpairs']:
                    table = self.dbTable(qh,tname)
                    if tname == 'hicpairs' and table and pcount[qh] == table.entryNum():
                        self.printLog('#SAVE','Skipping saving of {0} table - no new pairs'.format(tname))
                        continue
                    elif tname == 'hicpairs': #!# Add hiczero=T/F
                        self.progLog('\r#DROP','Dropping zero-overlaps from {0}...'.format(tname))
                        table.dropEntriesDirect('pairs',[0])
                    if table: table.saveToFile()
            return True
        except:
            self.errorLog('%s.synBadHiCMap error' % self.prog())
            return False
#########################################################################################################################
    def hicPair(self,regions,qh,pdb=None,ptype=''):  ### Returns the hicpairs entry for pair of regions.
        '''
        Returns the hicpairs entry for pair of regions.
        >> regions:list [flank1,flank2]
        >> qh:str qry or hit
        >> pdb:Table hicpairs table (will grab if None)
        >> ptype:str [''] = Pair type for entry
        :return: entry
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if regions[0] == regions[1]: return None
            if not pdb: pdb = self.dbTable(qh,'hicpairs')
            bamdir = {'qry':self.getStr('HiCDir1'),'hit':self.getStr('HiCDir2')}[qh]
            if bamdir.lower() in ['','none']: bamdir = '{0}.{1}flanks/'.format(pdb.baseFile(),qh)
            regions.sort()
            pentry = pdb.data((regions[0],regions[1]))
            if pentry: return pentry
            if self.getBool('FullHiC'): return None
            regfile = '%s%s.id' % (bamdir,regions[0])
            if regions[0] in self.dict['FlankMap'][qh]: regfile = '%s%s.id' % (bamdir,self.dict['FlankMap'][qh][regions[0]])
            idlist1 = string.split(open(regfile,'r').read())
            regfile = '%s%s.id' % (bamdir,regions[1])
            if regions[1] in self.dict['FlankMap'][qh]: regfile = '%s%s.id' % (bamdir,self.dict['FlankMap'][qh][regions[1]])
            idlist2 = string.split(open(regfile,'r').read())
            pairs = rje.listIntersect(idlist1,idlist2)
            pentry = {'region1':regions[0],'region2':regions[1],'id1':len(idlist1),'id2':len(idlist2),'pairs':len(pairs),'type':ptype,'wtscore':0.0,'score':0.0}
            if pentry['pairs']:
                pentry['score'] = 2.0 * pentry['pairs'] / (pentry['id1']+pentry['id2'])
                pentry['wtscore'] = (0.5 * pentry['pairs'] / pentry['id1']) + (0.5 * pentry['pairs'] / pentry['id2'])
            elif min(pentry['id1'],pentry['id2']) < 1:
                pentry['score'] = 'NA'
                pentry['wtscore'] = 'NA'
            elif self.getBool('FullHiC'): return None
            pdb.addEntry(pentry)
            return pentry
        except:
            self.errorLog('%s.hicPair error' % self.prog())
            return None
#########################################################################################################################
    def mapFlanksToGenome(self,qh):   ### Runs GABLAM of flanks against assembly to map old flanks onto new assembly.
        '''
        Runs GABLAM of flanks against assembly to map old flanks onto new assembly.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            gabbase = '{0}.{1}.flanks'.format(basefile,qh)
            wanted = ['.local.tdt','.map.tdt']
            maptable = '{0}.flanks.map'.format(qh)
            oldflanks = {'qry':'MapFlanks1','hit':'MapFlanks2'}[qh]
            if not self.getStrLC(oldflanks):
                self.printLog('#FLANKS','No {0}=FILE flank mapping provided'.format(oldflanks.lower()))
                return None
            elif not rje.exists(self.getStr(oldflanks)):
                raise IOError('{0}={1} not found!'.format(oldflanks.lower(), self.getStr(oldflanks)))

            ### ~ [1] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('{0} GABLAM mapping'.format(qh),line='-')
            if not self.force() and rje.checkForFiles(wanted,basename=gabbase,log=self.log,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Flanks local hits mapping found: skipping GABLAM (force=F)')
                mdb = db.addTable(name=maptable,mainkeys=['GapFlank'],ignore=[])
                return mdb
            ## ~ [1a] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            loctable = '{0}.local.tdt'.format(gabbase)
            flankfas = '{0}.{1}.flanks.fasta'.format(basefile,qh)
            #!# Add more intelligent minloclen filter based on contig size?
            if self.force() or not rje.exists(loctable):
                gabcmd = ['seqin={0}'.format(flankfas),'searchdb={0}'.format(self.getStr(oldflanks)),'mapper=minimap','minlocid=99','minloclen=0','basefile={0}'.format(gabbase),'uniqueout=F']
                gabobj = gablam.GABLAM(self.log,self.cmd_list+gabcmd+['debug=F'])
                gabobj.gablam()
            ## ~ [2b] Load tables and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mdb = db.addEmptyTable(maptable,['GapFlank','OldFlank'],['GapFlank'])
            ldb = db.addTable(loctable,mainkeys=['Qry','Hit','AlnNum'],name=loctable,ignore=[],expect=True)
            ldb.dataFormat({'AlnNum':'int','Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            #i# Qry is the new flank. Hit is the old flank. Should not need to consider the Strand.
            for entry in ldb.sortedEntries('Length',reverse=True):
                if entry['Qry'] in mdb.data(): continue
                mdb.addEntry({'GapFlank':entry['Qry'],'OldFlank':entry['Hit']})
            mdb.saveToFile()
            return mdb
        except: self.errorLog('%s.mapFlanksToGenome error' % self.prog()); return None
#########################################################################################################################
    ### <7> ### Assembly Map Methods                                                                     #
#########################################################################################################################
    #!# Will want to move these functions out to rje_assemblies.py to make available for other tools.
    #i# Assembly maps are stored in self.list['qry'] and self.list['hit']
    #i# Form: ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
    #i# Saved to a *.txt file with one line per sequence and a "seqname = " prefix
    #i# Also save the contig sequences themselves to a matching *.contigs.fasta file.
    #i# Have methods for:
    # - (1) Generating the *.txt file from the map, renaming where required.
    # - (2) Taking a *.txt file and *.contigs.fasta file and making a new fasta file.
    #i# Edit this assembly map with the corrections rather than the actual assembly -> then make new fasta
    #i# Can re-map flanks and thus define any new assembly in this context too.
#########################################################################################################################
    def makeAssemplyMap(self,qh):   ### Populates self.list[qh] with the assembly map.
        '''
        Populates self.list[qh] with the assembly map. Assembly map will be in the form:
        ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            qhbase = '{0}.{1}'.format(basefile,qh)
            self.list[qh] = []  # Assembly map
            self.headLog('{0} Assembly Map generation'.format(qh),line='-')
            ### ~ [1] Generate Map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.dbTable(qh,'contigs')
            gdb = self.dbTable(qh,'gap')
            gdb.index('GapFlank5')
            cx = 0.0; ctot = cdb.entryNum()
            seqname = None
            seqmap = []
            for entry in cdb.entries(sorted=True):
                self.progLog('\r#MAP','Generating assembly map: {0:.1f}%'.format(cx/ctot)); cx += 100.0
                if seqmap and seqname != entry['seqname']:
                    #self.bugPrint('\n{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                    self.list[qh] += seqmap
                    seqmap = []
                # Generate seqmap
                seqname = entry['seqname']
                if not seqmap: seqmap = ['|']
                seqmap += [entry['flank5'],'>',entry['name'],'>',entry['flank3']]
                gtype = entry['synbad'].split('-')[1]
                if gtype == 'End':
                    seqmap += ['|']
                else:
                    seqmap += [':{0}:{1}:'.format(gtype,gdb.indexEntries('GapFlank5',entry['flank3'])[0]['gaplen'])]
            if seqmap:
                #self.bugPrint('\n{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                self.list[qh] += seqmap
            self.printLog('\r#MAP','Generated {0} assembly map from {1} contigs.'.format(qh,rje.iStr(ctot)))
            return True
        except:
            self.errorLog('%s.makeAssemplyMap error' % self.prog())
            return False
#########################################################################################################################
    def OLDmakeAssemplyMap(self,qh):   ### Populates self.list[qh] with the assembly map and saves contigs fasta.
        '''
        Populates self.list[qh] with the assembly map and saves contigs fasta. Assembly map will be in the form:
        ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand. This should make it
        easy to identify the right things for edits, but also make some nice assembly map strings if needed.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            qhbase = '{0}.{1}.'.format(basefile,qh)
            self.list[qh] = []  # Assembly map
            #i# 'contigs.fasta' = Fasta file of the contigs, to be used for correction
            #   'contigs.txt' = SeqName = Map for each scaffold
            #   'contigs.tdt' = Table of contigs = Made during HiC mapping and flank extraction
            #       -> seqname,start,end,ctglen,name,flank5,flank3,synbad
            wanted = ['contigs.fasta','contigs.txt','contigs.tdt']
            self.headLog('{0} Assembly Map generation'.format(qh),line='-')
            maptxt = '{0}contigs.txt'.format(qhbase)
            ### ~ [1] Check/Load Map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.force() and rje.checkForFiles(wanted,basename=qhbase,log=self.log,cutshort=False,ioerror=False,missingtext='Not found.'):
                self.printLog('#MAP','{0} assembly map found: loading (force=F)'.format(qh))
                try:
                    for line in open(maptxt,'r').readlines():
                        #self.bugPrint(line)
                        data = string.split(line)
                        if not data: continue
                        if data[1] == '=': self.list[qh] += data[2:]
                    return True
                except:
                    self.errorLog('Problem processing {0}contigs.txt: will regenerate'.format(qhbase))
                    self.list[qh] = []
            ### ~ [2] Generate Map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdb = self.dbTable(qh,'contigs')
            gdb = self.dbTable(qh,'gap')
            gdb.index('GapFlank5')
            cx = 0.0; ctot = cdb.entryNum()
            seqname = None
            seqmap = []
            rje.backup(self,maptxt)
            OUT = open(maptxt,'w')
            for entry in cdb.entries(sorted=True):
                self.progLog('\r#MAP','Generating assembly map: {0:.1f}%'.format(cx/ctot))
                if seqmap and seqname != entry['seqname']:
                    #self.bugPrint('{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                    OUT.write('{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                    self.list[qh] += seqmap
                    seqmap = []
                # Generate seqmap
                seqname = entry['seqname']
                if not seqmap: seqmap = ['|']
                seqmap += [entry['flank5'],'>',entry['name'],'>',entry['flank3']]
                gtype = entry['synbad'].split('-')[1]
                if gtype == 'End':
                    seqmap += ['|']
                else:
                    seqmap += [':{0}:{1}:'.format(gtype,gdb.indexEntries('GapFlank5',entry['flank3'])[0]['gaplen'])]
            if seqmap:
                #self.bugPrint('{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                OUT.write('{0} = {1}\n'.format(seqname, ' '.join(seqmap)))
                self.list[qh] += seqmap
            OUT.close()
            self.printLog('\r#MAP','Generated {0} assembly map -> {1}'.format(qh,maptxt))

            return True
        except:
            self.errorLog('%s.corrections error' % self.prog())
            return False
#########################################################################################################################
    def getContigSeqObj(self,qh,reload=False):   ### Returns contigs.fasta seqlist object
        '''
        Returns contigs.fasta seqlist object.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqobj = '{0}.contigs'.format(qh)
            if seqobj in self.obj and self.obj[seqobj] and self.obj[seqobj].seqNum() and not reload: return self.obj[seqobj]
            ### ~ [1] Load contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fasin = '{0}.{1}.fasta'.format(self.baseFile(),seqobj)
            seqcmd = self.cmd_list + ['seqin={0}'.format(fasin),'summarise=T']
            self.obj[seqobj] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
            if not self.obj[seqobj].seqNum(): raise IOError('Failed to load sequences from {0}'.format(fasin))
            return self.obj[seqobj]
        except:
            self.errorLog('%s.getContigSeqObj error' % self.prog()); raise
#########################################################################################################################
    def saveAssemblyMaps(self,qh,mapname='map'):    ### Saves map text and fasta files for qry or hit data
        '''
        Saves map text and fasta files for qry or hit data.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newacc = self.getStrLC({'qry':'NewAcc1','hit':'NewAcc2'}[qh])
            if newacc:
                newacc = self.getStr({'qry':'NewAcc1','hit':'NewAcc2'}[qh])
            mapout = '{0}.{1}.{2}.txt'.format(self.baseFile(),qh,mapname)
            if not self.outputAssemblyMap(self.list[qh],mapout,newacc): return False
            fasout = rje.baseFile(mapout) + '.fasta'
            contigs = self.getContigSeqObj(qh)
            if not self.fastaFromAssemblyMap(contigs,mapout,fasout): return False
            return True
        except:
            self.errorLog('%s.saveAssemblyMaps error' % self.prog()); raise
#########################################################################################################################
    def outputAssemblyMap(self,maplist,mapout,newacc=None):   ### Generate assembly text file using assembly map.
        '''
        Generate assembly text file using assembly map. Will regenerate names and change gap lengths if required.
        This file will then be used with a contigs fasta file to generate an updated assembly fasta output.
        >> maplist:list = ['|',Flank1,'>',CtgName,'>',Flank2,':SynBad:GapLen:',Flank1 ... ,'|'] with '<' for -ve Strand.
        >> mapout:str = output file name for assembly
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Assembly Map output',line='-')
            #i# Sort scaffolds by length.
            #i# If not newacc, will use the first seqname from that assembly map sequence.
            #i# Where there is a clash, will add a .X counter.
            # gaplen=INT : Set new standardised gap length for SynBad assembly output [500]
            newgaplen = max(0,self.getInt('GapLen'))

            ### ~ [1] Generate scaffolds from map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Setup scaffolds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scaffolds = []  # List of (scafflen, tempname, mapstr) for each scaffold.
            maptxt = ' '.join(maplist)
            #self.bugPrint(maptxt)
            scaffmaps = maptxt.split('|')
            namecount = {}
            sx = 0.0; stot = len(scaffmaps)
            #self.bugPrint('Assembly map split into {0} scaffold chunks'.format(stot))
            for scaff in scaffmaps:
                smap = scaff.split()
                if not scaff or not smap:
                    stot -= 1
                    continue
                self.progLog('\r#SCAFF','Parsing scaffolds from assembly map: {0:.1f}%'.format(sx/stot)); sx += 100.0
                #self.bugPrint(scaff)
                #self.bugPrint(smap[0])
                scaffname = '.'.join(smap[0].split('.')[:-1])
                #self.bugPrint(scaffname)
                scafflen = 0
                #i# Calculate scaffold length, fix gap lengths and check formatting
                #i# Each cycle should be Flank1,'>',CtgName,'>',Flank2 then either ':SynBad:GapLen:' or end
                i = 0
                while i < len(smap):
                    mapel = ' '.join(smap[i:i+6])
                    try:
                        #!# Move check for coherent map elements into a function -> check edits too
                        #i# Check Dirn
                        if smap[i+1] != smap[i+3] or smap[i+1] not in '<>': raise ValueError('Assembly map strand formatting error!')
                        contig = smap[i+2]
                        cspan = map(int,contig.split('.')[-1].split('-'))
                        scafflen += (cspan[1] - cspan[0] + 1)
                        #i# Gap?
                        i += 5
                        if i >= len(smap): break
                        part = smap[i]
                        if part[:1] != ':': raise ValueError('Assembly map gap formatting error!')
                        gap = part.split(':')
                        if len(gap) != 4: raise ValueError('Assembly map gap split formatting error!')
                        gaplen = int(gap[2])
                        if newgaplen and newgaplen != gaplen:
                            gap[2] = str(gaplen)
                            smap[i] = ':'.join(gap)
                            gaplen = newgaplen
                        scafflen += gaplen
                        i += 1
                    except:
                        self.errorLog('Assembly map generation error')
                        raise ValueError('Problem with map element: {0}'.format(mapel))
                scaffstr = ' '.join(smap)
                scaffolds.append((scafflen,scaffname,scaffstr))
                if scaffname not in namecount: namecount[scaffname] = 0
                namecount[scaffname] += 1
            self.printLog('\r#SCAFF','Parsing {0} scaffolds from assembly map complete'.format(rje.iStr(stot)))
            if len(scaffolds) != stot: raise ValueError('Scaffold count mismatch!')

            ## ~ [1b] Sort scaffolds by length and rename if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maxx = max(namecount.values())
            for name in namecount.keys():
                if namecount[name] == 1: namecount.pop(name)
            scaffolds.sort()
            sorted = scaffolds
            scaffolds = []
            while sorted:
                if newacc: newname = '{0}{1}'.format(newacc,rje.preZero(len(sorted),len(scaffolds)))
                else:
                    newname = sorted[0][1]
                    if newname in namecount:
                        accx = rje.preZero(namecount[newname],maxx)
                        namecount[newname] -= 1
                        newname = '{0}{1}'.format(sorted[0][1],accx)
                scaff = sorted.pop(0)
                scaffolds.append((newname,scaff[0],scaff[2]))
            self.printLog('#SCAFF','{0} scaffolds sorted and renamed.'.format(rje.iLen(scaffolds)))

            ## ~ [1c] Sort scaffolds by name and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            scaffolds.sort()
            sx = 0.0
            rje.backup(self,mapout)
            OUT = open(mapout,'w')
            for scaff in scaffolds:
                self.progLog('\r#MAP','Outputting assembly map: {0:.1f}%'.format(sx/stot)); sx += 100.0
                #self.bugPrint('{0} ({1}) = {2}\n'.format(scaff[0],scaff[1],scaff[2]))
                OUT.write('{0} ({1}) = {2}\n'.format(scaff[0],scaff[1],scaff[2]))
            OUT.close()
            self.printLog('\r#MAP','Outputted assembly map -> {0}'.format(mapout))

            return True
        except:
            self.errorLog('%s.outputAssemblyMap error' % self.prog())
            return False
#########################################################################################################################
    def fastaFromAssemblyMap(self,contigs,maptxt,fasout):   ### Generate assembly fasta from assembly map and contig sequences.
        '''
        Generate assembly fasta from assembly map and contig sequences.
        >> contigs:SeqList = contig sequences.
        >> maptxt:str = Assembly map text file.
        >> fasout:str = output file name.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Assembly Map fasta generation',line='-')
            ctgdict = {}    # {name:sequence}
            for seq in contigs.seqs():
                ctgdict[contigs.shortName(seq)] = contigs.getSeq(seq)[1]
            self.printLog('#CTG','Built sequence dictionary for {0} contigs.'.format(rje.iStr(contigs.seqNum())))
            rje.backup(self,fasout)

            ### ~ [1] Load Map and output fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            FASOUT = open(fasout,'w'); sx = 0
            for mapline in open(maptxt,'r').readlines():
                self.progLog('\r#MAPFAS','Converting assembly map to fasta: {0} scaffolds'.format(rje.iStr(sx)))
                mapdata = mapline.split(' = ')
                if len(mapdata) < 2: continue
                elif len(mapdata) > 2: raise ValueError('Input mapping format error!')
                scaffname = mapdata[0].split()[0]
                scafflen = rje.matchExp('\((\d+)\)',mapdata[0].split()[-1])
                if scafflen: scafflen = int(scafflen[0])
                scaffseq = ''
                smap = mapdata[1].split()
                smap = string.join(smap)
                smap = string.split(smap)
                #self.deBug('\n"{0}" "{1}" "{2}" ... "{3}" "{4}"'.format(smap[0],smap[1],smap[2],smap[-2],smap[-1]))
                i = 0; cx = 0
                while i < len(smap):
                    if not smap[i] or not rje.matchExp('(\S)',smap[i]): i += 1; continue
                    try:
                        #i# Check Dirn
                        if smap[i+1] != smap[i+3] or smap[i+1] not in '<>': raise ValueError('Assembly map formatting error!')
                        ctgseq = ctgdict[smap[i+2]]
                        if smap[i+1] == '<': ctgseq = rje_sequence.reverseComplement(ctgseq)
                        scaffseq += ctgseq; cx += 1
                        #i# Gap?
                        i += 5
                        if i >= len(smap): break
                        if not smap[i] or not rje.matchExp('(\S)',smap[i]): i += 1; continue
                        if not smap[i].startswith(':'):
                            self.warnLog('Map assembly element should be gap but is not: "{0}"'.format(smap[i]))
                            i += 1
                            continue
                        gap = smap[i].split(':')
                        if len(gap) != 4: raise ValueError('Assembly map formatting error!')
                        gaplen = int(gap[2])
                        scaffseq += 'N' * gaplen
                        i += 1
                        #self.bugPrint('{0}: "{1}"'.format(i,' '.join(smap[i-6:i])))
                    except:
                        #self.bugPrint('\n{0}: "{1}"'.format(i,' '.join(smap[i-6:i])))
                        self.errorLog('Assembly map processing error')
                        raise ValueError('Problem with map element: "{0}"'.format(' '.join(smap[i:i+6])))
                #i# Check scaffold length
                if scafflen and scafflen != len(scaffseq):
                    self.warnLog('Scaffold {0} length mismatch: map says {1} bp but sequence parsed is {2} bp'.format(scaffname,rje.iStr(scafflen),rje.iLen(scaffseq)))
                FASOUT.write('>{0} {1} bp ({2} contigs)\n{3}\n'.format(scaffname,len(scaffseq),cx,scaffseq)); sx += 1
            self.printLog('\r#MAPFAS','Converted assembly map to fasta: {0} scaffolds -> {1}'.format(rje.iStr(sx),fasout))

            return True
        except:
            self.errorLog('%s.fastaFromAssemblyMap error' % self.prog())
            return False
#########################################################################################################################
    def couldConflict(self,qh,entry): # Whether a correction could conflict with another
        '''
        Whether a correction could conflict with another.
        >> qh:str = qry/hit
        >> correction: dict = Corrections table entry #['seqname','start','end','flank1','flank2','edit']
        :return: True/False whether this edit could cause conflicts
        '''
        #i# Return False if edit is a single contig
        if entry['flank1'] == entry['flank2']: return False
        flank1 = entry['flank1']
        flank2 = entry['flank2']
        amap = self.list[qh]
        mapi = [amap.index(flank1), amap.index(flank2)]
        mapi.sort()
        editchunk = amap[mapi[0]:mapi[1]+1]
        if len(editchunk) == 5: return False
        #i# ELse, return True
        return True
#########################################################################################################################
    def checkSubMap(self,submap,gap5='',gap3=''):   ### Checks whether a submap list is a legitimate chunk
        '''
        Checks whether a submap list is a legitimate chunk.
        >> submap:list = list of assembly map elements
        >> gap5:str = optional upstream gap position to check
        >> gap3:str = optional downstream gap position to check
        '''
        if gap5 and gap5[:1] not in ['|',':']: raise ValueError('Upstream element not an end or gap')
        if gap3 and gap3[:1] not in ['|',':']: raise ValueError('Downstream element not an end or gap')
        if len(submap) % 6 != 5:
            self.debug(submap.join())
            raise ValueError('Length of assembly submap inconsistent with whole contigs')
        return True
#########################################################################################################################
    def mapCorrections(self):  ### Loads corrections table, edits sequences and outputs corrected assembly
        '''
        Loads corrections table, edits sequences and outputs corrected assembly.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seqObjSetup()
            seqobj = {}
            seqobj['qry'] = self.obj['SeqList1']
            seqobj['hit'] = self.obj['SeqList2']
            base = {'qry':rje.baseFile(self.getStr('Genome1'),strip_path=True),'hit':rje.baseFile(self.getStr('Genome2'),strip_path=True)}
            basefile = self.baseFile(strip_path=True)
            db = self.db()
            goodedits = ['invert']
            for qh in ['qry','hit']:
                cdb = self.dbTable(qh,'corrections')
                if cdb: db.deleteTable(cdb)
                cdb = db.addTable(mainkeys=['seqname','start','end','edit'],name='{0}.corrections'.format(qh),replace=True,ignore=[],expect=True)

            ### ~ [1] Check edits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Check for conflicting edits, i.e. those sharing a flank? Might want to have an explicit ordering added.
            #!# Add start and end positions to table for correct sorting. Use this for overlaps? (See older code.)
            #!# Make sure it is reported if edits are skipped. (Filter edit types first and add commandline option.)
            for qh in ['qry','hit']:
                cdb = self.dbTable(qh,'corrections')
                prev = None
                for entry in cdb.entries(sorted=True):
                    if entry['edit'] not in goodedits:
                        entry['synbad'] = 'Not implemented'
                        continue
                    if not self.couldConflict(qh,entry): continue
                    if prev and prev['end'] >= entry['start'] and prev['seqname'] == entry['seqname']:
                        self.warnLog('Dropped correction {0} {1}-{2} due to overlapping edit region.'.format(entry['seqname'],entry['start'],entry['end']))
                        entry['synbad'] = 'Blocked'
                        continue
                    prev = entry
                cdb.indexReport('seqname')
                cdb.indexReport('edit')
                for etype in cdb.index('edit'):
                    if etype not in goodedits:
                        self.warnLog('Edit type "{0}" not implemented.')

            ### ~ [2] Update assembly map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
                cdb = self.dbTable(qh,'corrections')
                ## ~ [2a] Inversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#INV','Processing {0} inversions'.format(qh))
                editx = 0
                for entry in cdb.indexEntries('edit','invert'):
                    if entry['synbad'] == 'Blocked': continue
                    if self.mapInversion(qh,entry['flank1'],entry['flank2']): entry['synbad'] = 'Fixed'
                    else: entry['synbad'] = 'Failed'
                    editx += 1
                self.printLog('\r#INV','Processed {0} inversions: {1} edits'.format(qh,editx))

            ### ~ [3] Output updated map and fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                self.saveAssemblyMaps(qh,mapname='synbad')
                cdb.saveToFile(backup=False)

            return True
        except:
            self.errorLog('%s.corrections error' % self.prog())
            return False
#########################################################################################################################
    def mapInversion(self,qh,flank1,flank2):    ### Inverts assembly map between flank1 and flank2
        '''
        Inverts assembly map between flank1 and flank2.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh]
            checklen = len(amap)
            ### ~ [1] Invert ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapi = [amap.index(flank1), amap.index(flank2)]
            mapi.sort()
            #!# Add check for correct orientation prior to inversion #!#
            #!# Add check for coherent map element
            invchunk = amap[mapi[0]:mapi[1]+1]
            try: self.checkSubMap(invchunk,amap[mapi[0]-1],amap[mapi[1]+1])
            except:
                self.errorLog('Inversion failed',quitchoice=False)
                return False
            invchunk.reverse()
            for i in range(len(invchunk)):
                if invchunk[i] == '>': invchunk[i] = '<'
                elif invchunk[i] == '<': invchunk[i] = '>'
            self.list[qh] = amap[:mapi[0]] + invchunk + amap[mapi[1]+1:]
            if len(self.list[qh]) != checklen: raise ValueError('Length of assembly map has changed during inversion')
            return True
        except:
            self.errorLog('%s.mapInversion error' % self.prog())
            return False
#########################################################################################################################
    def mapDeletion(self,qh,flank1,flank2):    ### Deletes assembly map between flank1 and flank2
        '''
        Deletes assembly map between flank1 and flank2.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            amap = self.list[qh]
            ### ~ [1] Delete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add check for correct orientation prior to inversion #!#
            #!# Add check for coherent map element
            mapi = [amap.index(flank1), amap.index(flank2)]
            mapi.sort()
            try: self.checkSubMap(amap[mapi[0]:mapi[1]+1],amap[mapi[0]-1],amap[mapi[1]+1])
            except:
                self.errorLog('Deletion failed',quitchoice=False)
                return False
            self.list[qh] = amap[:mapi[0]] + amap[mapi[1]+1:]
            return True
        except:
            self.errorLog('%s.mapDeletion error' % self.prog())
            return False
#########################################################################################################################
### End of SECTION II: SynBad Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: SynBad(mainlog,['basefile=synbad']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
