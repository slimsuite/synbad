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
Version:      0.7.0
Last Edit:    15/04/21
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
    * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different contigs, 1+ of which is not a chromosome scaffold.
    * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
    * `Inv` = `Inversion` = Flanking hits are on alternative strands.
    * `InvBrk` = `Inversion Breakpoint` = Reclassified Inversion that appears to be out of place.
    * `InvDupFix` = `Fixed Inversion Duplication` = Inversion that becomes a Duplication when inverted. (See `*.corrections.tdt`)
    * `InvFix` = `Fixed Inversion` = Inversion that becomes collinear when inverted. (See `*.corrections.tdt`)
    * `Long` = `Long Syntenic` = Gap flanked by collinear Qry/Hit pairs but distance is greater than `maxsynspan=INT`  (default 25 kb).
    * `Null` = No mapping between genomes for that gap.
    * `Span` = `Spanned` = Any gaps with Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
    * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 25 kb).
    * `Tran` = `Translocation` = `SynSpan` indicates matches are on different contigs.
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
    gablam=X        : Optional prefix for GABLAM search [defaults to basefile=X]
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
    ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
import rje, rje_obj, rje_db, rje_rmd, rje_seqlist, rje_sequence
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
    # [ ] : Add hicbam file loading and flank/join assessements.
    # [ ] : Add final output where Qry and Hit are replaced for Genome1/2 output.
    # [ ] : Add modified fasta output following edits.
    # [ ] : Build LocusFixer chassis and use for re-assembling duplication gap regions.
    # [ ] : Execute inversions and output updated fasta file.
    # [ ] : Empirically establish the best values for maxoverlap, maxsynskip and maxsynspan.
    # [ ] : Option to add DepthCharge gaps where there might be misassemblies.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SynBad', '0.7.0', 'April 2021', '2020')
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
#i# For now, InvFix are not in goodgaps to avoid nested inversions etc. Recommend re-running on fixed fasta.
goodgaps = ['Syn', 'Aln', 'Span', 'Div', 'Ins', 'Long']
skipgaps = ['Syn', 'Aln', 'Span', 'InvFix', 'InvDupFix', 'Long', 'Div', 'Dup', 'Ins']
fraggaps = ['Brk', 'Inv', 'InvBrk', 'Frag', 'Tran']
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
    - GABLAM=X        : Optional prefix for GABLAM search [defaults to basefile=X]
    - GapMode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    - Genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
    - Genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
    - PAF1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    - PAF2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]

    Bool:boolean
    - BestPair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - Fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
    - Update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]

    Int:integer
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

    Dict:dictionary    

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
        self.strlist = ['BAM1','BAM2','Chr1','Chr2','GABLAM','GapMode','Genome1','Genome2','PAF1','PAF2']
        self.boollist = ['BestPair','DocHTML','Fragment','Update']
        self.intlist = ['MaxOverlap','MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan','SpannedFlank','SynReadSpan']
        self.numlist = ['MinLocID']
        self.filelist = []
        self.listlist = ['FragTypes','Reads1','Reads2','ReadType1','ReadType2']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'GapMode':'gapspan'})
        self.setBool({'BestPair':False,'Fragment':False,'Update':False})
        self.setInt({'MaxOverlap':500,'MaxSynSkip':4,'MaxSynSpan':25000,'MinLocLen':1000,'MinReadSpan':1,'SpannedFlank':0,'SynReadSpan':5})
        self.setNum({'MinLocID':50.0})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.list['FragTypes'] = fraggaps
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
                self._cmdReadList(cmd,'str',['Chr1','Chr2','GABLAM','GapMode'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['BAM1','BAM2','Genome1','Genome2','PAF1','PAF2'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BestPair','DocHTML','Fragment','Update'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxOverlap','MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan','SpannedFlank','SynReadSpan'])   # Integers
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
        * `Tran` = `Translocation` = `SynSpan` indicates matches are on different contigs.
        * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different contigs, 1+ of which is not a chromosome scaffold.
        * `Term` = `Terminal` = Gap is between a local alignment and the end of the query sequence.
        * `Span` = `Spanned` = Any gaps with Aligned or Syntenic rating that are spanned by at least `synreadspan=INT` reads.
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
        gablam=X        : Optional prefix for GABLAM search [defaults to basefile=X]
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
        ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        bestpair=T/F    : Whether to restrict the paired output to the top scaffold pairs [False]
        update=T/F      : Whether to reload compressed qry and hit tables but re-run additional compression [False]
        force=T/F       : Whether to force regeneration of SynBad results tables [False]
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

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
    def docHTML(self):  ### Generate the Diploidocus Rmd and HTML documents.                                        # v0.1.0
        '''Generate the Diploidocus Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            info = self.log.obj['Info']
            prog = '%s V%s' % (info.program,info.version)
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
            #!# Replace this with documentation text?
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
    ### <3> ### Main SynBad Methods                                                                                     #
#########################################################################################################################
    def diploidocusGapRun(self):    ### Runs Diploidocus and returns the gap tables
        '''
        Runs Diploidocus and returns the gap tables
        :return: (cdb1,cdb2)
        '''
        try:  ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            dcmd = ['bam=', 'paf=', 'reads=', 'readtype=ont']
            # !# Add checking of seqnames read in for gaps and local hits and warn if none match #!#

            ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DIPLOIDOCUS GAP ANALYSIS', line='=')
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
                    cdb1.setStr({'Name': 'qrygap'})
                    self.db().list['Tables'].append(cdb1)
                else:
                    cdb1 = db.addTable('%s.checkpos.tdt' % base1, mainkeys=['seqname', 'start', 'end'], name='qrygap', ignore=[], expect=True)
            if not cdb1:
                if rundip:
                    self.errorLog('Diploidocus checkpos run appears to have failed', printerror=False, quitchoice=self.i() >= 0)
                else:
                    self.printLog('#READS', 'No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base1))
                if not rje.checkForFiles(filelist=['.gaps.tdt'], basename=base1, log=self.log):
                    seqcmd = self.cmd_list + cmd1 + ['summarise=T', 'gapstats=T', 'raw=F', 'dna=T']
                    self.obj['SeqList1'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
                cdb1 = db.addTable('%s.gaps.tdt' % base1, mainkeys=['seqname', 'start', 'end'], name='qrygap', ignore=[], expect=True)
                cdb1.addField('Span0', evalue=0)
            # checkpos1 = '{0}.checkpos.tdt'.format(base1)
            # cdb1 = db.addTable(checkpos1,mainkeys=['seqname','start','end'],name='qrygap',ignore=[],expect=True)
            cdb1.dataFormat({'seqlen': 'int', 'start': 'int', 'end': 'int', 'Span0': 'int'})
            cdb1.addField('GapSpan', evalue='.')
            cdb1.addField('SynSpan', evalue='.')
            cdb1.addField('SynBad', evalue='Null')
            cdb1.index('seqname')

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
                    cdb2.setStr({'Name': 'hitgap'})
                    self.db().list['Tables'].append(cdb2)
                else:
                    cdb2 = db.addTable('%s.checkpos.tdt' % base2, mainkeys=['seqname', 'start', 'end'], name='hitgap', ignore=[], expect=True)
            if not cdb2:
                if rundip:
                    self.errorLog('Diploidocus checkpos run appears to have failed', printerror=False, quitchoice=self.i() >= 0)
                else:
                    self.printLog('#READS', 'No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base2))
                if not rje.checkForFiles(filelist=['.gaps.tdt'], basename=base2, log=self.log):
                    seqcmd = self.cmd_list + cmd2 + ['summarise=T', 'gapstats=T', 'raw=F', 'dna=T']
                    self.obj['SeqList2'] = rje_seqlist.SeqList(self.log, seqcmd + ['autoload=T', 'seqmode=file', 'autofilter=F'])
                    # dip2.cmd_list.append('gapstats')    # Try to run automatically if possible
                    # seqin = dip2.seqinObj()
                    # if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=base2,log=self.log):
                    #     seqin.setBool({'Raw':False,'GapStats':True,'DNA':True})
                    #     seqin.str['SeqType'] = 'dna'
                    #     seqin.summarise()
                cdb2 = db.addTable('%s.gaps.tdt' % base2, mainkeys=['seqname', 'start', 'end'], name='hitgap', ignore=[], expect=True)
                cdb2.addField('Span0', evalue=0)
            # checkpos2 = '{0}.checkpos.tdt'.format(base2)
            # cdb2 = db.addTable(checkpos2,mainkeys=['seqname','start','end'],name='hitgap',ignore=[],expect=True)
            cdb2.dataFormat({'seqlen': 'int', 'start': 'int', 'end': 'int', 'Span0': 'int'})
            cdb2.addField('GapSpan', evalue='.')
            cdb2.addField('SynSpan', evalue='.')
            cdb2.addField('SynBad', evalue='Null')
            cdb2.index('seqname')

            return(cdb1,cdb2)
        except: self.errorLog('%s.diploidocusGapRun error' % self.prog()); raise
#########################################################################################################################
    def runGABLAM(self):   ### Runs GABLAM of two assemblies
        '''
        Runs GABLAM of two assemblies and loads results into database tables.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            wanted = ['qry.tdt','hit.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=None,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad qry.tdt and hit.tdt files found: skipping GABLAM (force=F)')
                return True
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#

            ### ~ [1] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GABLAM SEARCH',line='=')
            gabbase = basefile
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
    def synBad(self):   ### Main SynBad run method
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#

            ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.headLog('DIPLOIDOCUS GAP ANALYSIS',line='=')
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
            (cdb1,cdb2) = self.diploidocusGapRun()

            ### ~ [2] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.headLog('GABLAM SEARCH',line='=')
            if not self.runGABLAM(): return False

            ### ~ [3] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.synBadMapping(cdb1,cdb2): return False
            if not self.synBadCompress(): return False
            if not self.tidyTables(): return False
            if not self.saveTables(): return False

            ## ~ [3a] Summary reports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
            base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
            self.headLog(base1,line='-')
            self.db('qrygap').indexReport('SynBad')
            self.headLog(base2,line='-')
            self.db('hitgap').indexReport('SynBad')

            #!# Move corrections and fragmentation outside this method and make self-sufficient (loading tables if missing)

            ### ~ [4] Update assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.corrections(): return False

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
    def synBadMapping(self,cdb1,cdb2):   ### Main SynBad gap mapping and compression
        '''
        Main SynBad gap mapping and compression.
        >> cdb1: Table of genome 1 gaps
        >> cdb2: Table of genome 2 gaps
        << True/False
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            base1 = rje.baseFile(self.getStr('Genome1'), strip_path=True)
            base2 = rje.baseFile(self.getStr('Genome2'), strip_path=True)
            wanted = ['qry.tdt','hit.tdt']
            if not self.force() and rje.checkForFiles(wanted,basename=db.baseFile()+'.',log=None,cutshort=True,ioerror=False,missingtext='Not found.'):
                self.printLog('#SKIP','Main SynBad qry.tdt and hit.tdt files found: skipping SynBad mapping (force=F)')
                return True

            ### ~ [1] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD GAP MAPPING',line='=')
            self.printLog('#NOTE','Genome1 ({0}) is always Qry and Genome2 ({1}) is always Hit'.format(base1,base2))
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
                gapdb = self.db('{0}gap'.format(qh.lower()))
                spanfield = 'Span{0:d}'.format(self.getInt('SpannedFlank'))
                spancheck = spanfield in gapdb.fields()
                if self.getInt('SynReadSpan') > 0 and not spancheck:
                    self.warnLog('synreadspan={0} but "{0}" not found in fields. Check spannedflank=INT.'.format(self.getInt('SynReadSpan'),spanfield))
                spancheck = spancheck and self.getInt('SynReadSpan') > 0
                if spancheck: gapdb.dataFormat({spanfield:'int'})
                #self.debug(gapdb)
                altdb = {cdb1:cdb2,cdb2:cdb1}[gapdb]
                locdb = self.db(qh.lower())
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
                        # * `Translocation` = `SynSpan` indicates matches are on different contigs.
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
            qrygapdb = self.db('qrygap')
            hitgapdb = self.db('hitgap')
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
            wanted = ['qryfull.tdt','qry.tdt','qrypairs.tdt','qrygap.tdt','hitfull.tdt','hit.tdt','hitpairs.tdt','hitgap.tdt']
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
                locdb = self.db(qh.lower())
                locdb.info['Name'] = '{0}full'.format(qh.lower())
                locdb.saveToFile()
                locdb.info['Name'] = qh.lower()
            ## ~ [0c] Add gaps to tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                alt = {'Qry':'Hit','Hit':'Qry'}[qh]
                gapdb = self.db('{0}gap'.format(qh.lower()))
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
                locdb = self.db(qh.lower())
                skipx = 0
                self.headLog('SynBad {0} quick compression: merge adjacent local hits'.format(qh),line='~')
                self.synBadCompressTable(qh,locdb,skipx,quick=True)

            ### ~ [2] Remove internal odd entries and merge flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('MaxSynSkip') > 0:
                for qh in ('Qry','Hit'):
                    alt = {'Qry':'Hit','Hit':'Qry'}[qh]
                    locdb = self.db(qh.lower())
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
    def saveTables(self):   ### SynBad tidying of gaps with clear fixes
        '''
        Saves the output tables, once all processing is complete.
        '''
        try:### ~ [1] Save gap and local tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SAVE SYNBAD TABLES',line='=')
            for qh in ('Qry','Hit'):
                self.db('{0}gap'.format(qh.lower())).saveToFile()
                locdb = self.db(qh.lower())
                locdb.dropField('AlnNum')
                outfields = locdb.fields()[0:]
                try:
                    outfields.remove('{0}5'.format(qh))
                    outfields.remove('{0}3'.format(qh))
                except:
                    pass    # Will need to update pairs tables!
                locdb.saveToFile(savefields=outfields)

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
                topdb = db.copyTable(self.db(qh.lower()),'{0}pairs'.format(qh.lower()),replace=True,add=True)
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
                topdb.saveToFile()

            return True
        except: self.errorLog('%s.saveTables error' % self.prog()); return False
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
            tdb = self.db().addEmptyTable('corrections',['seqname','start','end','edit','details'],['seqname','start','end','edit'])

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
                locdb = self.db(qh.lower())
                gapdb = self.db('{0}gap'.format(qh.lower()))
                skipx = self.getInt('MaxSynSkip') + 2
                self.headLog('SynBad {0} translocation/breakpoint assessment'.format(qh),line='~')
                self.synBadDivergence(qh,locdb,gapdb,skipx,quick=True)

            ### ~ [2] Inversions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Inversions consist of a region flanked by two Inv QryGaps. In the simplest case, there are no intervening
            #i# gaps, but the following gaps can also be skipped: Syn, Span, Dup. (Aln should not appear!)
            for qh in ('Qry','Hit'):
                locdb = self.db(qh.lower())
                gapdb = self.db('{0}gap'.format(qh.lower()))
                #skipx = self.getInt('MaxSynSkip') + 2
                self.headLog('SynBad {0} inversion assessment'.format(qh),line='~')
                self.synBadInversions(qh,locdb,gapdb,0,quick=True)
                self.synBadInvertedDuplications(qh,locdb,gapdb,0,quick=True)

            ### ~ [3] Update hitgaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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

            #i# ~ [!] Breakpoints ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# When a Brk gap is encountered, look for something that would fit in the gap. If something is found,
            #i# and has Tran/Brk flanks, replace them with 'Mis' ("Misassembly").


            ### ~ [!] Translocations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Find regions that are flanked by two Trans gaps with the same Hit. Then find another Trans gap that is
            #i# adjacent and close enough to the flanked region to move.


            #i# ~ [!] Dodgy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Identify dodgy regions that have non-syntenic gaps and duplications etc. without a clear fix.


            #!# Sort this out so it is done properly!
            tdb.saveToFile(append=not self.force())
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
            gapdb = self.db('{0}gap'.format(alt.lower()))
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
        #># Generate a new fragments table that skips through Syn, Aln, Span, InvFix, InvDupFix, Long and Div gaps. Sort by size?
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
            gapdb = self.db('{0}gap'.format(qh.lower()))
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
            syntable = '{0}blocks'.format(qh.lower())
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
            tdb = self.db('corrections')
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
                self.deBug('\n{0} -| {1} -> {2} |- {3} ?'.format(flankstart,invstart,invend,flankend))
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
                replace = toskip[1:-1]
                invstart = replace[0][qstart]
                invend = replace[-1][qend]
                for lentry in replace:
                    lentry['Strand'] = {'+':'-','-':'+','.':'.'}[lentry['Strand']]
                    tmplen = lentry[qend] - lentry[qstart]
                    lentry[qstart] = invstart + (invend - lentry[qend])
                    lentry[qend] = lentry[qstart] + tmplen
                replace.reverse()
                self.bugPrint(len(entries))
                entries = entries[:i1+2] + replace + entries[i2-1:]
                self.debug(len(entries))
                #i# Make sure that this inversion is documentation -> will need to check/execute later (or reverse)
                self.printLog('#FIX','Invert {0}:{1}-{2}'.format(merge1[qh],invstart,invend))
                tdb.addEntry({'seqname':merge1[qh],'start':invstart,'end':invend,'edit':'invert','details':'Inversion in-place.'})
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
            tdb = self.db('corrections')
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
                self.deBug('\n{0} -| {1} -> {2} |- {3} ?'.format(flankstart,invstart,invend,flankend))

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
                    replace = toskip[1:-1]
                    invstart = replace[0][qstart]
                    invend = replace[-1][qend]
                    for lentry in replace:
                        lentry['Strand'] = {'+':'-','-':'+','.':'.'}[lentry['Strand']]
                        tmplen = lentry[qend] - lentry[qstart]
                        lentry[qstart] = invstart + (invend - lentry[qend])
                        lentry[qend] = lentry[qstart] + tmplen
                    replace.reverse()
                    self.bugPrint(len(entries))
                    entries = entries[:i1+2] + replace + entries[i2-1:]
                    self.debug(len(entries))
                    #i# Make sure that this inversion is documentation -> will need to check/execute later (or reverse)
                    self.printLog('#FIX','Invert {0}:{1}-{2}'.format(merge1[qh],invstart,invend))
                    tdb.addEntry({'seqname':merge1[qh],'start':invstart,'end':invend,'edit':'invert','details':'Inversion in-place.'})
            locdb.remakeKeys()
            self.printLog('\r#INV','Assessed complex %s table inversions: %s gaps reclassified' % (tname,rje.iStr(dx)))
            if dx: self.warnLog('NOTE: Fixed inversions will not carry over into %s table. (Fix assemblies and re-run SynBad.)' % alt)

            return dx
        except: self.errorLog('%s.synBadInvertedDuplications error' % self.prog()); raise
#########################################################################################################################
    ### <5> ### SynBad Output Methods                                                                                   #
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
            db = self.db()
            cdb = db.addTable(mainkeys=['seqname','start','end','edit'],name='corrections',replace=True,ignore=[],expect=True)
            cdb.dataFormat({'start':'int','end':'int'})
            goodedits = ['invert']

            ### ~ [1] Check edits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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

            #!# This workflow will need editing when going beyond inversions: need to be careful of clashes.
            ### ~ [2] Perform edits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for qh in ['qry','hit']:
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
            qrygap = self.db('qrygap')
            if not qrygap:
                qrygap = db.addTable('%s.gaps.tdt' % base1,mainkeys=['seqname','start','end'],name='qrygap',ignore=[],expect=True)
            qrygap.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            qryfrag = db.copyTable(qrygap,'qryfrag',replace=True,add=True)
            hitgap = self.db('hitgap')
            if not hitgap:
                hitgap = db.addTable('%s.gaps.tdt' % base2,mainkeys=['seqname','start','end'],name='hitgap',ignore=[],expect=True)
            hitgap.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            hitfrag = db.copyTable(hitgap,'hitfrag',replace=True,add=True)
            ### ~ [1] Filter gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # * `Aln` = `Aligned` = Gap is found in the middle of a local alignment to the Hit
            # * `Syn` = `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
            # * `Ins` = `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
            # * `Brk` = `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
            # * `Dup` = `Duplication` = Overlapping flanking hits on the same strand.
            # * `Inv` = `Inversion` = Flanking hits are on alternative strands.
            # * `Tran` = `Translocation` = `SynSpan` indicates matches are on different contigs.
            # * `Frag` = `Fragmentation` = `SynSpan` indicates matches are on different contigs, 1+ of which is not a chromosome scaffold.
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
            for frag in ('qryfrag','hitfrag'):
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
