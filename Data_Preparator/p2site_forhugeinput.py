#!/usr/bin/python

# import global modules
import os
import sys
import argparse
import logging
import numpy
import re
import pandas
import multiprocessing

# Module metadata variables
__author__ = "Jose Rodriguez"
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "1.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"

# Global var
NUM_CPUs = multiprocessing.cpu_count()
# NUM_CPUs = 20

class corrector:
    '''
    Make the calculations in the workflow
    '''
    # The list of variant positions for PTMs
    PTM_POS = (-1, 0, 1)

    def __init__(self, infile, indb, search, decoy=None):
        # create an index from the fasta sequence, if apply
        self.indb = self._read_fasta(indb, decoy)
        # create the report with the relationship table
        # self.data = numpy.genfromtxt(infile, delimiter="\t", dtype="string", comments=None, skip_header=1)
        self.data = pandas.read_csv(infile, delimiter="\t", na_values=['NA'], header=None)
        # type of search
        self.search = search
        # header for the output
        self.rst_header = ['[FASTAProteinDescription]', '[Sequence]', '[Tags]']

    def _read_fas(self, fp, dec):
        indb = {}
        discard = False
        for line in fp:
            line = line.rstrip()
            if line:
                if dec and line.startswith(">"+dec):
                    discard = True
                elif line.startswith(">"):
                    discard = False
                    fid = line
                    indb[fid] = ''
                elif not discard:
                    indb[fid] += line
        return indb

    def _read_fasta(self, infile, decoy):
        '''
        Read a fasta file creating a hash variable with the identifiers
        Only with UniProt.
        Discard the DECOY sequences
        '''
        indb = {}
        with open(infile) as fp:
            indb = self._read_fas(fp, decoy)
        return indb

    def _delete_delta_mass(self, seq):
        '''
        Delete the delta mass from the sequence
        '''
        s = None
        if self.search == "open":
            s = re.sub('\[[^\]]*\]', '', seq)
        else:
            s = re.sub('[^a-zA-Z]', '', seq)
        return s

    def _find_ptm_pos(self, seq):
        '''
        Find the list of ptm positions (-1,0,1) from the peptide:
            ptm in the first aa => we only return two postions
            ptm in the last aa => we only return two postions
            ptm in the middle => we only return three postions
        '''
        pos = []
        if self.search == "open":
            for m in re.finditer('\[', seq):
                pos.append(m.start())
        else:
            i = 0 # fix the increased shift
            for m in re.finditer('[^a-zA-Z]', seq):
                pos.append(m.start()-i)
                i += 1
        return pos

    def _find_peptide(self, seq, pep):
        '''
        Find the peptide in the protein sequence.
        get the list of (start,end) peptides jumping by two
        '''
        # for each index of peptide in the sequence
        pos = []
        for m in re.finditer(pep, seq):
            pos.append( m.start() )
        return pos


    def _get_pep_site(self, pdesc, sequence, peptide, pos):
        '''
        Extract the multiple sites for each PTM
        '''
        # first delete [delta_mass]
        pep_seq = self._delete_delta_mass(peptide)
        # get the length of peptide (without delta_mass)
        pep_len = len(pep_seq)
        # find the list of ptm positions (-1,0,1) from the peptide
        #   ptm in the first aa => we only return two postions
        #   ptm in the last aa => we only return two postions
        #   ptm in the middle => we only return three postions
        pep_pos = self._find_ptm_pos(peptide)
        # find the peptide in the protein sequence
        # get the list of (start,end) peptides jumping by two
        seq_pos = self._find_peptide(sequence, pep_seq)
        # for each list of peptide matches in the sequence
        # get the range of positions and residues
        results = []
        for i in self.PTM_POS:
            for seq_p in seq_pos:
                for pep_p in pep_pos:
                    if ( pep_p+i >= 1 and pep_p+i <= pep_len ):
                        ptm_p = seq_p + pep_p + i
                        ptm_r = sequence[ptm_p-1]
                        results.append( (str(ptm_p)+"_"+ptm_r+"_"+pdesc, peptide, str(i)) )
        return results

    def _get_peptide_sites(self, row):
        '''
        Extract the multiple sites for each PTM
        '''
        # extract columns
        pid = row[0]
        ptm = row[1]
        # delete marks
        pid = re.sub('^\"',"", pid)
        pid = re.sub('\"$',"", pid)
        # get the protein sequence
        if pid in self.indb:
            sequence = self.indb[pid]
            rst = self._get_pep_site(pid, sequence, ptm, self.PTM_POS)
            # print( rst )
            # print( pandas.DataFrame(rst) )
            # print( "-------" )
        else:
            rst = [("XXX_NaN_"+pid, ptm, "NaN")]
            # rst = pandas.DataFrame( [("XXX_NaN_"+pid, ptm, "NaN")] )

        return rst
    
    def get_peptide_sites(self):
        '''
        Extract the multiple sites for each PTM
        '''
        # get the peptide sites for every row
        results = apply_by_multiprocessing(self.data, self._get_peptide_sites, axis=1, workers=NUM_CPUs)

        return results

    def write_to_file(self, results, outfile):
        '''
        Print to CSV
        '''
        # function that update the file with row data
        def _print_rst(r, f):
            rst = "\n".join("{}\t{}\t{}".format(*el) for el in r) + "\n"
            f.write(rst)

        # print header
        rst = "\t".join(self.rst_header) + "\n"
        f = open(outfile, "w")
        f.write(rst)
        f.close()
        # update the file with the row data using pandas apply
        f = open(outfile, "a")
        results.apply(_print_rst, args=(f,) )
        f.close()

def _apply_df(args):
    df, func, kwargs = args
    return df.apply(func, **kwargs)

def apply_by_multiprocessing(df, func, **kwargs):
    workers = kwargs.pop('workers')
    pool = multiprocessing.Pool(processes=workers)
    result = pool.map(_apply_df, [(d, func, kwargs)
            for d in numpy.array_split(df, workers)])
    pool.close()
    return pandas.concat(result)

def main(args):
    '''
    Main function
    '''
    # extract temporal working directory...
    if args.tmpdir:
        tmpdir = args.tmpdir
    # otherwisae, get directory from input files
    else:
        tmpdir = os.path.dirname(os.path.realpath(args.relfile))+"/tmp"
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    logging.info('create corrector object')
    c = corrector(args.relfile, args.indb, args.search, args.decoy)

    logging.info('calculate the unique protein')
    p = c.get_peptide_sites()
    
    logging.info('print output file')
    c.write_to_file(p, tmpdir+"/p2q_rpos_v2.tsv")


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Create the relationship table for the positions of PTM in the protein sequence',
        epilog='''Examples:
        python  src/SanXoT/p2site.py 
          -f test/p2site/Human_jul14.curated.fasta
          -r test/p2site/p2q_rels_open.tsv
          -t test/p2site/open
          -s open
          -d DECOY

        python  src/SanXoT/p2site.py
          -f test/p2site/Mouse_April27_2016_contam.fasta
          -r test/p2site/p2q_rels_close.tsv
          -t test/p2site/close
          -d DECOY
        ''',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r',  '--relfile',  required=True, help='Peptide2Protein relationship table. Separeted in tabular')
    parser.add_argument('-f',  '--indb',     required=True, help='Database file with the protein sequences (FASTA format). At the moment, only applied for UniProt files')
    parser.add_argument('-d',  '--decoy',    help='Tag that identify the DECOY sequences. If not giving, we use all the sequences')
    parser.add_argument('-s',  '--search',   choices=['open', 'close'], default="close", help='Type of search')
    parser.add_argument('-t',  '--tmpdir',   help='Output directory')
    parser.add_argument('-l',  '--logfile',  help='Output file with the log tracks')
    parser.add_argument('-vv', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # set-up logging
    scriptname = os.path.splitext( os.path.basename(__file__) )[0]

    # add logging handler, formatting, etc.
    # the logging is info level by default
    # add filehandler
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    if args.logfile:
        # logfile = os.path.dirname(os.path.realpath(args.relfile)) + "/"+ scriptname +".log"
        ch = logging.FileHandler(args.logfile)
    else:
        ch = logging.StreamHandler()
    ch.setFormatter( logging.Formatter('%(asctime)s - %(levelname)s - '+scriptname+' - %(message)s', '%m/%d/%Y %I:%M:%S %p') )
    logger.addHandler(ch)

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')