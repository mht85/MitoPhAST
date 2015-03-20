##########################
#dr.mark.schultz@gmail.com
#github schultzm
#date 18/02/15
#requires python and biopython to run, can be installed by typing on the command line:
#'sudo pip install biopython'
##########################

#import the modules
from Bio import Entrez
import argparse

#set up the arguments
parser = argparse.ArgumentParser(description = "Will pull genbank files specified by accession numbers and put them into your \'pwd\' (present working directory).  To get help, type: \'python pull_genomes_as_genbank.py -h\'.  Example usage, type: \'python pull_genomes_as_genbank.py -a AB071201 GQ223286 GQ292768 KC107808 AB859775 -e user@xa.ya.za -s .gb\'")
parser.add_argument('-a', '--accession_ids', nargs ='+', help = 'Accession numbers.  If more than one accession, separate list members by whitespace', required=True)
parser.add_argument('-e', '--user_email', help = 'User email address', required = False)
#the default file extension is set in the set_suffix() function [end of script], not in the next line in argparse
parser.add_argument('-s', '--file_suffix', help = 'File extension, defaults to \'.gb\'.  Enter alternative extensions as required (e.g., \'.gbk\' or \'.gff\').', required = False)
args = parser.parse_args()

Entrez.email = args.user_email

# accession id works, returns genbank format, looks in the 'nucleotide' database:
def get_genbank(extension):
	for i in args.accession_ids:
		try:
			handle=Entrez.efetch(db='nucleotide', id=i, rettype='gbwithparts', retmode="text")
# 			if handle != None:
			with open(i+'.'+extension, 'w') as output_file:
				print 'downloading '+i+' to '+i+'.'+extension
				output_file.write(handle.read())
		except:
			print 'Accession number '+i+' not found'
			continue

		handle.close()

#sets up the file extension (suffix)
def set_suffix(suffix):
	if suffix == None:
		extension = 'gb'
		get_genbank(extension)
	else:
		extension = suffix.replace('.', '')
		get_genbank(extension)

#executes the set_suffix() function, which in turn executes the get_genbank() function
set_suffix(args.file_suffix)

length = len(args.accession_ids)
print '\nDone.   Check \'pwd\' for', length, 'new output files.\n'
