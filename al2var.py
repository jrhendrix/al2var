'''
FILE:	al2var.py
AUTHOR:	J.R. Hendrix
URL: 	http://stronglab.org
DESC:	This script aligns paired-end reads to a referece assembly.
		Reports alignment rate
		Reports variant (estimated error) rate
RECQ:	bcftools, bowtie2, samtools

ACKNOWLEDGMENTS:
		Sean Beagle - Provided basic alignment workflow
'''

# IMPORT FROM PYTHON STANDARD LIBRARY
import argparse
import logging
import operator
import os
import subprocess
import sys

from subprocess import Popen, PIPE

# IMPORT FROM PROJECT LIBRARY
#from handler import Dir, File, Fasta, Fastq


class Dir:
	""" Base class for system directories """

	def __init__(self, path):
		self._path = None
		self.path = path

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isdir(value):
			self._path = value
		else:
			raise NotADirectoryError(value)

	@property
	def dirname(self):
		return self.path.strip("/").split("/")[-1]

	@property
	def children(self):
		children = [Dir(os.path.join(self.path, subdir)) 
			for subdir in os.listdir(self.path) 
			if os.path.isdir(os.path.join(self.path, subdir))]
		if len(children) > 0:
			return children
		else:
			return None

	@property
	def files(self):
		files = [File(os.path.join(self.path, file))
			for file in os.listdir(self.path)
			if os.path.isfile(os.path.join(self.path, file))]
		if len(files) > 0:
			return files
		else:
			return None

	def join(self, *args):
		return os.path.join(self.path, *args)

	def make_subdir(self, *args):
		""" Makes recursive subdirectories from 'os.path.join' like arguments """
		subdir = self.join(*args)
		return self.make(subdir)

	@classmethod
	def make(cls, path):
		try:
			os.makedirs(path)
			return cls(path)
		except FileExistsError:
			return cls(path)

	def __repr__(self):
		return self.path
	
	

class File:
	""" Base class for all file-types """

	def __init__(self, path, file_type=None):
		self._path = None
		self.path = path
		self.file_type = file_type

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isfile(value):
			self._path = value
		else:
			raise FileNotFoundError(value)

	@property
	def dir(self):
		return Dir(os.path.dirname(self.path))

	@property
	def filename(self):
		return os.path.basename(self.path)

	@property
	def file_prefix(self):
		return self.filename.split(".")[0]

	@property
	def extension(self):
		return self.filename.split(".")[-1]
	

class Fasta(File):
	"""File Type"""
	extensions = ('fasta', 'fa', 'fna', 'faa')

	def __init__(self, file):
		super().__init__(file, file_type="fasta")

		self.stats = {
			'contigs': None, 
			'bp': None, 
			'avg_contig_len': None, 
			'std_contig_len': None, 
			'nt_freq': None,
			'aa_freq': None}

	@property
	def bp(self, stat='bp'):
		if self.stats[stat] is None:
			self.get_stats()
		return self.stats[stat]

	@property
	def nt_freq(self, stat='nt_freq'):
		if self.stats[stat] is None:
			self.get_stats_nt()
		return self.stats[stat]

	@property
	def aa_freq(self, stat='aa_freq'):
		if self.stats[stat] is None:
			self.get_stats_aa()
		return self.stats[stat]
	
	
	def get_stats_nt(self):
			n = []
			nt_freq = {'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N': 0, '-': 0}
			for record in SeqIO.parse(self.path, self.file_type):
				n.append(len(record))
				for nt in record.seq:
					nt_freq[nt.upper()] += 1

			self.stats['bp'] = sum(n)
			self.stats['contigs'] = len(n)
			self.stats['nt_freq'] = nt_freq

	def get_stats_aa(self):
			n = []
			aa_freq = {
			'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 
			'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0, 
			'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0,
			'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0,
			'-': 0}
			for record in SeqIO.parse(self.path, self.file_type):
				n.append(len(record))
				for aa in record.seq:
					aa_freq[aa.upper()] += 1

			self.stats['bp'] = sum(n)
			self.stats['contigs'] = len(n)
			self.stats['aa_freq'] = aa_freq


class Fastq(File):
	""" Create FastQ file objects """
	extensions = ('fastq', 'fq', 'fastq.gz', 'fq.gz')

	def __init__(self, file):
		super().__init__(file, file_type="fastq")
		self._stats = {'set': False, 
		'read_count': None, 
		'read_lengths': None,
		'read_seqs': None,
		'read_quals': None, 
		'avg_read_len': None, 
		'std_read_len': None}

	@property
	def reads(self, stat='reads'):
		self.get_stats()
		return self._stats[stat]

	@property
	def read_lengths(self, stat='read_lengths'):
		self.get_stats()
		return self._stats[stat]

	@property
	def avg_read_len(self, stat='avg_read_len'):
		self.get_stats()
		return self._stats[stat]

	@property
	def std_read_len(self, stat='std_read_len'):
		self.get_stats()
		return self._stats[stat]

	def get_stats(self):
		if not self._stats['set']:
			with gzip.open(self.path, "rt") as handle:
				n = [len(record) for record in SeqIO.parse(handle, self.file_type)]
				self._stats['reads'] = len(n)
				self._stats['avg_read_len'] = mean(n)
				self._stats['std_read_len'] = stdev(n, self._stats['avg_read_len'])
				self._stats['read_lengths'] = n
				self._stats['set'] = True
	

# INITIATE LOG
LOG = logging.getLogger('log_file')


def configure(args, ref_id, samp_id):
	''' Configure used directories and log file '''
	global OUTDIR, LOG, REFDIR

	# CREATE OUTPUT DIRECTORIES
	basedir = Dir(args.output_path)
	OUTDIR = basedir.make_subdir(args.output_directory)
	basename = '_X_'.join((ref_id, samp_id))
	REFDIR = OUTDIR.make_subdir("indexes", "indexes")

	# INITIATE LOG FILE
	LOG.setLevel(logging.DEBUG)
	#formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
	formatter = logging.Formatter('%(message)s')
	LOG_File = OUTDIR.join("al2var.log")
	file_handler = logging.FileHandler(LOG_File)
	file_handler.setFormatter(formatter)

	LOG.addHandler(file_handler)

	intro = 'Running Al2Var script\n'
	LOG.info(intro)

	LOG.info('CONFIGURATION...')
	LOG.info(f'Aligning {samp_id} to {ref_id}.\n')
	LOG.info(f'OUTDIR {OUTDIR.path}')

	return basename, LOG_File


def make_bowtie_db(args, reference, ref_id):
	''' Set up reference directory '''

	LOG.info("BUILDING BOWTIE DATABASE...")
	#ref_path = '/'.join((REFDIR.path, REFDIR.dirname))
	#scratch = REFDIR.make_subdir("bowtie2", REFDIR.dirname)
	#ref_path = ref.path


	# MAKE BOWTIE DB
	try:
		#outpath = '/'.join((REFDIR.path, ref_id))
		command = ["bowtie2-build", reference.path, REFDIR.path]
		process = subprocess.run(command, capture_output=True)
		if process.returncode == 0:
			LOG.info(f"... DONE")
		if process.returncode == 1:
			LOG.info(f"... FAILED : {process.stderr}")
			exit(process.stderr)

	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)



def write_to_log(args, std):
	''' Interate over std input and write to log '''

	for line in std:
		l = ''.join(('\t', str(line).split("'")[1][:-2]))
		LOG.info(l)

def run_subprocess(args, command):
	''' Runs a command as subprocess '''

	process = subprocess.Popen(command, stdout=PIPE, stderr=PIPE)
	stdout = process.stdout
	stderr = process.stderr

	write_to_log(args, stdout)
	write_to_log(args, stderr)

	process.wait()


def run_pipeline(args, program, basename, reference, read_arr):
	
	# DIRECTORIES
	bam_dir = OUTDIR.make_subdir("bam")
	fq_dir = OUTDIR.make_subdir("unconc")
	index_dir = OUTDIR.make_subdir("indexes")
	mpileup_dir = OUTDIR.make_subdir("mpileup")
	sam_dir = OUTDIR.make_subdir("sam")
	vcf_dir = OUTDIR.make_subdir("vcf")

	# FILES
	sam_file = sam_dir.join(f"{basename}.sam")
	fq_file = fq_dir.join(f"{basename}.fq")
	bam_file = bam_dir.join(f"{basename}_000.bam")
	sorted_bam_file = bam_file.replace("_000.bam", "_SORTED.bam")
	mpileup_file = mpileup_dir.join(f"{basename}.mpileup")
	vcf_file = vcf_dir.join(f"{basename}_000.vcf")
	#filtered_vcf_file = vcf_file.replace("_000.vcf", "_FILTERED.vcf")
	var_vcf_file = vcf_file.replace("_000.vcf", "_var.vcf")

	
	# MAP ASSEMBLY TO REFERENCE (CREATE SAM FILE)
	LOG.info("MAPPING TO REFERENCE GENOME...")
	try:
		num_mm = str(args.num_mismatch)
		len_seed = str(args.length_seed)

		LOG.info(f'\tNumber of mismatches allowed:\t{num_mm}')
		LOG.info(f'\tLength of seed sequence:     \t{len_seed}')

		if program == 'bowtie2' and len(read_arr) == 2:
			pair1 = read_arr[0]
			pair2 = read_arr[1]

			command = [
			"bowtie2", "-x", REFDIR.path, 
			"-1", pair1.path, "-2", pair2.path,
			"--un-conc", fq_file, 
			"-S", sam_file,
			"-N", num_mm,
			"-L", len_seed]
		'''
		elif program == 'bowtie2' and len(read_arr) == 1:
			reads = read_arr[0]
			command = [
			"bowtie2", "-x", REFDIR.path, "-U", 
			reads.path, 
			"-S", sam_file,
			"-N", num_mm,
			"-L", len_seed]	
		'''

		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)
	


	
	# CONVERT SAM --> BAM
	try:
		LOG.info("CONVERTING SAM FILE TO BAM FILE...")
		command = ['samtools', 'view', sam_file, '-u', '-o', bam_file]
		#print(command)
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)
	
	# SORT BAM FILE
	try:
		LOG.info("SORTING BAM FILE...")
		command = ['samtools', 'sort', bam_file, '-o', sorted_bam_file]
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)


	# INDEX SORTED BAM FILE
	try:
		LOG.info("INDEXING SORTED BAM FILE...")
		command = ['samtools', 'index', sorted_bam_file]
		run_subprocess(args, command)
		LOG.info('... DONE')	
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)	
	
	
	# CONVERT SORTED BAM --> MPILEUP
	try:
		LOG.info("CONVERTING SORTED BAM FILE TO MPILEUP FILE...")
		command = [
			'samtools', 'mpileup', '-u', '-t', 'DP', '-t', 'AD', '-t', 'ADF', '-t', 'ADR',
			'-f', reference.path, '-o', mpileup_file, sorted_bam_file]
		#print(command)
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)

	
	# CONVERT MPILEUP --> VCF
	try:
		LOG.info("CONVERTING MPILEUP FILE TO VCF FILE...")
		command = ['bcftools', 'call', '-c', mpileup_file, '-o', vcf_file]
		#print(command)
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)
	

	# EXTRACT VARS
	try:
		LOG.info("EXTRACTING VARS FROM VCF FILE...")

		command1 = ['awk', '$5 != "."', vcf_file]
		command2 = ['awk', '!/^#/']

		# Run LHS of pipe
		process1 = subprocess.Popen(command1, stdout=subprocess.PIPE, shell=False)
		stdout = process1.stdout
		stderr = process1.stderr
		#process1.wait()
		LOG.info(f'\tFinished process 1')

		# Run RHS of pipe
		process2 = subprocess.Popen(command2, stdin=stdout, stdout=subprocess.PIPE, shell=False)
		process1.stdout.close()
		LOG.info(f'\tRan process 2')
		#process1.wait()
		stdout = process2.stdout
		stderr = process2.stderr

		# Write VAR data
		f = open(var_vcf_file, 'w')
		f2 = open(vcf_file, 'r')
		for line in f2:
			if line.startswith('#'):
				f.write(line)
			else:
				break
		f2.close()

		numVar = 0
		if stdout != None:
			for line in stdout:
				l = ''.join((str(line).split("'")[1][:-2], '\n'))
				l = l.replace('\\t', '\t')
				numVar = numVar + 1
				f.write(l)
		f.close()

		#process2.stdout.close()
		#process2.wait()
		LOG.info(f"\tNUMBER OF VARIANTS: {numVar}")
		LOG.info('... DONE')
		

	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)


	# CALCULATE VARIANT RATE
	try:
		LOG.info("CALCULATING VARIANT RATE...")
		# Count genome length
		f = open(args.reference, 'r')
		bases = 0
		for line in f:
			if line.startswith('>'):
				continue
			bases = bases + len(line)
		f.close()

		# Calculate rate: vars/genome = x/100,000
		v = round((numVar/bases)*100000, 3)
		LOG.info(f"\tGENOME LENGTH: {bases}")
		LOG.info(f"\tVARIANT RATE: {v}/100,000 bases")
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)
	

	if args.cleanup:
		dirList = [fq_dir, mpileup_dir, sam_dir, index_dir]
		if not args.keep_bam:
			dirList.append(bam_dir)

		for d in dirList:
			cmd = ' '.join(('rm -r', d.path))
			os.system(cmd)

	return(v)


def return_align_rate(args, LOG_File):
	''' Extracts and returns the Bowtie reported Overall Alignment Rate '''

	command = ['grep', 'alignment', LOG_File]

	(stdout, stderr) = subprocess.Popen(command, stdout=PIPE, stderr=PIPE).communicate()
	alignRate = str(stdout).replace("b'", "")
	alignRate = alignRate.replace('\\t', '\t')
	alignRate = alignRate.replace("\\n'", '\n')

	return alignRate


def report_stats(args, aR, vR):
	fname = '.'.join((args.savename, 'txt'))
	outfile = '/'.join((OUTDIR.path, fname))
	f1 = open(outfile, 'w')

	if aR is not None:
		val = str(aR).strip().split(' ')[0]
	else:
		val = 'Could not be calculated'
	entry = ' '.join(('Alignment rate:', val))+'\n'
	f1.write(entry)

	if vR is not None:
		val = ''.join((str(vR), '/100,000 bases'))
	else:
		val = 'Could not be calculated'
	entry = ' '.join(('Variant rate:', val))
	f1.write(entry)

	f1.close()
	

def bowtie2(args, command):
	# SET UP REFERENCE
	reference = Fasta(args.reference)
	ref_id = reference.filename.split('.')[0]

	# HANDLE PAIRED READS
	pair1 = Fastq(args.pair1)
	pair2 = Fastq(args.pair2)
	samp_id = pair1.filename.split('-pair1')[0]
	sample_arr = [pair1, pair2]

	# CONFIGURE
	basename, LOG_File = configure(args, ref_id, samp_id)
	LOG.info(f'REFERENCE : {reference.filename}')
	LOG.info(f'SAMPLE FWD: {pair1.filename}')
	LOG.info(f'SAMPLE RVS: {pair2.filename}\n')

	# SET UP REFERENCE
	make_bowtie_db(args, reference, ref_id)

	# RUN ALIGNMENT PIPELINE DATABASE
	var_rate = run_pipeline(args, 'bowtie2', basename, reference, sample_arr)

	# GET ALIGNMENT RATE
	align_rate = return_align_rate(args, LOG_File)
	
	# REPORT STATS TO REPORT FILE
	report_stats(args, align_rate, var_rate)

def main(program):
	command = 'Command: %s' % ' '.join(sys.argv)
	print(f'Running: ', command)
	#configure(args, command)
		
	args.func(args, command)


if __name__== "__main__":

	cwd = os.getcwd()

	parser = argparse.ArgumentParser(description='program description')
	subparsers = parser.add_subparsers(dest="cmd", help='available actions')
	#subparsers = parser.add_subparsers(title="build", dest="Build a de bruijn graph", help='available actions')
	#subparsers = parser.add_subparsers(title="stat", dest="Get basic statistics", help='available actions')
	subparsers.required = True

	# PARSER : ROOT
	__version__ = "0.0.0"
	parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))



	# PARSER : ROOT
	'''
	parent_parser = argparse.ArgumentParser(prog='al2var', add_help=False)
	parent_parser.add_argument('-b', '--keep_bam', help='Do no remove bam files during cleanup step. Only usable with -c flag', default=False, action='store_true')
	parent_parser.add_argument('-c', '--cleanup', help='Enable automatic cleanup of intermediary files', default=False, action='store_true')
	parent_parser.add_argument('-l', '--length_seed', default=22, help='Sets length of seed. Smaller values increase sensitivity and runtime', type=int)
	#parent_parser.add_argument('-mk_off', '--make_db_off', default=False, action='store_true', help='Disable make Bowtie database function.')
	#parent_parser.add_argument('-mod', '--module_load', default=False, action='store_true', help='Enable module load function.')
	parent_parser.add_argument('-n', '--num_mismatch', default=0, help='Number of mismatches allowed in the seed sequence. Setting to 1 increases sensitivity and runtime', choices=['0','1'])
	parent_parser.add_argument('-o', '--output_directory', default='out_al2var', help='Prefix of output directory', type=str)
	parent_parser.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	#parent_parser.add_argument('-ref_dir', '--reference_directory', default=None, help='Reference Directory', type=str)
	parent_parser.add_argument('-r', '--reference', default=None, help='Reference sequence', type=str)
	#parent_parser.add_argument('-rate', '--report_rate', default=False, action='store_true', help='Return the alignment rate')
	parent_parser.add_argument('-s', '--savename', default='report', help='Savename for report file.')	
	subparsers = parent_parser.add_subparsers(help='sub-command help')
	'''
	
	# DEFINE SUBPARSERS
	parser_bowtie2 = subparsers.add_parser('bowtie2')
	parser_bowtie2.set_defaults(func=bowtie2)

	# PARSER : BOWTIE2
	#bowtie = subparsers.add_parser('bowtie2', help='Align paired-end Illumina reads to reference', parents=[parent_parser])
	parser_bowtie2.add_argument('-1', '--pair1', help='fwd trimmed paired-end reads to assemble', required=True)
	parser_bowtie2.add_argument('-2', '--pair2', help='rev trimmed paired-end reads to assemble', required=True)
	parser_bowtie2.add_argument('-b', '--keep_bam', help='Do no remove bam files during cleanup step. Only usable with -c flag', default=False, action='store_true')
	parser_bowtie2.add_argument('-c', '--cleanup', help='Enable automatic cleanup of intermediary files', default=False, action='store_true')
	parser_bowtie2.add_argument('-l', '--length_seed', default=22, help='Sets length of seed. Smaller values increase sensitivity and runtime', type=int)
	parser_bowtie2.add_argument('-n', '--num_mismatch', default=0, help='Number of mismatches allowed in the seed sequence. Setting to 1 increases sensitivity and runtime', choices=['0','1'])
	parser_bowtie2.add_argument('-o', '--output_directory', default='out_al2var', help='Prefix of output directory', type=str)
	parser_bowtie2.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parser_bowtie2.add_argument('-r', '--reference', default=None, help='Reference sequence', type=str)
	parser_bowtie2.add_argument('-s', '--savename', default='report', help='Savename for report file.')	



	args = parser.parse_args()

	main(args)










