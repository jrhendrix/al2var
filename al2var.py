'''
FILE:	bar_script.py
AUTHOR:	J.R. Hendrix
URL: 	http://stronglab.org
DESC:	This script aligns reads (paired-end short, single-end short, and log)
		to a referece assembly.

ACKNOWLEDGMENTS:
		Sean Beagle - Provided basic alignment workflow
		This pipeline has been passed through the department - many people to acknowledge
		Unclear who the original author was
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
	LOG_File = OUTDIR.join("wrapper.log")
	file_handler = logging.FileHandler(LOG_File)
	file_handler.setFormatter(formatter)

	LOG.addHandler(file_handler)

	intro = 'Running Bowtie2 Alignment Runner script\n'
	LOG.info(intro)

	LOG.info('CONFIGURATION...')
	LOG.info(f'Aligning {samp_id} to {ref_id}.\n')
	LOG.info(f'OUTDIR {OUTDIR.path}')

	return basename


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



def run_pipeline(args, basename, reference, read_arr):
	
	# DIRECTORIES
	bam_dir = OUTDIR.make_subdir("bam")
	fq_dir = OUTDIR.make_subdir("unconc")
	#index_dir = OUTDIR.make_subdir("indexes")
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
	filtered_vcf_file = vcf_file.replace("_000.vcf", "_FILTERED.vcf")
	var_vcf_file = vcf_file.replace("_000.vcf", "_var.vcf")

	
	# MAP ASSEMBLY TO REFERENCE (CREATE SAM FILE)
	LOG.info("MAPPING TO REFERENCE GENOME...")
	try:
		num_mm = str(args.num_mismatch)
		len_seed = str(args.length_seed)

		LOG.info(f'\tNumber of mismatches allowed:\t{num_mm}')
		LOG.info(f'\tLength of seed sequence:     \t{len_seed}')

		if len(read_arr) == 2:
			pair1 = read_arr[0]
			pair2 = read_arr[1]

			command = [
			"bowtie2", "-x", REFDIR.path, 
			"-1", pair1.path, "-2", pair2.path,
			"--un-conc", fq_file, 
			"-S", sam_file,
			"-N", num_mm,
			"-L", len_seed]
		else:
			reads = read_arr[0]
			command = [
			"bowtie2", "-x", REFDIR.path, "-U", 
			reads.path, 
			"-S", sam_file,
			"-N", num_mm,
			"-L", len_seed]	
		print(command)
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)



	# CONVERT SAM --> BAM
	try:
		LOG.info("CONVERTING SAM FILE TO BAM FILE...")
		command = ['samtools', 'view', sam_file, '-u', '-o', bam_file]
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
			'samtools', 'mpileup', '-u', '-t', 'DP', '-t', 'AD', '-t', 'ADF', '-t', 'ADR', '-f',
			reference.path, sorted_bam_file, '-o', mpileup_file]
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)


	# CONVERT MPILEUP --> VCF
	try:
		LOG.info("CONVERTING MPILEUP FILE TO VCF FILE...")
		command = ['bcftools', 'call', '-c', mpileup_file, '-o', vcf_file]
		run_subprocess(args, command)
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)


	# FILTER VCF
	try:
		LOG.info("FILTERING VCF FILE...")
		script = "/Strong/proj/shared_code/filter_cf_file_4x_depth.pl"
		command = ['perl', script, '-i', vcf_file, '-o', filtered_vcf_file]

		# Run subprocess
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
		f = open(args.reference_sequence, 'r')
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
	
	# COUNTING VARs
	'''
	count = 0
	try:
		LOG.info("COUNTING VARs...")
		f = open(var_vcf_file, 'r')
		for l in f:
			if l.startswith('#'):
				continue
			else:
				count = count + 1
		f.close()
		LOG.info(f"\tNUMBER OF VARIANTS: {numVar}")
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)
	'''
	# COUNTING VARs
	'''
	TODO: Update command to ignore lines containing comments
	try:
		LOG.info("COUNTING VARs...")

		command = ['wc', '-l', var_vcf_file]

		# GET NUMBER OF VARs
		(stdout, stderr) = subprocess.Popen(command, stdout=PIPE, stderr=PIPE).communicate()
		numVar = str(stdout).split()[0].replace("b'", "")

		LOG.info(f"\tNUMBER OF VARIANTS: {numVar}")
		LOG.info('... DONE')
	except Exception as e:
		LOG.warning(f"... FAILED : {e}")
		exit(e)
	'''



def return_align_rate(args):
	''' Extracts and returns the Bowtie reported Overall Alignment Rate '''

	command = ['grep', 'alignment', LOG_File]

	(stdout, stderr) = subprocess.Popen(command, stdout=PIPE, stderr=PIPE).communicate()
	alignRate = str(stdout).replace("b'", "")
	alignRate = alignRate.replace('\\t', '\t')
	alignRate = alignRate.replace("\\n'", '\n')

	return alignRate


def main(program):
	cwd = os.getcwd()

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='bar_script', add_help=False)
	parent_parser.add_argument('-l', '--length_seed', default=22, help='Sets length of seed. Smaller values increase sensitivity and runtime', type=int)
	parent_parser.add_argument('-mk_off', '--make_db_off', default=False, action='store_true', help='Disable make Bowtie database function.')
	parent_parser.add_argument('-mod', '--module_load', default=False, action='store_true', help='Enable module load function.')
	parent_parser.add_argument('-n', '--num_mismatch', default=0, help='Number of mismatches allowed in the seed sequence. Setting to 1 increases sensitivity and runtime', choices=['0','1'])
	parent_parser.add_argument('-o', '--output_directory', default='bowtie', help='Prefix of output directory', type=str)
	parent_parser.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parent_parser.add_argument('-ref_dir', '--reference_directory', default=None, help='Reference Directory', type=str)
	parent_parser.add_argument('-ref', '--reference_sequence', default=None, help='Reference sequence', type=str)
	parent_parser.add_argument('-rate', '--report_rate', default=False, action='store_true', help='Return the alignment rate')
	subparsers = parent_parser.add_subparsers(help='sub-command help')

	# PARSER : BOWTIE2_UNPAIRED
	bowUn = subparsers.add_parser('single', help='Align FASTQ to reference alignment', parents=[parent_parser])
	bowUn.add_argument('-samp', '--sample', help='file containing unpaired reads to be aligned (FASTQ)', required=True)

	# PARSER : BOWTIE2_PAIRED
	bowPair = subparsers.add_parser('paired', help='Align FASTQ to reference alignment', parents=[parent_parser])
	bowPair.add_argument('-1', '--pair1', help='fwd trimmed paired-end reads to assemble', required=True)
	bowPair.add_argument('-2', '--pair2', help='rev trimmed paired-end reads to assemble', required=True)

	args = parent_parser.parse_args()

	# DETERMINE REFERENCE
	'''
	if args.make_db_off == False:
		if args.reference_sequence != None:
			reference = Fasta(args.reference_sequence)
			ref_id = reference.filename.split('.')[0]
		else:
			args.make_db_off = True
	if args.make_db_off == True:
		if args.reference_directory != None:
			ref_dir = Dir(args.reference_directory)
			ref_id = ref_dir.dirnames
		else:
			print('ERROR. Specify either a reference sequence or directory.')
			exit()	

	'''
	reference = Fasta(args.reference_sequence)
	ref_id = reference.filename.split('.')[0]


	# HANDLE PAIRED READS
	if program == 'paired':
		pair1 = Fastq(args.pair1)
		pair2 = Fastq(args.pair2)
		samp_id = pair1.filename.split('-pair1')[0]

		# CONFIGURE
		basename = configure(args, ref_id, samp_id)
		LOG.info(f'REFERENCE : {reference.filename}')
		LOG.info(f'SAMPLE FWD: {pair1.filename}')
		LOG.info(f'SAMPLE RVS: {pair2.filename}\n')

		# COLLECT SAMPLES
		sample_arr = [pair1, pair2]

	# HANDLE UNPAURED READS
	if program == 'single':
		reads = Fastq(args.sample)
		samp_id = reads.filename.split(".")[0]

		# CONFIGURE
		basename = configure(args, ref_id, samp_id)
		LOG.info(f'REFERENCE: {reference.filename}')
		LOG.info(f'SAMPLE   : {reads.filename}\n')

		# COLLECT SAMPLES
		sample_arr = [reads]


	# SET UP REFERENCE
	# TODO: allow the option to turn off make db for each.
	make_bowtie_db(args, reference, ref_id)


	# RUN ALIGNMENT PIPELINE
	run_pipeline(args, basename, reference, sample_arr)

	if args.report_rate == True:
		align_rate = return_align_rate(args)
		LOG.info(f'ALIGNMENT RATE: {align_rate}')
	


if __name__== "__main__":
	main(sys.argv[1])










