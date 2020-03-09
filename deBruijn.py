class Node(object):
	def __init__(self, kmer, ineo=None, exeo=None, count=1, branching=False):
		self.kmer = kmer
		self.ineo = ineo if ineo is not None else list() # list of strings!
		self.exeo = exeo if ineo is not None else list() # "
		self.count = count # number of occurences in the ENTIRE file
		self.branching = branching
	def __eq__(self, other):
		return True if self.kmer == other.kmer else False
	def __hash__(self):
		return hash((self.kmer))
	def __str__(self):
		return self.kmer

class deBruijn(object):
	def __init__(self):
		self.nodes = dict()
	def add(self, node):
		if node.kmer in self.nodes:
			self.nodes[node.kmer].count += 1
		else:
			self.nodes[node.kmer] = node

	def is_branching(self, node):
		if len(self.nodes[node.kmer].ineo) > 1 or len(self.nodes[node.kmer].exeo) > 1:
			return True
		return False

	def update(self, A, B):
		if B.kmer not in self.nodes[A.kmer].exeo:
			self.nodes[A.kmer].exeo.append(B.kmer)
		if B.kmer not in self.nodes:
			self.add(B)
		else:
			self.nodes[B.kmer].count += 1
			if A.kmer not in self.nodes[B.kmer].ineo:
				self.nodes[B.kmer].ineo.append(A.kmer)
			if self.is_branching(A):
				self.nodes[A.kmer].branching = True
			if self.is_branching(B):
				self.nodes[B.kmer].branching = True
			return
		if self.is_branching(A):
			self.nodes[A.kmer].branching = True
		if self.is_branching(B):
			self.nodes[B.kmer].branching = True
		self.nodes[B.kmer].ineo.append(A.kmer)

def assemble(source, dB):
	curr_mer = source # string only
	contig = source
	while dB.nodes[curr_mer].branching == False and dB.nodes[curr_mer].exeo: # branching and sinks
		curr_mer = dB.nodes[curr_mer].exeo[0] # should only have one!
		contig += curr_mer[-1]
	return contig

def remove_branching_nodes(dB):
	bad_guys = []
	for kmer in dB.nodes:
		if dB.nodes[kmer].branching:
			for ine in dB.nodes[kmer].ineo:
				dB.nodes[ine].exeo.remove(kmer)
			for exe in dB.nodes[kmer].exeo:
				dB.nodes[exe].ineo.remove(kmer)
			bad_guys.append(kmer)
	for kmer in bad_guys:
		del dB.nodes[kmer]

def contigs(dB):
	remove_branching_nodes(dB)
	sources = [kmer for kmer in dB.nodes if len(dB.nodes[kmer].ineo) == 0]
	contigs = []
	print([source for source in sources])
	for source in sources:
		contigs.append(assemble(source, dB))
	with open('all_good_contigs80.txt', 'w') as fptr:
		for c in contigs:
			fptr.write(c + '\n')
	return contigs

def filter_contigs(min_size, contigs):
	filtered = [c for c in contigs if len(c) >= min_size]
	with open('output_contigs80.txt', 'w') as fptr, open('contig_lengths80.txt', 'w') as gptr:
		for c in filtered:
			fptr.write(c + '\n')
			gptr.write(str(len(c)) + '\n')

def goodReads(seqs, dB):
	goodies = [seq for seq in seqs]
	for seq in seqs:
		seq_kmers = get_kmers(seq)
		for k in seq_kmers:
			if dB.nodes[k].count == 1 and seq in goodies:
				goodies.remove(seq)
	with open('good_reads80.txt', 'w') as fptr:
		for goodie in goodies:
			fptr.write(goodie + '\n')
	return goodies

def de_bruijn(k, seqs):
	dB = deBruijn()
	for seq in seqs:
		i = 0
		prev = None
		while i <= len(seq)-k:
			curr = Node(seq[i:k+i])
			if prev is not None: # past first kmer
				dB.update(prev, curr)
			else:
				dB.add(curr)
			prev = curr
			i += 1
	return dB

def get_kmers(seq):
	kmers = []; i = 0
	while i <= len(seq)-k:
		kmers.append(seq[i:k+i])
		i += 1
	return kmers

def read_sequences(file_path):
	with open(file_path, 'r') as fptr:
		seqs = [s.rstrip() for s in fptr.readlines()]
	return seqs

if __name__ == '__main__':
	k = 10
	contig_size = 100
	seqs = read_sequences('sequence_reads.txt')
	dB = de_bruijn(k, seqs)
	goodies = goodReads(seqs, dB)
	good_dB = de_bruijn(k, goodies)
	contigs = contigs(good_dB)
	filter_contigs(contig_size, contigs)
