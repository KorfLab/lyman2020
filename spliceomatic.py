import random

def isoformer(don, acc, samples=1000, min_intron=35, min_exon=15):

	# generate transcripts
	isoforms = {}
	introns = {}
	itotal = 0
	for n in range(samples):
		tds = [i for i in don if random.random() < don[i]]
		tas = [i for i in acc if random.random() < acc[i]]
		used = {}
		splice = []
		last_acc = 0
		for d in tds:
			edist = d - last_acc
			if edist < min_exon: continue
			for a in tas:
				idist = a - d;
				if idist >= min_intron and a not in used and d not in used:
					intron = (d, a)
					itotal += 1
					if intron not in introns: introns[intron] = 1
					else:                     introns[intron] += 1
					splice.append(intron)
					used[d] = True
					used[a] = True
					last_acc = a
		
		t = tuple(splice)
		if t not in isoforms: isoforms[t] = 1
		else:                 isoforms[t] += 1

	for i in introns: introns[i] /= itotal

	# get site usage
	du, au = {}, {}
	for iso in isoforms:
		count = isoforms[iso]
		for intron in iso:
			d, a = intron
			if d not in du: du[d] = count
			else:           du[d] += count
			if a not in au: au[a] = count
			else:           au[a] += count
	for v in du: du[v] /= samples
	for v in au: au[v] /= samples

	return isoforms, introns, du, au


