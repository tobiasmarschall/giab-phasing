
samples = {
	'AJ': ['HG004', 'HG003', 'HG002'],
}
samples2accession = {
	'HG002': 'NA24385',
	'HG003': 'NA24149',
	'HG004': 'NA24143',
}
fam2child = {
	'AJ': 'HG002',
}
sample2fam = {}
for fam, sample_list in samples.items():
	for sample in sample_list:
		sample2fam[sample] = fam

ref_versions = ['GRCh38','hg19']
callset_version = '3.3.2'

chromosomes = {
	'GRCh38': [ 'chr{}'.format(i) for i in range(1,23) ],
	'hg19': [ '{}'.format(i) for i in range(1,23) ],
}


def locate_ref(wildcards):
	m = {
		'hg19': '/MMCI/TM/scratch/ref/hg19/human_g1k_v37.fasta',
		'GRCh38': '/MMCI/TM/scratch/ref/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa',
	}
	return m[wildcards.ref_version]

ruleorder: download_giab_ftp > index_bam

rule master:
	input:
		#expand('ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/{sample}_{ref_version}/{sample}_{ref_version}_{endings}',
		#sample=samples['AJ'], ref_version=ref_versions, endings= ['phased_possorted_bam.bam', 'phased_possorted_bam.bam.bai', 'phased_variants.vcf.gz', 'phased_variants.vcf.gz.tbi']),
		#expand('release/{callset_version}/{ref_version}/{sample}.vcf.gz', callset_version=[callset_version], ref_version=ref_versions, sample=samples['AJ']),
		expand('ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.{ending}', ending=['bam', 'bam.bai', 'bam.md5sum']),
		#expand('10x/{ref_version}/{sample}.vcf.gz', ref_version=ref_versions, sample=samples['AJ']),
		'whatshap-merged/GRCh38/giab-3.3.2/10x/trio/AJ.vcf.gz',
		'whatshap-merged/hg19/giab-3.3.2/10x/trio/AJ.vcf.gz',
		'whatshap-merged/hg19/giab-3.3.2/pacbio/trio-childreads/AJ.vcf.gz',
		'whatshap-merged/hg19/giab-3.3.2/pacbio/single-childonly/AJ.vcf.gz',
		expand('ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg002_gr37_{chromosome}.bam', chromosome=chromosomes['hg19']),
		

rule download_giab_ftp:
	output: 'ftp/{file}'
	log: 'ftp/{file}.wgetlog'
	resources: download=1
	shell: 'wget --output-file={log} -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/{wildcards.file}'


def locate_release(wildcards):
	m = {
		('HG002','GRCh38','3.3.2'): ['ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf' + e for e in ['_triophased.vcf.gz', '_triophased.vcf.gz.tbi', '_noinconsistent.bed']],
		('HG004','GRCh38','3.3.2'): ['ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf' + e for e in ['.vcf.gz', '.vcf.gz.tbi', '_noinconsistent.bed']],
		('HG003','GRCh38','3.3.2'): ['ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf' + e for e in ['.vcf.gz', '.vcf.gz.tbi', '_noinconsistent.bed']],
		('HG002','hg19','3.3.2'): ['ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf' + e for e in ['_triophased.vcf.gz', '_triophased.vcf.gz.tbi', '_noinconsistent.bed']],
		('HG004','hg19','3.3.2'): ['ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh37/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf' + e for e in ['.vcf.gz', '.vcf.gz.tbi', '_noinconsistent.bed']],
		('HG003','hg19','3.3.2'): ['ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf' + e for e in ['.vcf.gz', '.vcf.gz.tbi', '_noinconsistent.bed']],
	}
	return m[(wildcards.sample, wildcards.ref_version, wildcards.callset_version)]


#rule link_release:
	#input: 
		#files=locate_release
	#output:
		#vcf='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.vcf.gz',
		#tbi='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.vcf.gz.tbi',
		#bed='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.callable.bed',
	#shell:
		#'''
		#cd release/{wildcards.callset_version}/{wildcards.ref_version}/sample-wise/full;
		#ln -s ../../../../../{input[0]} {wildcards.sample}.vcf.gz;
		#ln -s ../../../../../{input[1]} {wildcards.sample}.vcf.gz.tbi;
		#ln -s ../../../../../{input[2]} {wildcards.sample}.callable.bed;
		#cd -
		#'''

#rule translate_release:
	#input: 
		#files=locate_release
	#output:
		#vcf='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.vcf.gz',
		#tbi='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.vcf.gz.tbi',
		#bed='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.callable.bed',
	#shell:
		#'''
		#zcat {input[0]} | sed -r '/^##contig/ s/ID=([0-9A-Z]+)/ID=chr\\1/g' | sed '/^[^#]/ s/^/chr/g' | bgzip > {output.vcf};
		#bcftools index --tbi {output.vcf}
		#cat {input[2]} | sed '/^[^#]/ s/^/chr/g' > {output.bed}
		#'''
		
rule translate_release:
	input: 
		files=locate_release
	output:
		vcf='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.vcf.gz',
		tbi='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.vcf.gz.tbi',
		bed='release/{callset_version}/{ref_version}/sample-wise/full/{sample}.callable.bed',
	shell:
		'''
		(bcftools view -h {input[0]} && bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n' {input[0]} ) | bgzip > {output.vcf};
		bcftools index --tbi {output.vcf}
		cp {input[2]} {output.bed}
		'''


rule merge_release:
	input:
		vcfs=lambda wc: ['release/{}/{}/sample-wise/full/{}.vcf.gz'.format(wc.callset_version, wc.ref_version, sample) for sample in samples[wc.family]],
	output:
		vcf='release/{callset_version}/{ref_version}/family-wise/full/{family}.vcf.gz',
		tbi='release/{callset_version}/{ref_version}/family-wise/full/{family}.vcf.gz.tbi',
	shell:
		'bcftools merge --output-type z --missing-to-ref -o {output.vcf} {input.vcfs} && bcftools index --tbi {output.vcf}'


rule intersect_callable_regions:
	input:
		beds=lambda wc: ['release/{}/{}/sample-wise/full/{}.callable.bed'.format(wc.callset_version, wc.ref_version, sample) for sample in samples[wc.family]],
	output:
		bed='release/{callset_version}/{ref_version}/family-wise/full/{family}.callable.bed',
	shell:
		'bedtools intersect -a {input.beds[0]} -b {input.beds[1]} | bedtools intersect -a - -b {input.beds[2]} > {output.bed}'


rule create_sample_name_file:
	output:
		samples_name='sample_names/{sample}'
	shell:
		'echo {wildcards.sample} > {output.samples_name}'


#def locate_10x_vcf(wildcards):
	#return 'ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/{0}_{1}/{0}_{1}_phased_variants.vcf.gz'.format(samples2accession[wildcards.sample],wildcards.ref_version)


rule reheader_10x_GRCh38:
	input:
		vcf=lambda wc: 'ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/{0}_GRCh38/{0}_GRCh38_phased_variants.vcf.gz'.format(samples2accession[wc.sample]),
		samples_name='sample_names/{sample}',
	output:
		vcf='10x/GRCh38/sample-wise/full/{sample}.vcf.gz',
		tbi='10x/GRCh38/sample-wise/full/{sample}.vcf.gz.tbi',
	shell:
		'(bcftools reheader -s {input.samples_name} {input.vcf} | gunzip | awk \'BEGIN {{OFS="\\t"}} $0 ~ /^#/ {{print}} ($0 !~ /^#/) && ($7=="PASS") {{$8="."; print}}\'  | bgzip > {output.vcf} ) && bcftools index --tbi {output.vcf}'

rule reheader_10x_hg19:
	input:
		vcf=lambda wc: 'ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/{0}_hg19/{0}_hg19_phased_variants.vcf.gz'.format(samples2accession[wc.sample]),
		samples_name='sample_names/{sample}',
	output:
		vcf='10x/hg19/sample-wise/full/{sample}.vcf.gz',
		tbi='10x/hg19/sample-wise/full/{sample}.vcf.gz.tbi',
	shell:
		'(bcftools reheader -s {input.samples_name} {input.vcf} | gunzip | awk \'BEGIN {{OFS="\\t"}} $0 ~ /^#/ {{print}} ($0 !~ /^#/) && ($7=="PASS") {{$8="."; print}}\'  | sed -r \'/^##contig/ s/ID=chr/ID=/g\' | sed \'/^[^#]/ s/^chr//g\' | bgzip > {output.vcf} ) && bcftools index --tbi {output.vcf}'

rule vcf_extract_chromosome:
	input:
		vcf='{dir}/{level}/full/{what}.vcf.gz',
		tbi='{dir}/{level}/full/{what}.vcf.gz.tbi',
	output:
		vcf='{dir}/{level,(sample|family)-wise}/by-chromosome/{what}.{chromosome}.vcf.gz',
		tbi='{dir}/{level,(sample|family)-wise}/by-chromosome/{what}.{chromosome}.vcf.gz.tbi',
	shell:
		'bcftools view --output-type z -o {output.vcf} {input.vcf} {wildcards.chromosome} && bcftools index --tbi {output.vcf}'


rule vcf_unphase_chromosome:
	input:
		vcf='{dir}/by-chromosome/{what}.{chromosome}.vcf.gz'
	output:
		vcf='{dir}/by-chromosome-unphased/{what}.{chromosome}.vcf.gz',
		tbi='{dir}/by-chromosome-unphased/{what}.{chromosome}.vcf.gz.tbi',
	shell:
		'~/scm/whatshap.to-run/bin/whatshap unphase {input.vcf} | bgzip >  {output.vcf} && bcftools index --tbi {output.vcf}'


rule vcf_extract_callable:
	input:
		vcf='release/{callset_version}/{ref_version}/family-wise/by-chromosome-unphased/{family}.{chromosome}.vcf.gz',
		bed='release/{callset_version}/{ref_version}/family-wise/full/{family}.callable.bed',
	output:
		vcf='release/{callset_version}/{ref_version}/family-wise/by-chromosome-unphased-callable/{family}.{chromosome}.vcf.gz',
		tbi='release/{callset_version}/{ref_version}/family-wise/by-chromosome-unphased-callable/{family}.{chromosome}.vcf.gz.tbi',
	shell:
		'bedtools intersect -header -a {input.vcf} -b {input.bed} | bgzip >  {output.vcf} && bcftools index --tbi {output.vcf}'


'release/{callset_version}/{ref_version}/family-wise/full/{family}.callable.bed',
rule merge_10x_vcfs:
	input:
		vcfs=lambda wc: ['10x/{}/sample-wise/full/{}.vcf.gz'.format(wc.ref_version, sample) for sample in samples[wc.family]],
	output:
		vcf='10x/{ref_version}/family-wise/full/{family}.vcf.gz',
		tbi='10x/{ref_version}/family-wise/full/{family}.vcf.gz.tbi',
	shell:
		'bcftools merge --output-type z -o {output.vcf} {input.vcfs} && bcftools index --tbi {output.vcf}'


rule create_ped:
	output:
		ped='ped/AJ.ped',
	shell:
		'echo "AJ HG002 HG003 HG004 1 1" > {output.ped}'


rule add_pacbio_readgroup:
	input:
		bam='ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam',
	output:
		bam='pacbio/hg19/HG002.bam',
	log:
		'pacbio/hg19/HG002.log'
	shell:
		'picard AddOrReplaceReadGroups I={input.bam} O={output.bam} PL=PacBio PU=unknown LB=unknown SM=HG002 > {log} 2>&1'


rule index_bam:
	input: 
		bam='{file}.bam'
	output:
		bai='{file}.bam.bai'
	shell:
		'samtools index {input.bam}'


# =========== WHATSHAP ===========
#rule whatshap_single:
	#input:
		#vcf=get_genotype_vcf,
		#phaseinput=get_phase_input,
		#ref=ref
	#output:
		#vcf='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.vcf',
		#corrected_gts='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.genotype-changes.tsv',
		#readlist='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.used-reads.tsv'
	#log: 'whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.vcf.log'
	#shell: '~/scm/whatshap.to-run/bin/whatshap phase --indels --distrust-genotypes --output-read-list {output.readlist} --changed-genotype-list {output.corrected_gts} --reference {input.ref} -o {output.vcf} {input.vcf} {input.phaseinput} > {log} 2>&1'


def get_genotype_vcf(wildcards):
	if wildcards.genotypes.startswith('giab-'):
		callset_version = wildcards.genotypes[5:]
		return 'release/{}/{}/family-wise/by-chromosome-unphased-callable/{}.{}.vcf.gz'.format(callset_version, wildcards.ref_version, wildcards.family, wildcards.chromosome)
	else:
		assert False


def get_phase_input(wildcards):
	if wildcards.phaseinput == '10x':
		return '10x/{}/family-wise/by-chromosome/{}.{}.vcf.gz'.format(wildcards.ref_version, wildcards.family, wildcards.chromosome)
	else:
		assert False


def get_phase_input_child(wildcards):
	if wildcards.phaseinput == 'pacbio':
		return 'pacbio/{}/{}.bam'.format(wildcards.ref_version, fam2child[wildcards.family])
	else:
		assert False


def get_phase_input_child_aux(wildcards):
	if wildcards.phaseinput == 'pacbio':
		return 'pacbio/{}/{}.bam.bai'.format(wildcards.ref_version, fam2child[wildcards.family])
	else:
		assert False


rule whatshap_trio:
	input:
		vcf=get_genotype_vcf,
		phaseinput=get_phase_input,
		ref=locate_ref,
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/{ref_version}/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.vcf.gz',
		recomb='whatshap/{ref_version}/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.recomb',
		readlist='whatshap/{ref_version}/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.used-reads.tsv'
	log: 'whatshap/{ref_version}/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.vcf.log'
	shell: '~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --indels --output-read-list {output.readlist} --reference {input.ref} --recombination-list {output.recomb} {input.vcf} {input.phaseinput}  2> {log} | bgzip > {output.vcf}'


rule whatshap_trio_childreads:
	input:
		vcf=get_genotype_vcf,
		phaseinput=get_phase_input_child,
		phaseinput_aux=get_phase_input_child_aux,
		ref=locate_ref,
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/{ref_version}/{genotypes}/{phaseinput}/trio-childreads/{family}.{chromosome}.vcf.gz',
		readlist='whatshap/{ref_version}/{genotypes}/{phaseinput}/trio-childreads/{family}.{chromosome}.used-reads.tsv',
	log: 'whatshap/{ref_version}/{genotypes}/{phaseinput}/trio-childreads/{family}.{chromosome}.vcf.log'
	shell: '~/scm/whatshap.to-run/bin/whatshap phase --max-coverage 45 --chromosome {wildcards.chromosome} --ped {input.ped} --indels --output-read-list {output.readlist} --reference {input.ref} {input.vcf} {input.phaseinput} 2> {log} | bgzip > {output.vcf}'


rule whatshap_child_only:
	input:
		vcf=get_genotype_vcf,
		phaseinput=get_phase_input_child,
		phaseinput_aux=get_phase_input_child_aux,
		ref=locate_ref
	output:
		vcf='whatshap/{ref_version}/{genotypes}/{phaseinput}/single-childonly/{family}.{chromosome}.vcf.gz',
		readlist='whatshap/{ref_version}/{genotypes}/{phaseinput}/single-childonly/{family}.{chromosome}.used-reads.tsv',
	params:
		child=lambda wc: fam2child[wc.family]
	log: 'whatshap/{ref_version}/{genotypes}/{phaseinput}/single-childonly/{family}.{chromosome}.vcf.log'
	shell: '~/scm/whatshap.to-run/bin/whatshap phase --chromosome {wildcards.chromosome} --sample {params.child} --indels --output-read-list {output.readlist} --reference {input.ref} {input.vcf} {input.phaseinput} 2> {log} | bgzip > {output.vcf}'


rule concat_vcfs:
	input: 
		vcfs=lambda wc: ['whatshap/{}/{}/{}/{}/{}.{}.vcf.gz'.format(wc.ref_version,wc.genotypes, wc.phaseinput, wc.what, wc.family, chromosome) for chromosome in chromosomes[wc.ref_version]]
	output:
		vcf='whatshap-merged/{ref_version}/{genotypes}/{phaseinput}/{what}/{family}.vcf.gz',
		tbi='whatshap-merged/{ref_version}/{genotypes}/{phaseinput}/{what}/{family}.vcf.gz.tbi',
	shell:
		'bcftools concat {input.vcfs} | bgzip > {output.vcf} && bcftools index --tbi {output.vcf}'
