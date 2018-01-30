
configfile : "config.yml"


rule all:
	input:
		[s+".krona" for s in config["SAMPLES"]]
		

rule clean :
	input:
		forward = config["RAW_DIR"] + "/{name}_1.fastq",
		reverse = config["RAW_DIR"] + "/{name}_2.fastq"
	output:
		forward = "{name}.clean_1.fastq",
		reverse = "{name}.clean_2.fastq",
		single   = "{name}.clean.single.fastq"
	log:
		"{name}.clean.log"
	shell:
		"sickle pe -f {input.forward} -r {input.reverse} -t sanger -o {output.forward} -p {output.reverse} -s {output.single} > {log}"


rule unzip:
	input:
		"{filename}.fastq.gz"
	output:
		"{filename}.fastq"
	shell:
		"gzip -kd {input}"
		

rule remove_trim:
	input: 
		"{name}.clean_{sens}.fastq" 
	output: 
		"{name}.trimmed_{sens}.fastq"
	log:
		"{name}.trimmed.log"
	shell:
		"cutadapt -g ^ATCGTCGTCGTAGGCTGCTC {input} -o {output} -e 0.1 > {log}"


rule bwa:
	input:
		forward="{name}.clean_1.fastq",
		reverse="{name}.clean_2.fastq",
	output:
		"{name}.host_mapping.sam"
	threads:
		128
	log:
		"{name}.host_mapping.log"
	shell:
		"bwa mem -t {threads} {config[HUMAN_REFERENCE]} {input.forward} {input.reverse} > {output} 2> {log}"
		

rule sam_to_bam:
	input: 
		"{name}.host_mapping.sam"
	output:
		sort     = temp("{name}.sort.host_mapping.sam"),
		bam      = "{name}.host_mapping.bam"

	shell:
		"samtools sort {input} > {output.sort};"
		"samtools view -b {output.sort} > {output.bam}"

rule filter_unmapped:
	input:
		"{name}.host_mapping.bam"
	output:
		"{name}.without_human.bam"
	shell:
		"samtools view -b -f 4 -o {output} {input}"


rule bam_to_fastq:
	input:
		"{name}.without_human.bam"
	output:
		forward = "{name}.without_human_1.fastq",
		reverse = "{name}.without_human_2.fastq"
	log:
		"{name}.bam_to_fastq.log"
	shell:
		"bedtools bamtofastq -i {input} -fq {output.forward} -fq2 {output.reverse} > {log}"

		
rule spades:
	input:
		forward="{name}.without_human_1.fastq",
		reverse="{name}.without_human_2.fastq"
	output:
		"{name}.spade_out/contigs.fasta"
	params:
		spade_dir = "{name}.spade_out"
	log:
		"{name}.spade.log"
	threads: 128
	shell:
		"spades.py --meta -o {params.spade_dir} -1 {input.forward} -2 {input.reverse} -t {threads} > {log}"


# rule centrifuge:
# 	input:
# 		"{name}.spade_out/contigs.fasta"
# 	output:
# 		report = "{name}.centrifuge.report",
# 		result = "{name}.centrifuge.result"
# 	params:
# 		index = config["CENTRIFUGE_INDEX"]
# 	threads:
# 		128
# 	log: 
# 		"{name}.log"
# 	shell:
# 		"centrifuge -f -x {params.index} -U {input}  -p {threads} -S {output.result} --report-file {output.report}"


rule centrifuge:
	input:
		R1 = "{name}.without_human_1.fastq",
		R2 = "{name}.without_human_2.fastq"
	output:
		report = "{name}.centrifuge.report",
		result = "{name}.centrifuge.result"
	params:
		index = config["CENTRIFUGE_INDEX"]
	threads:
		128
	log: 
		"{name}.log"
	shell:
		"centrifuge -q -x {params.index} -1 {input.R1} -2 {input.R2} -p {threads} -S {output.result} --report-file {output.report}"



rule krona_report:
	input:	
		"{name}.centrifuge.result"
	output:
		"{name}.krona"

	shell:
		"centrifuge-kreport -x {config[CENTRIFUGE_INDEX]} {input} > {output}"
