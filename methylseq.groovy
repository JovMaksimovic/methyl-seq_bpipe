// Bpipe for analysing bisulfite sequencing methyl-seq data
// Written by: Jovana Maksimovic
// Last modified on: 21/01/2015

REF = "/mnt/storage/shared/genomes/hg19/fasta/single_fasta_file/" // Location of bisulfite converted genome

def get_name(filename) {

    name = filename.split("_")[0]
    num = filename.split("_")[1]

    return(name + "_R" + num + ".fastq.gz")

}

sra_to_fastq = {
  // Convert SRA files to fastq format 
  transform(".sra") to ("_1.fastq.gz","_2.fastq.gz"){
    exec """
      fastq-dump --split-3 --gzip $input.sra;
    """
  }
}

fastqc = {
    // Run FASTQC to generate QC metrics for the reads  
    output.dir = "fastqc"
    transform(".fastq.gz",".fastq.gz") to ("_fastqc.zip","_fastqc.zip"){
        exec """
          fastqc --quiet -o ${output.dir} $input1.gz $input2.gz
        """
    }
}

trim_reads = {
  // Trim "clip5" bp from 5' end of both read pairs
  // Trim "clip3" bp from 3' end of both read pairs adapter/quality trimming has been performed
  clip5 = "10"
  clip3 = "5"

  transform("(.*)_1.fastq.gz","(.*)_2.fastq.gz") to ("\$1_1_val_1.fq.gz","\$1_2_val_2.fq.gz") {
        exec """
            trim_galore --paired --clip_R1 ${clip5} --clip_R2 ${clip5} --three_prime_clip_R1 ${clip3} --three_prime_clip_R2 ${clip3}
	    --fastqc_args "--quiet --outdir ./fastqc" $input1.gz $input2.gz
        """
  }
}

rename_fastq = {
  // Rename fastq files; this helps next stage tu run correctly
  newname1 = get_name(input1)
  newname2 = get_name(input2)

  transform("(.*)_1_val_1.fq.gz","(.*)_2_val_2.fq.gz") to ("\$1_R1.fastq.gz","\$1_R2.fastq.gz"){
    exec """
      cp $input1.gz $newname1;
      cp $input2.gz $newname2;
    """
  }
}

map_reads = {
  // Map read pairs using bismark and bowtie2
  transform("gz_bismark_bt2_pe.bam") {
    exec """
      bismark --bowtie2 --bam -n 1 $REF -1 $input1.gz -2 $input2.gz
    """
  }
}

dedupe = {
  // Remov duplicates
  transform("bam") to ("deduplicated.bam"){
    exec """
      deduplicate_bismark --paired --bam $input.bam
    """
  }
}

extract_methylation = {
  // Get methylation calls at each CpG
  output.dir = "methylation_calls"
  transform("bam") to ("CpG_report.txt"){
    exec """
      bismark_methylation_extractor -p --no_overlap --merge_non_CpG --no_header --bedGraph --cytosine_report --genome_folder $REF
      -o ${output.dir} $input.bam;
    """
  }
}

get_covered_CpGs = {
  // Create files containing only CpGs that are covered by 1+ reads
  filter("covered") {
    exec """
      awk '\$1 ~ /chr[12]?[0-9]\$/ && (\$4 + \$5 != 0) {printf(\"%s\\t%s\\t%s\\t%s\\t%s\\n\",\$1,\$2,\$3,\$4+\$5,\$4)}'
      $input.txt > $output.txt
    """
  }
}

sort_bam = {
  // Sort bam file
  filter("sorted") {
         exec "samtools sort $input.bam $output.prefix"
  }
}

index_bam = {
  // Index bam file
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

run {
  // Run the methyl-seq pipeline
  "%.sra" * [sra_to_fastq + fastqc + trim_reads + rename_fastq + map_reads] +
  "%.bam" * [dedupe + [extract_methylation + get_covered_CpGs, sort_bam + index_bam]]
}



