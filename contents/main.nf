/*
 * Copyright (c) 2014-2021 Genome Research Ltd
 *
 * Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
 *
 * This file is part of cgpPindel.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * 1. The usage of a range of years within a copyright statement contained within
 * this distribution should be interpreted as being equivalent to a list of years
 * including the first and last year specified and all consecutive years between
 * them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
 * 2009, 2011-2012’ should be interpreted as being identical to a statement that
 * reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
 * statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
 * identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
 * 2009, 2010, 2011, 2012’.
 */

nextflow.enable.dsl=2

def helpMessage() {
    log.info """
    Usage:

      nextflow run https://github.com/cancerit/cgpPindel -entry <entry-point> --help

    Available entry points:
      - np_generation
      - pindel_pl

    """.stripIndent()
}

def helpNpMessage() {
    log.info """
    Usage:

    The typical command for running the pipeline is as follows:

      nextflow run cgppindel-np -entry np_generation --bams <file_list.txt> --genomefa <genome.fa> --exclude <csv_contigs> --badloci <highDepth.bed.gz> [--range]

    Mandatory arguments:

      --bams      Text list of input BAM/CRAM files to be used during normal panel creation.
      --genomefa  Genome fasta file, expects co-located '*.fa.fai'.
      --exclude   CSV list of contigs to completely exclude from processing '%' as wildcard.
      --badloci   Regions to skip as tabix indexed bed, traditionally regions of high depth.

    Optional arguments:

      --range     Boolean, presence triggers range based normal panel generation [false]
      --outdir    Folder to write output to [\$PWD/results].
      --help      Show this usage.

    """.stripIndent()
}

def helpPindelMessage() {
    log.info """
    Usage:

      nextflow run https://github.com/cancerit/cgpPindel -entry pindel_pl --help

  Required parameters:
    --outdir      Folder to output result to.
    --genomefa    Path to reference genome file *.fa[.gz]
    --tumour      Tumour BAM/CRAM file (co-located index and bas files)
    --normal      Normal BAM/CRAM file (co-located index and bas files)
    --simrep      Full path to tabix indexed simple/satellite repeats.
    --filter      VCF filter rules file (see FlagVcf.pl for details)
    --genes       Full path to tabix indexed coding gene footprints.
    --unmatched   Full path to tabix indexed gff3 of unmatched normal panel
                   - see pindel_np_from_vcf.pl
  Optional
    --seqtype    Sequencing protocol, expect all input to match [WGS]
    --assembly   Name of assembly in use
                  -  when not available in BAM header SQ line.
    --species    Species
                  -  when not available in BAM header SQ line.
    --exclude    Exclude this list of ref sequences from processing, wildcard '%'
                  - comma separated, e.g. NC_007605,hs37d5,GL%
    --badloci    Tabix indexed BED file of locations to not accept as anchors
                  - e.g. hi-seq depth from UCSC
    --skipgerm   Don't output events with more evidence in normal BAM.
    --softfil    VCF filter rules to be indicated in INFO field as soft flags
    --limit      When defined with '-cpus' internally thread concurrent processes.
                  - requires '-p', specifically for pindel/pin2vcf steps
    --debug      Don't cleanup workarea on completion.
    --apid       Analysis process ID (numeric) - for cgpAnalysisProc header info
                  - not necessary for external use

    """.stripIndent()
}

/*
 * Creates a fake BAM to use as the "tumour"
 */
process create_fake_bam {
    input:
        tuple path('genome.fa'), path('genome.fa.fai')

    output:
        tuple path('empty.bam'), path('empty.bam.bai')

    script:
        """
        set -o pipefail
        (samtools dict genome.fa | cut -f 1-4 && echo -e '@RG\tID:1\tSM:FAKE') | samtools view -S -bo empty.bam -
        samtools index empty.bam
        """
}

process np_pindel {
    input:
        tuple path('genome.fa'), path('genome.fa.fai')
        tuple path('badloci.bed.gz'), path('badloci.bed.gz.tbi')
        val exclude
        // the tuple from last process
        tuple path(fake_bam), path(fake_bai)
        // the input channel
        each wt_bam

    output:
        path 'result/*.vcf.gz'

    // to ensure best cpu use actually run 3 commands
    // can't do this (easily) as workflow due to the way the wrapper code works
    script:
        def applyExclude = exclude != 'NO_EXCLUDE' ? "-e $exclude" : ''
        """
        ## setup the species/assembly variables
        SP=\$(samtools view -H ${wt_bam} | grep -m 1 '^@SQ' | perl -ane 'if(m/SP:([^\t]+)/){ print \$1;} else {print q{UNKNOWN};}')
        AS=\$(samtools view -H ${wt_bam} | grep -m 1 '^@SQ' | perl -ane 'if(m/AS:([^\t]+)/){ print \$1;} else {print q{UNKNOWN};}')
        ## process input 1 (tumour)
        pindel.pl -noflag -o result -sp \${SP} -as \${AS} \
        -r genome.fa \
        -t ${fake_bam} \
        -n ${wt_bam} \
        ${applyExclude} \
        -b badloci.bed.gz \
        -c ${task.cpus} \
        -process input -i 1
        ## process input 2 (normal)
        pindel.pl -noflag -o result -sp \${SP} -as \${AS} \
        -r genome.fa \
        -t ${fake_bam} \
        -n ${wt_bam} \
        ${applyExclude} \
        -b badloci.bed.gz \
        -c ${task.cpus} \
        -process input -i 2
        ## do everything else
        pindel.pl -noflag -o result -sp \${SP} -as \${AS} \
        -r genome.fa \
        -t ${fake_bam} \
        -n ${wt_bam} \
        ${applyExclude} \
        -b badloci.bed.gz \
        -c ${task.cpus}
        """
}

process np_creation {
    input:
        path '*.vcf.gz'
        val outdir
        val range

    output:
        // allow for bed or gff3, need to determine how to handle in params
        path 'pindel_range_np.*.gz*'

    publishDir path: "$outdir", mode: 'move', overwrite: true, enabled: true

    script:
        def applyRange = range ? ' -r ' : ''
        """
        pindel_np_from_vcf.pl -o pindel_range_np -s NORMAL ${applyRange} *.vcf.gz
        """
}

process pindel {
    input:
        tuple path('genome.fa'), path('genome.fa.fai')
        tuple path('badloci.bed.gz'), path('badloci.bed.gz.tbi')
        val exclude
        tuple path('mt.bam'), path('mt.bam.bai')
        tuple path('wt.bam'), path('wt.bam.bai')
        val species
        val assembly
        val seqtype
        val outdir

    output:
        tuple path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: vcf
        tuple path('*_mt.bam'), path('*_mt.bam.md5'), path('*_mt.bam.bai'), emit: mt_out
        tuple path('*_wt.bam'), path('*_wt.bam.md5'), path('*_wt.bam.bai'), emit: wt_out

    publishDir path: "$outdir", mode: 'copy', overwrite: true, enabled: true

    // to ensure best cpu use actually run 3 commands
    // can't do this (easily) as workflow due to the way the wrapper code works
    script:
        def applySpecies = species != 'NO_SPECIES' ? "-sp $species" : ''
        def applyAssembly = assembly != 'NO_ASSEMBLY' ? "-as $assembly" : ''
        def applyExclude = exclude != 'NO_EXCLUDE' ? "-e $exclude" : ''
        """
        pindel.pl -noflag -o result \
        ${applySpecies} \
        ${applyAssembly} \
        ${applyExclude} \
        -r genome.fa \
        -t mt.bam \
        -n wt.bam \
        -b badloci.bed.gz \
        -st ${seqtype} \
        -c ${task.cpus}
        # easier to link the files than use "publishDir saveAs:"
        ln -f result/*.vcf.gz* .
        ln -f result/*_mt.bam* .
        ln -f result/*_wt.bam* .

        """
}

process pindel_flag {
    input:
        val unmatched_ext
        path filter
        tuple path('genes.bed.gz'), path('genes.bed.gz.tbi')
        tuple path("unmatched.${unmatched_ext}.gz"), path("unmatched.${unmatched_ext}.gz.tbi")
        tuple path('simplerepeats.bed.gz'), path('simplerepeats.bed.gz.tbi')
        tuple path('input.vcf.gz'), path('input.vcf.gz.tbi')
        // optional
        path softfil
        val apid
        val outdir

    output:
        tuple path('*.pindel.flagged.vcf.gz'), path('*.pindel.flagged.vcf.gz.tbi')

    publishDir path: "$outdir", mode: 'move', overwrite: true, enabled: true

    script:
        def applySoft = softfil.name != 'NO_FILE' ? "-sr $softfil" : ''
        def applyProcId = apid != 'NO_PROCESS' ? "-p $apid" : ''
        """
        MT_NAME=\$(tabix -H input.vcf.gz | grep '^##SAMPLE=<ID=TUMOUR' | perl -ne 'm/SampleName=([^>]+)/; print \$1;')
        WT_NAME=\$(tabix -H input.vcf.gz | grep '^##SAMPLE=<ID=NORMAL' | perl -ne 'm/SampleName=([^>]+)/; print \$1;')
        FlagVcf.pl -r ${filter} -a genes.bed.gz -u unmatched.${unmatched_ext}.gz -s simplerepeats.bed.gz -i input.vcf.gz -o \${MT_NAME}_vs_\${WT_NAME}.pindel.flagged.vcf ${applySoft} ${applyProcId}
        bgzip -c \${MT_NAME}_vs_\${WT_NAME}.pindel.flagged.vcf > \${MT_NAME}_vs_\${WT_NAME}.pindel.flagged.vcf.gz
        tabix -p vcf \${MT_NAME}_vs_\${WT_NAME}.pindel.flagged.vcf.gz
        """
}

workflow {
    if ( params.help ) {
        helpMessage()
        exit 0
    }
}

workflow subwf_pindel_pl {
    // this is a sub workflow, it has no idea about params or "help"
    take:
        unmatched_ext
        genome
        badloci
        exclude
        mt
        wt
        species
        assembly
        outdir
        filter
        genes
        unmatched
        simrep
        // optional
        softfil
        apid
        seqtype

    main:
        pindel(
            genome,
            badloci,
            exclude,
            mt,
            wt,
            species, assembly, seqtype,
            outdir
        )
        pindel_flag(
            unmatched_ext,
            filter,
            genes,
            unmatched,
            simrep,
            pindel.out.vcf,
            softfil,
            apid,
            outdir
        )

    emit:
        pindel_flag.out
}

workflow pindel_pl {
    // This is the workflow for direct use as a stand-alone
    // Show help message
    if ( params.help ) {
        helpPindelMessage()
        exit 0
    }
    log.info """\
    cgpPindel:pindel_pl - NF PIPELINE
    ==================================
    genomefa : ${params.genomefa}
    exclude  : ${params.exclude}
    badloci  : ${params.badloci}
    outdir   : ${params.outdir}
    @TODO
    """
    .stripIndent()

    main:
        // setup tuples for index inclusion
        mt = tuple file(params.tumour), file("${params.tumour}.bai")
        wt = tuple file(params.normal), file("${params.normal}.bai")
        badloci = tuple file(params.badloci), file("${params.badloci}.tbi")
        genome = tuple file(params.genomefa), file("${params.genomefa}.fai")
        genes = tuple file(params.genes), file("${params.genes}.tbi")
        unmatched = tuple file(params.unmatched), file("${params.unmatched}.tbi")
        simrep = tuple file(params.simrep), file("${params.simrep}.tbi")

        unmatched_ext = params.unmatched.contains('gff3') ? 'gff3' : 'bed'

        subwf_pindel_pl(
            unmatched_ext,
            genome,
            badloci,
            params.exclude,
            mt,
            wt,
            params.species,
            params.assembly,
            params.outdir,
            params.filter,
            genes,
            unmatched,
            simrep,
            params.softfil,
            params.apid,
            params.seqtype
        )
}

workflow subwf_np_gen {
    // this is a sub workflow, it has no idea about params or "help"
    take:
        genome
        badloci
        exclude
        bam_ch
        outdir
        range

    main:
        // process
        create_fake_bam(genome)

        // process
        np_pindel(
            genome,
            badloci,
            exclude,
            create_fake_bam.out,
            bam_ch
        )

        np_creation(np_pindel.out.collect(), outdir, range)
    emit:
        np_creation.out
}

workflow np_generation {
    // This is the workflow for direct use as a stand-alone
    // Show help message
    if ( params.help ) {
        helpNpMessage()
        exit 0
    }
    log.info """\
    cgpPindel:np_generation - NF PIPELINE
    ==================================
    bams     : ${params.bams}
    genomefa : ${params.genomefa}
    exclude  : ${params.exclude}
    badloci  : ${params.badloci}
    outdir   : ${params.outdir}
    """
    .stripIndent()

    main:
        // setup tuples for index inclusion
        badloci = tuple file(params.badloci), file("${params.badloci}.tbi")
        genome = tuple file(params.genomefa), file("${params.genomefa}.fai")

        // Create the input channel for BAM/CRAMs
        Channel
            .fromPath(params.bams)
            .splitText()
            .map{it -> file(it.trim())}
            .set { bam_ch }
        // run the sub workflow
        subwf_np_gen(
            genome,
            badloci,
            params.exclude,
            bam_ch,
            params.outdir,
            params.range
        )
}
