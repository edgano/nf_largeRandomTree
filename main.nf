#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */

 // input sequences to align in fasta format
params.seqs = "/users/cn/lmansouri/PROJECTS/BENCHFAM_28.0_2018/BENCHFAM_FOR_REGRESSIVE/FILES_FOR_REGRESSIVE/combined_seqs/PF13193.fa"

// output directory
params.output = "${baseDir}/results"

log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Output directory (DIRECTORY)                   : ${params.output}
         """
         .stripIndent()

// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  .view()
  .into { seqsCh; seqs2 }
}

process generateTree {
    conda './environment.yml'

    tag "${id}"
    publishDir "${params.output}/trees", mode: 'copy', overwrite: true

    input:
      set val(id), file(seqs) from seqsCh

    output:
     set val(id), file("${id}.rndTree.dnd") into reformatSeqsOut

    script:
    """
    #!/usr/bin/env Rscript
    fastaFile <- seqinr::read.fasta(file = "$seqs")

    #how many fasta sequences
    length(fastaFile)

    #get the seqs IDs
    seqNames<-seqinr::getName(fastaFile)

    #define numbers of leave
    n <- length(fastaFile)

    #produce the tree
    rt<-ape::rmtree(1, n, rooted = TRUE, tip.label = seqNames, br = runif)

    #convert & save to newick
    ape::write.tree(rt, file = "${id}.rndTree.dnd", append = FALSE,digits = 2, tree.names = FALSE)

    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
