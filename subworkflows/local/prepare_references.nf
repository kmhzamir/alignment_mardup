//
// Prepare reference files
//

include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_MT_SHIFT            } from '../../modules/nf-core/bwamem2/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../modules/nf-core/cat/cat/main'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT               } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD_MT_SHIFT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT                } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_PREPROCESSINTERVALS as GATK_PREPROCESS_WGS   } from '../../modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_PREPROCESSINTERVALS as GATK_PREPROCESS_WES   } from '../../modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_SHIFTFASTA as GATK_SHIFTFASTA                } from '../../modules/nf-core/gatk4/shiftfasta/main'
include { GET_CHROM_SIZES                                    } from '../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_EXTRACT_MT              } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MT_SHIFT          } from '../../modules/nf-core/samtools/faidx/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_GENOME      } from '../../modules/nf-core/sentieon/bwaindex/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_MT_SHIFT    } from '../../modules/nf-core/sentieon/bwaindex/main'
include { TABIX_BGZIPTABIX as TABIX_PBT                      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_DBSNP                         } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PT                            } from '../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {
    take:
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_target_bed      // channel: [mandatory for WES] [ path(bed) ]


    main:
        ch_versions    = Channel.empty()
        ch_tbi         = Channel.empty()
        ch_bgzip_tbi   = Channel.empty()
        ch_bwa         = Channel.empty()
        ch_sentieonbwa = Channel.empty()

        // Genome indices
        BWA_INDEX_GENOME(ch_genome_fasta).index.set{ch_bwa}
        BWAMEM2_INDEX_GENOME(ch_genome_fasta)
        SENTIEON_BWAINDEX_GENOME(ch_genome_fasta).index.set{ch_sentieonbwa}
        SAMTOOLS_FAIDX_GENOME(ch_genome_fasta, [[],[]])
        GATK_SD(ch_genome_fasta)
        ch_fai = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        GET_CHROM_SIZES( ch_fai )
        ch_genome_fasta.map { meta, fasta -> return [meta, fasta, [], [] ] }
            .set {ch_rtgformat_in}

        // MT indices
        ch_fai = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        SAMTOOLS_EXTRACT_MT(ch_genome_fasta, ch_fai)

        //bed indices
        TABIX_PT(ch_target_bed).tbi.set { ch_tbi }
        TABIX_PBT(ch_target_bed).gz_tbi.set { ch_bgzip_tbi }

        // Generate bait and target intervals
        GATK_BILT(ch_target_bed, GATK_SD.out.dict).interval_list
        GATK_ILT(GATK_BILT.out.interval_list)
        GATK_ILT.out.interval_list
            .collect{ it[1] }
            .map { it ->
                meta = it[0].toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                return [[id:meta], it]
            }
            .set { ch_bait_intervals_cat_in }

        // Gather versions
        ch_versions = ch_versions.mix(BWA_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SENTIEON_BWAINDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(GATK_SD.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(GATK_BILT.out.versions)

    emit:
        genome_bwa_index      = Channel.empty().mix(ch_bwa, ch_sentieonbwa).collect()            // channel: [ val(meta), path(index) ]
        genome_bwamem2_index  = BWAMEM2_INDEX_GENOME.out.index.collect()                         // channel: [ val(meta), path(index) ]
        genome_chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()                              // channel: [ path(sizes) ]
        genome_fai            = ch_fai                                                           // channel: [ val(meta), path(fai) ]
        genome_dict           = GATK_SD.out.dict.collect()                                       // channel: [ path(dict) ]
        target_bed            = Channel.empty().mix(ch_tbi, ch_bgzip_tbi).collect()              // channel: [ val(meta), path(bed), path(tbi) ]
        target_intervals      = GATK_BILT.out.interval_list.map{ meta, inter -> inter}.collect() // channel: [ path(interval_list) ]
        versions              = ch_versions                                                      // channel: [ path(versions.yml) ]

}
