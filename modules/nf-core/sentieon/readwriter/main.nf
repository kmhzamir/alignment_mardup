process SENTIEON_READWRITER {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sentieon:202308.01--h43eeafb_0' :
        'biocontainers/sentieon:202308.01--h43eeafb_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.${format}"),                     emit: output
    tuple val(meta), path("*.${index}") ,                     emit: index
    tuple val(meta), path("*.${format}"), path("*.${index}"), emit: output_index
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The following code sets LD_LIBRARY_PATH in the script-section when the module is run by Singularity.
    // That turned out to be one way of overcoming the following issue with the Singularity-Sentieon-containers from galaxy, Sentieon (LD_LIBRARY_PATH) and the way Nextflow runs Singularity-containers.
    // The galaxy container uses a runscript which is responsible for setting LD_PRELOAD properly. Nextflow executes singularity containers using `singularity exec`, which avoids the run script, leading to the LD_LIBRARY_PATH/libstdc++.so.6 error.
    if (workflow.containerEngine == 'singularity') {
        fix_ld_library_path = 'LD_LIBRARY_PATH=/usr/local/lib/:\$LD_LIBRARY_PATH;export LD_LIBRARY_PATH'
    } else {
        fix_ld_library_path = ''
    }

    def args            = task.ext.args   ?: ''
    def args2           = task.ext.args2  ?: ''
    def input_str       = input.sort().collect{"-i $it"}.join(' ')
    def reference       = fasta ? "-r $fasta" : ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    format              = input.extension == "bam" ? "bam"     : "cram"
    index               = format          == "bam" ? "bam.bai" : "cram.crai"
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''
    """
    if [ "\${#SENTIEON_LICENSE_BASE64}" -lt "1500" ]; then  # If the string SENTIEON_LICENSE_BASE64 is short, then it is an encrypted url.
        export SENTIEON_LICENSE=\$(echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d)
    else  # Localhost license file
        # The license file is stored as a nextflow variable like, for instance, this:
        # nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
        export SENTIEON_LICENSE=\$(mktemp)
        echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d > \$SENTIEON_LICENSE
    fi

    if  [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_auth_data_base64} ]; then
        # If sentieon_auth_mech_base64 and sentieon_auth_data_base64 are non-empty strings, then Sentieon is mostly likely being run with some test-license.
        export SENTIEON_AUTH_MECH=\$(echo -n "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_AUTH_DATA=\$(echo -n "${sentieon_auth_data_base64}" | base64 -d)
        echo "Decoded and exported Sentieon test-license system environment variables"
    fi

    $fix_ld_library_path

    sentieon \\
        driver \\
        -t $task.cpus \\
        $reference \\
        $args \\
        $input_str \\
        --algo ReadWriter \\
        $args2 \\
        ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    // The following code sets LD_LIBRARY_PATH in the script-section when the module is run by Singularity.
    // That turned out to be one way of overcoming the following issue with the Singularity-Sentieon-containers from galaxy, Sentieon (LD_LIBRARY_PATH) and the way Nextflow runs Singularity-containers.
    // The galaxy container uses a runscript which is responsible for setting LD_PRELOAD properly. Nextflow executes singularity containers using `singularity exec`, which avoids the run script, leading to the LD_LIBRARY_PATH/libstdc++.so.6 error.
    if (workflow.containerEngine == 'singularity') {
        fix_ld_library_path = 'LD_LIBRARY_PATH=/usr/local/lib/:\$LD_LIBRARY_PATH;export LD_LIBRARY_PATH'
    } else {
        fix_ld_library_path = ''
    }

    def prefix = task.ext.prefix ?: "${meta.id}"
    format = input.extension == "bam" ? "bam"     : "cram"
    index  = format          == "bam" ? "bam.bai" : "cram.crai"
    """
    $fix_ld_library_path

    touch ${prefix}.${format}
    touch ${prefix}.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
