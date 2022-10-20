//
// Will parse the user provided barcode regex and/or pattern
//

include { CREATE_REGEX } from '../../modules/local/create_regex'

workflow CREATE_REGEX_INFO {
    take:
    barcode_regex
    barcode_pattern
    barcode_lengths
    umi_lengths
    fixed_seqs

    main:
    CREATE_REGEX ( barcode_regex,
                   barcode_pattern,
                   barcode_lengths,
                   umi_lengths,
                   fixed_seqs )
        .regex_pattern_file
        .map { parse_regex_info(it) }
        .set { regex }

    emit:
    regex
}

// Function to get list of [ umi_tools_barcode, barcode_regex, barcode_length, umi_length ]
def parse_regex_info(regex_file) {
    def regex_info = [:]
    file(regex_file).withReader {
        String line
        while ( line = it.readLine() ) {
            String[] split_line
            split_line = line.split('\t')

            if (split_line[0] == 'REGEX') {
                regex_info.regex = split_line[1]

            } else if (split_line[0] == 'UMI_TOOLS') {
                regex_info.umi_tools = split_line[1]

            } else if (split_line[0] == 'BC_LENGTH') {
                regex_info.bc_length = split_line[1]

            } else if (split_line[0] == 'UMI_LENGTH') {
                regex_info.umi_length = split_line[1]
            }
        }
    }

    return regex_info
}
