import re
import getopt
import sys
import datetime

col_seq, att_seq, att_ass = '\t', ';', '='
att_transcript_name = "Parent"
att_transcript_seq = ":"
att_rank_key = 'rank'
att_exon_key = 'exon_id'
att_exon_alias = 'exon_id,Name'
freature_gene_name = "gene,ncRNA_gene"
freature_exon_name = "exon"
replace_exon_key = "exon"
replace_intron_key = "intron"
output_order = "Parent,Name,constitutive,ensembl_end_phase,ensembl_phase,exon_id,rank"

def usage():
    print("""
 USAGE: python generate_intron_annotation_from_gff3.py [options]
 For example:
  python generate_intron_annotation_from_gff3.py -i gene_annotation.gff3 -o intron_annotation.gff3 -g gene,ncRNA_gene -e exon -p Parent &

 options:
 -i: the gene annotation with gff3 format.
 -o: the intron annotation with gff3 format (default: intron_annotation.gff3).
 -g: the gene name of gff3 in 3rd column (default: gene).
 -e: the exon name of gff3 in 3rd column (default: exon).
 -p: the parent's attribute name in 9rd column of gff3 (default: Parent).
 -h: print this page.
 """)
    sys.exit()


def obtainParameter():
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:g:e:p:")
    output_file = "intron_annotation.gff3"

    if opts:
        for op, value in opts:
            if op == "-i":
                input_file = value
            elif op == "-o":
                output_file = value
            elif op == "-g":
                freature_gene_name = value
            elif op == "-e":
                freature_exon_name = value
            elif op == "-p":
                att_transcript_name = value
            elif op == "-h":
                usage()
    else:
        usage()
    return(input_file, output_file)

############


def covert_str(exon):
    new_atts = ''
    for att in output_order:
        if att in exon[-1]:
            new_atts = new_atts + att + att_ass + exon[-1][att] + att_seq
    output_str = col_seq.join(exon[:-1]) + col_seq + new_atts + '\n'
    return(output_str)


def write_gene(whole_gene, o):
    for transcript in whole_gene:
        if transcript:
            o.writelines("###\n")
            for exon in transcript:
                out = covert_str(exon)
                out = out.replace(replace_exon_key, replace_intron_key)
                o.writelines(out)
        else:
            pass

def remove_duplicate_minus_strand_intron_name(whole_gene):
    intron_loci = {}
    gene_len = len(whole_gene)
    for m in range(gene_len):
        tran_len = len(whole_gene[gene_len-m-1])
        for n in range(tran_len):
            exon = whole_gene[gene_len-m-1][tran_len-n-1]
            intron_loci_str = '_'.join([exon[0], exon[3], exon[4], exon[6]])
            if intron_loci_str not in intron_loci:
                intron_loci[intron_loci_str] = exon[-1][att_exon_key]
            else:
                for name in att_exon_alias:
                    whole_gene[gene_len-m-1][tran_len-n-1][-1][name] = intron_loci[intron_loci_str]
    return(whole_gene)

def remove_duplicate_plus_strand_intron_name(whole_gene):
    intron_loci = {}
    for m in range(len(whole_gene)):
        for n in range(len(whole_gene[m])):
            exon = whole_gene[m][n]
            intron_loci_str = '_'.join([exon[0], exon[3], exon[4], exon[6]])
            if intron_loci_str not in intron_loci:
                intron_loci[intron_loci_str] = exon[-1][att_exon_key]
            else:
                for name in att_exon_alias:
                    whole_gene[m][n][-1][name] = intron_loci[intron_loci_str]
    return(whole_gene)

def remove_duplicate_intron_name(whole_gene):
    intron_loci = {}
    if len(whole_gene[0])>0:
        print(whole_gene)
        strand_symbol = whole_gene[0][0][6]
        if strand_symbol == '+':
            whole_gene = remove_duplicate_plus_strand_intron_name(whole_gene)
        elif strand_symbol == '-':
            whole_gene = remove_duplicate_minus_strand_intron_name(whole_gene)
        else:
            pass
    else:
        pass
    return(whole_gene)


def generate_intron_name(whole_gene):
    for m in range(len(whole_gene)):
        for n in range(len(whole_gene[m])):
            exon = whole_gene[m][n]
            intron_name = exon[-1][att_transcript_name].split(att_transcript_seq)[1] + replace_intron_key + exon[-1][att_rank_key]
            whole_gene[m][n][-1][att_exon_key] = intron_name
            for name in att_exon_alias:
                whole_gene[m][n][-1][name] = whole_gene[m][n][-1][att_exon_key]
    return(whole_gene)


def generate_minus_strand_intron(transcript):
    tran_len = len(transcript)
    for n in range(1, tran_len):
        transcript[-n][4] = str(int(transcript[-n][3]) - 1)
        transcript[-n][3] = str(int(transcript[-n-1][4]) + 1)
    transcript = transcript[1:]
    return(transcript)

def generate_plus_strand_intron(transcript):
    for n in range(len(transcript) - 1):
        transcript[n][3] = str(int(transcript[n][4]) + 1)
        transcript[n][4] = str(int(transcript[n+1][3]) - 1)
    transcript = transcript[:-1]
    return(transcript)

def get_transcript_strand(transcript):
    strand_symbol = ''
    for exon in transcript:
        if not strand_symbol:
            strand_symbol = exon[6]
        elif strand_symbol != exon[6]:
            print('there are different strand in:', transcript)
            exit()
        else:
            pass
    return(strand_symbol)


def generate_intron_per_transcript(transcript):
    strand_symbol = get_transcript_strand(transcript)
    if strand_symbol == '+':
        transcript = generate_plus_strand_intron(transcript)
    elif strand_symbol == '-':
        transcript = generate_minus_strand_intron(transcript)
    else:
        pass
    return(transcript)


def have_gene(whole_gene):
    if whole_gene:
        return(True)

def write_gene_into_file(whole_gene, o):
    if have_gene(whole_gene):
        for n in range(len(whole_gene)):
            whole_gene[n] = generate_intron_per_transcript(whole_gene[n])
        whole_gene = generate_intron_name(whole_gene)
        whole_gene = remove_duplicate_intron_name(whole_gene)
        write_gene(whole_gene, o)
    else:
        pass


def add_transcript_exon(whole_gene, line):
    whole_gene[-1].append(line)
    return(whole_gene)

def add_transcript(whole_gene, line):
    whole_gene.append([line])
    return(whole_gene)

def different_transcript(whole_gene, line):
    gene_last_transcript_name = whole_gene[-1][-1][-1][att_transcript_name]
    new_transcript_name = line[-1][att_transcript_name]
    if gene_last_transcript_name != new_transcript_name:
        return(True)
    else:
        return(False)


def is_new_transcript(whole_gene, line):
    if not whole_gene:
        return(True)
    elif different_transcript(whole_gene, line):
        return(True)
    else:
        return(False)


def split_line_last_column(line):
    line[-1] = line[-1].split(att_seq)
    per_att_dict = {}
    for n in range(len(line[-1])):
        tmp_atts = line[-1][n].split(att_ass)
        per_att_dict[tmp_atts[0]] = tmp_atts[1]
    line[-1] = per_att_dict
    return(line)


def add_exon(whole_gene, line):
    line = split_line_last_column(line)
    if is_new_transcript(whole_gene, line):
        whole_gene = add_transcript(whole_gene, line)
    else:
        whole_gene = add_transcript_exon(whole_gene, line)
    return(whole_gene)


def is_new_exon(line, freature_exon_name):
    if line[2] in freature_exon_name:
        return(True)
    else:
        return(False)


def is_new_gene(line, freature_gene_name):
    new_gene_symbol = False
    if line[2] in freature_gene_name:
        new_gene_symbol = True
    return(new_gene_symbol)


def write_into_file(line, o):
    residue = line.strip().replace('#', '').replace(' ', '')
    if residue:
        o.writelines(line)


def not_annotation_line(line):
    if line[0] != '@' and line[0] != '#':
        return(True)
    else:
        return(False)


def get_introns_per_gene_from_gff3(input_file, output_file):
    whole_gene = []
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                if not_annotation_line(line):
                    line = line.strip().split(col_seq)
                    if is_new_gene(line, freature_gene_name):
                        write_gene_into_file(whole_gene, o)
                        whole_gene = []
                    elif is_new_exon(line, freature_exon_name):
                        whole_gene = add_exon(whole_gene, line)
                    else:
                        pass
                else:
                    write_into_file(line, o)
            write_gene_into_file(whole_gene, o)


if __name__ == '__main__':

    input_file, output_file = obtainParameter()

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)

    freature_gene_name = freature_gene_name.split(',')
    freature_exon_name = freature_exon_name.split(',')
    output_order = output_order.split(',')
    att_exon_alias = att_exon_alias.split(',')
    get_introns_per_gene_from_gff3(input_file, output_file)

    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))
