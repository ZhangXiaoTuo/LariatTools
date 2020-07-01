import re, getopt, sys, datetime

def usage():
    print("""
USAGE: python generate_intron_annotation_from_gff3.py [options]
For example:
  python generate_intron_annotation_from_gff3.py -i gene_annotation.gff3 -o intron_annotation.gff3 -g gene,ncRNA_gene -e exon -p Parent &

options:
-i: the gene annotation with gff3 format.
-o: the intron annotation with gff3 format (default: intron_annotation.gff3).
-g: the gene name in 3rd column of gff3 (default: gene).
-e: the exon name in 3rd column of gff3 (default: exon).
-p: the parent's attribute name in 9rd column of gff3 (default: Parent).
-h: print this page.
""")
    sys.exit()

def obtainParameter():
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:g:e:")
    output_file = "intron_annotation.gff3"
    gene_name = "gene"
    exon_name = "exon"
    attribute_name = "Parent"
    if opts:
        for op, value in opts:
            if op == "-i":
                print("-i:", value.split('/')[-1])
                input_file = value
            elif op == "-g":
                print('-g:', value.split('/')[-1])
                gene_name = value
            elif op == "-e":
                print('-e:', value.split('/')[-1])
                exon_name = value
            elif op == "-p":
                print('-p:', value.split('/')[-1])
                attribute_name = value
            elif op == "-o":
                print('-o:', value.split('/')[-1])
                output_file = value
            elif op == "-h":
                usage()
    else:
        usage()
    return(input_file, output_file, gene_name, exon_name, attribute_name)

############
def combine_list(intron_list):
    for n in range(len(intron_list[-1])):
        intron_list[-1][n] = '='.join(intron_list[-1][n])
    intron_list[-1] = ';'.join(intron_list[-1])
    intron_list = '\t'.join(intron_list) + '\n'
    return(intron_list)

def write_gene(whole_gene, o):
    for m in range(len(whole_gene)):
        if whole_gene[m]:
            o.writelines("###\n")
            for n in range(len(whole_gene[m])):
                out = combine_list(whole_gene[m][n])
                o.writelines(out)

def remove_duplicate_intron_name_per_gene(whole_gene):
    intron_loci = {}
    for a in range(len(whole_gene)):
        for b in range(len(whole_gene[a])):
            locus_list = [whole_gene[a][b][0], whole_gene[a][b][3], whole_gene[a][b][4], whole_gene[a][b][6]]
            locus = '_'.join(locus_list)
            if locus not in intron_loci:
                intron_loci[locus] = whole_gene[a][b][-1]
            else:
                whole_gene[a][b][-1][1:3] = intron_loci[locus][1:3]
    return(whole_gene)

def replace_exon_to_intron(exon):
    exon[2] = exon[2].replace('exon', 'intron')
    exon_attribute = exon[-1]
    exon_attribute = [exon_attribute[0], exon_attribute[1], exon_attribute[5], exon_attribute[6]]
    for a in range(len(exon_attribute)):
        for b in range(len(exon_attribute[a])):
            exon_attribute[a][b] = exon_attribute[a][b].replace('exon', 'intron')
    return(exon)

def get_gene_strand(whole_gene):
    strand_symbol = ''
    for m in range(len(whole_gene)):
        for n in range(len(whole_gene[m])):
            if not strand_symbol:
                strand_symbol = whole_gene[m][n][6]
            elif strand_symbol != whole_gene[m][n][6]:
                print('there are different strand in:', whole_gene)
                exit()
    return(strand_symbol)

def reverse_negative_strand_gene(whole_gene):
    strand_symbol = get_gene_strand(whole_gene)
    if strand_symbol == '-':
        whole_gene = whole_gene[::-1]
        for n in range(len(whole_gene)):
            whole_gene[n] = whole_gene[n][::-1]
    return(whole_gene)

def generate_intron_per_transcript(transcript):
    for n in range(len(transcript) - 1):
        transcript[n][3] = str(int(transcript[n][4]) + 1)
        transcript[n][4] = str(int(transcript[n+1][3]) - 1)
        transcript[n] = replace_exon_to_intron(transcript[n])
    transcript = transcript[:-1]
    return(transcript)

def have_gene(whole_gene):
    if whole_gene:
        return(True)

def write_gene_into_file(whole_gene, o):
    if have_gene(whole_gene):
        whole_gene = reverse_negative_strand_gene(whole_gene)
        for n in range(len(whole_gene)):
            whole_gene[n] = generate_intron_per_transcript(whole_gene[n])
        whole_gene = remove_duplicate_intron_name_per_gene(whole_gene)
        whole_gene = reverse_negative_strand_gene(whole_gene)
        write_gene(whole_gene, o)
        #exit()

def add_transcript_exon(whole_gene, line):
    transcript_num = len(whole_gene) - 1
    whole_gene[transcript_num].append(line)
    return(whole_gene)

def add_transcript(whole_gene, line):
    whole_gene.append([line])
    return(whole_gene)

def different_transcript(whole_gene, line):
    transcript_num = len(whole_gene) - 1
    exon_num = len(whole_gene[transcript_num]) - 1
    gene_last_transcript_name = whole_gene[transcript_num][exon_num][-1][0][1]
    new_transcript_name = line[-1][0][1]
    if gene_last_transcript_name != new_transcript_name:
        return(True)

def is_new_transcript(whole_gene, line):
    if not whole_gene:
        return(True)
    elif different_transcript(whole_gene, line):
        return(True)
    else:
        return(False)

def split_line_last_column(line):
    line[-1] = line[-1].split(';')
    for n in range(len(line[-1])):
        line[-1][n] = line[-1][n].split('=')
    return(line)

def add_exon(whole_gene, line):
    line = split_line_last_column(line)
    if is_new_transcript(whole_gene, line):
        whole_gene = add_transcript(whole_gene, line)
    else:
        whole_gene = add_transcript_exon(whole_gene, line)
    return(whole_gene)

def is_new_exon(line, exon_name):
    if line[2] == exon_name:
        return(True)
    else:
        return(False)

def is_new_gene(line, gene_name):
    new_gene_symbol = False
    for name in gene_name:
        new_gene = new_gene or line[2] == name
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
    

def get_introns_per_gene_from_gff3(input_file, output_file, gene_name, exon_name, attribute_name):
    whole_gene = []
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                if not_annotation_line(line):
                    line = line.strip().split('\t')
                    if is_new_gene(line, gene_name):
                        write_gene_into_file(whole_gene, o)
                        whole_gene = []
                    elif is_new_exon(line, exon_name):
                        whole_gene = add_exon(whole_gene, line)
                else:
                    write_into_file(line, o)
            write_gene_into_file(whole_gene, o)


if __name__ == '__main__':

    input_file, output_file, gene_name, exon_name, attribute_name = obtainParameter()

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)

    gene_name = gene_name.split(',')
    get_introns_per_gene_from_gff3(input_file, output_file, gene_name, exon_name, attribute_name)


    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))

