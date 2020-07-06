import re, getopt, sys, datetime
from sys import path
path.append(r'../')
from base.getParm import obtain_parameters


############
# 需要注意gene为空和transcript为空的情况

def write_gene_into_file(gene_dict, transcript_dict, gene_order, output_file):
    with open(output_file, 'w') as o:
        for gene in gene_order:
            if gene in gene_dict:
                for transcript in gene_dict[gene]:
                    if transcript in transcript_dict:
                        for exon in transcript_dict[transcript]:
                            o.writelines(exon)


def replace_exon_to_intron(exon, replace_strings):
    exon = exon.replace(replace_strings[0], replace_strings[1])
    return(exon)

def join_splited_annotation_line(splited_line, Separator, attribute_assignment_operator, retained_attributes):
    new_line = splited_line[:-1]
    attributes = ''
    for attribute in retained_attributes:
        if attribute in splited_line[-1]:
            attributes = attributes + attribute + attribute_assignment_operator + splited_line[-1][attribute] + Separator
    new_line.append(attributes)
    new_line = '\t'.join(new_line) + '\n'
    return(new_line)

def split_annotation_line(line, Separator, attribute_assignment_operator):
    line = line.strip().split('\t')
    attributes_dict = {}
    attributes_list = line[-1].split(Separator)
    for n in range(len(attributes_list)):
        if attributes_list[n]:
            att = attributes_list[n].split(attribute_assignment_operator)
            attributes_dict[att[0]] = att[1]
    line[-1] = attributes_dict
    return(line)

def unique_intron_id_per_gene(gene_dict, transcript_dict, Separator, attribute_assignment_operator, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings):
    for gene in gene_dict:
        unique_locu_dict = {}
        for transcript in gene_dict[gene]:
            for n in range(len(transcript_dict[transcript])):
                exon = transcript_dict[transcript][n]
                splited_line = split_annotation_line(exon, Separator, attribute_assignment_operator)
                unique_locu = '_'.join((splited_line[0], splited_line[3], splited_line[4], splited_line[6]))
                if exon_id_attribute in splited_line[-1]:
                    if unique_locu in unique_locu_dict:
                        splited_line[-1][exon_id_attribute] = unique_locu_dict[unique_locu]
                    else:
                        unique_locu_dict[unique_locu] = splited_line[-1][exon_id_attribute]
                transcript_dict[transcript][n] = join_splited_annotation_line(splited_line, Separator, attribute_assignment_operator, retained_attributes)
                transcript_dict[transcript][n] = replace_exon_to_intron(transcript_dict[transcript][n], replace_strings)
    return(transcript_dict)


def generate_intron_attributes(exon, Separator, attribute_assignment_operator, node_connection_keyword, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings):
    exon = exon.strip().split('\t')
    atts = exon[-1].split(Separator)
    atts_dict = {}
    for n in range(len(atts)):
        att = atts[n].split(attribute_assignment_operator)
        atts_dict[att[0]] = att[1]
    if (node_connection_keyword in atts_dict) and (exon_rank_attribute in atts_dict):
        if exon_id_attribute in atts_dict:
            atts_dict[exon_id_attribute] = atts_dict[node_connection_keyword] + '.' + replace_strings[0] + atts_dict[exon_rank_attribute]
        # 'Name' 最好从参数处输入
        elif 'Name' in atts_dict:
            atts_dict['Name'] = atts_dict[node_connection_keyword] + '.' + replace_strings[0] + atts_dict[exon_rank_attribute]
    atts_strs = ''
    for att in retained_attributes:
        atts_strs = atts_strs + att + attribute_assignment_operator + atts_dict[att] + Separator
    exon[-1] = atts_strs
    exon = '\t'.join(exon) + '\n'
    return(exon)

def assignment_intron_attributes_per_transcript(transcript_list, Separator, attribute_assignment_operator, node_connection_keyword, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings):
    for n in range(len(transcript_list)):
        transcript_list[n] = generate_intron_attributes(transcript_list[n], Separator, attribute_assignment_operator, node_connection_keyword, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings)
    return(transcript_list)


def generate_intron_locu(transcript_list):
    for n in range(len(transcript_list)):
        transcript_list[n] = transcript_list[n].strip().split('\t')
    for n in range(len(transcript_list) - 1):
        transcript_list[n][3] = str(int(transcript_list[n][4]) + 1)
        transcript_list[n][4] = str(int(transcript_list[n + 1][3]) - 1)
    for n in range(len(transcript_list)):
        transcript_list[n] = '\t'.join(transcript_list[n]) + '\n'
    return(transcript_list[:-1])

def reverse_negative_strand_transcript(tmp_list, column, negative_symbol, column_separator):
    is_negative = False
    for element in tmp_list:
        ele_symbol = element.split(column_separator)[column]
        is_negative = is_negative or ele_symbol == negative_symbol
    if is_negative:
        tmp_list = tmp_list[::-1]
    return(tmp_list)

def generate_intron_per_transcript(transcript_list):
    transcript_list = reverse_negative_strand_transcript(transcript_list, 6, '-', '\t')
    transcript_list = generate_intron_locu(transcript_list)
    transcript_list = reverse_negative_strand_transcript(transcript_list, 6, '-', '\t')
    return(transcript_list)

def generate_intron_region(gene_dict, transcript_dict, Separator, attribute_assignment_operator, node_connection_keyword, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings):
    for transcript_id in transcript_dict:
        transcript_dict[transcript_id] = generate_intron_per_transcript(transcript_dict[transcript_id])
        transcript_dict[transcript_id] = assignment_intron_attributes_per_transcript(transcript_dict[transcript_id], Separator, attribute_assignment_operator, node_connection_keyword, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings)
    transcript_dict = unique_intron_id_per_gene(gene_dict, transcript_dict, Separator, attribute_assignment_operator, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings)
    return(gene_dict, transcript_dict)


def add_affiliation_into_dict(gene_affiliation, gene_dict, transcript_dict, gene_order, transcript_id, line):
    if transcript_id in gene_affiliation:
        gene_id = gene_affiliation[transcript_id]
        if gene_id in gene_dict:
            if transcript_id not in gene_dict[gene_id]:
                gene_dict[gene_id].append(transcript_id)
        else:
            gene_dict[gene_id] = [transcript_id]
            gene_order.append(gene_id)
        if transcript_id in transcript_dict:
            transcript_dict[transcript_id].append(line)
        else:
            transcript_dict[transcript_id] = [line]
    else:
        print("Error: " + transcript_id + " in " + line + " don't have appear")
    return(gene_affiliation, gene_dict, transcript_dict, gene_order)
    
def is_exon_annotation(line, exon_keyword):
    is_exon = line.strip().split('\t')[2] == exon_keyword
    if is_exon:
        return(True)
    else:
        return(False)

def build_affiliation_map(gene_affiliation, self_id, parent_id):
    if self_id:
        if  self_id not in gene_affiliation:
            gene_affiliation[self_id] = parent_id
    return(gene_affiliation)

def obtain_id_from_attributes(attributes, root_note_attribute, node_connection_keyword):
    self_id, parent_id = '', ''
    if root_note_attribute in attributes:
        self_id = attributes[root_note_attribute]
    if node_connection_keyword in attributes:
        parent_id = attributes[node_connection_keyword]
    return(self_id, parent_id)

def obtain_affiliation(gene_affiliation, gene_dict, transcript_dict, gene_order, line, Separator, attribute_assignment_operator, root_note_attribute, node_connection_keyword, exon_keyword):
    splited_line = split_annotation_line(line, Separator, attribute_assignment_operator)
    self_id, parent_id = obtain_id_from_attributes(splited_line[-1], root_note_attribute, node_connection_keyword)
    gene_affiliation = build_affiliation_map(gene_affiliation, self_id, parent_id)
    if is_exon_annotation(line, exon_keyword):
        gene_affiliation, gene_dict, transcript_dict, gene_order = add_affiliation_into_dict(gene_affiliation, gene_dict, transcript_dict, gene_order, parent_id, line)
    return(gene_affiliation, gene_dict, transcript_dict, gene_order)


def write_into_file(line, output_file):
    with open(output_file, 'w') as o:
        o.writelines(line)

def not_annotation_line(line):
    if line[0] != '@' and line[0] != '#':
        return(True)
    else:
        return(False)

def split_parms(parameters):
    for parm in parameters:
        if ',' in parameters[parm]:
            parameters[parm] = parameters[parm].split(',')
    return(parameters)

def format_parameters():
    usage_str = """Usage: python generate_intron_annotation_from_gff3.py [options]
    options:
    -i [str]  the gene annotation file name with gff3 format.
    -o [str]  the intron annotation file name with gff3 format.
    -e [str]  the Exon feature name in column 3 of gff3 (default: exon).
    -s [str]  the Separator of attributes in column 9 of gff3 (default: ;).
    -a [str]  the Assignment operator of attributes in column 9 of gff3 (default: =).
    -r [str]  the Root node of attribute name in column 9 of gff3 (default: ID).
    -c [str]  the node Connection keyword of attribute name in column 9 of gff3 (default: Parent).
    -E [str]  the Exon id of attribute name in column 9 of gff3 (default: exon_id).
    -R [str]  the exon Rank of attribute name in column 9 of gff3 (default: rank).
    -n [str]  the Names of attribute name that needs to be retained in column 9 of gff3 (default: Parent,Name,constitutive,ensembl_end_phase,ensembl_phase,exon_id,rank).
    -S [str]  the Strings needs to replace in gff3 (default 'exon' to be replace by 'intron': exon,intron).
    """

    input_parms = {'i':['input_file', ''],
        'o':['output_file', ''],
        'e':['exon_keyword', 'exon'],
        's':['Separator', ';'],
        'a':['attribute_assignment_operator', '='],
        'r':['root_note_attribute', 'ID'],
        'c':['node_connection_keyword', 'Parent'],
        'E':['exon_id_attribute', 'exon_id'],
        'R':['exon_rank_attribute', 'rank'],
        'n':['retained_attributes', 'Parent,Name,constitutive,ensembl_end_phase,ensembl_phase,exon_id,rank'],
        'S':['replace_strings', 'exon,intron']
        }

    # no_input_parms is dictory like input_parms, and set the value of index 1 as '' or False.
    no_input_parms = {}

    parameters = obtain_parameters(usage_str, input_parms, no_input_parms)
    return(parameters)

if __name__ == '__main__':
    parameters = format_parameters()
    parameters = split_parms(parameters)
    for parm in parameters:
        locals()[parm] = parameters[parm]

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)

    gene_affiliation, gene_dict, transcript_dict, gene_order = {}, {}, {}, []
    with open(input_file, 'r') as f:
        for line in f:
            if not_annotation_line(line):
                gene_affiliation, gene_dict, transcript_dict, gene_order = obtain_affiliation(gene_affiliation, gene_dict, transcript_dict, gene_order, line, Separator, attribute_assignment_operator, root_note_attribute, node_connection_keyword, exon_keyword)
            else:
                write_into_file(line, output_file)
        gene_dict, transcript_dict = generate_intron_region(gene_dict, transcript_dict, Separator, attribute_assignment_operator, node_connection_keyword, exon_id_attribute, exon_rank_attribute, retained_attributes, replace_strings)
        write_gene_into_file(gene_dict, transcript_dict, gene_order, output_file)


    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))

