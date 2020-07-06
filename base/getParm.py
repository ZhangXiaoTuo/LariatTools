import getopt, sys

def print_parms(key, parm_dict, parms):
    if parm_dict[key]:
        if len(parm_dict[key]) > 1:
            parms[parm_dict[key][0]] = parm_dict[key][1]
            print('-' + key + ': ' + parm_dict[key][0] + ' is ' + parm_dict[key][1])
        else:
            print(key, ':', parm_dict[key])
    return(parms)

def merge_parms(input_parms, no_input_parms):
    parms = {}
    for parm in input_parms:
        parms = print_parms(parm, input_parms, parms)
    for parm in no_input_parms:
        parms = print_parms(parm, no_input_parms, parms)
    return(parms)

def usage(usage_str):
    print(usage_str)
    sys.exit()

def get_parms_str(input_parms, no_input_parms):
    parms_str = 'h'
    for parm in no_input_parms:
        parms_str = parms_str + parm
    for parm in input_parms:
        parms_str = parms_str + parm + ':'
    return(parms_str)

def add_help_parm(usage_str):
    help_str = '-h    print this page.'
    usage_str = usage_str + help_str
    return(usage_str)


def obtain_parameters(usage_str, input_parms, no_input_parms, *key):
    usage_str = add_help_parm(usage_str)
    parms_str = get_parms_str(input_parms, no_input_parms)
    opts, args = getopt.getopt(sys.argv[1:], parms_str)
    if opts:
        for op, value in opts:
            for input_parm in input_parms:
                tmp_parm = '-' + input_parm
                if op == tmp_parm : input_parms[input_parm][1] = value
            for no_input_parm in no_input_parms:
                tmp_parm = '-' + no_input_parm
                print(no_input_parms)
                if op == tmp_parm: no_input_parms[no_input_parm][1] = True
            if op == "-h": usage(usage_str)
    else:
        usage(usage_str)
    for input_parm in input_parms:
        if not input_parms[input_parm][1]:
            print("Error: please input the parameter '-" + input_parm + "'")
            usage(usage_str)
    parms = merge_parms(input_parms, no_input_parms)
    return(parms)
