import sys, datetime
from getParm import obtain_parameters

def format_parameters():
    usage_str = """Usage: python self_program_name.py [options]
    options:
    -i [str]  the input file name.
    -o [str]  the output file name.
    -m    output interMediate file or not (default: False).
    """

    input_parms = {'i':['input_file', ''],
        'o':['output_file', ''],
        }

    # no_input_parms is dictory like input_parms, and set the value of index 1 as '' or False.
    no_input_parms = {'m':['intermediate', 'False']}

    parameters = obtain_parameters(usage_str, input_parms, no_input_parms)
    return(parameters)

if __name__ == '__main__':
    startTime = datetime.datetime.now()
    print('Start Time:', startTime)

    parameters = format_parameters()

    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))
