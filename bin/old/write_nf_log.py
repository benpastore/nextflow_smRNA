
import argparse

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", required = True )
    
    return parser.parse_args()


def summary(input) : 

    output = open("NF-run-info.log", 'w')
    file = ""
    for arg in input : 
        file += f"{arg}\n"
    
    output.write(file)
    output.close()

def main() : 

    args = get_args() 
    x = args.input.replace('"', "").replace("[", "").replace("]", "")
    input = [ i for i in x.split(",")]
    print(input)
    summary(input)

if __name__ == "__main__" : 

    main()
