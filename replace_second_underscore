import sys

def replace_second_underscore(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    for line in lines:
        parts = line.split('_', 2)
        if len(parts) > 2:
            line = parts[0] + '_' + parts[1] + ' ' + parts[2]
        sys.stdout.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: python replace_second_underscore.py <input_file>"
        sys.exit(1)

    input_file = sys.argv[1]
    replace_second_underscore(input_file)
