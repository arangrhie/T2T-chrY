import sys
import re

def convert_string(input_string):
    items =   re.findall(r'([><\[][^\]><[]*\]?)', input_string)
    # Process each item
    result = []
    for item in items:
        if item.startswith('>'):
            result.append(item.lstrip('<>') + '+')
        elif item.startswith('<'):
            result.append(item.lstrip('<>') + '-')
        elif item.endswith(']'):
            result.append(item)
        else :
            error = "Invalid input string: " + item
            print(error)
            return

    
    # Join the processed items into a single string
    return ','.join(result)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python convert_path.py consensus_paths.txt")
        print("consensus_paths.txt Verkko output in 6-layoutContigs")
        sys.exit(1)

    input_file = sys.argv[1]
    with open(input_file, 'r') as fi:
        for line in fi:
            # print(f":: DEBUG :: {line.strip()}")
            tokens = line.strip()
            tokens = tokens.split('\t')
            # Skip the first line
            if tokens[0] == "name":
                print(f"{line.strip()}")
                continue
            path=tokens[1]
            output = convert_string(path)
            print(f"{tokens[0]}\t{output}\t{tokens[2]}")
            

