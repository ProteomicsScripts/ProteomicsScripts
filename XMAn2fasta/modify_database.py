import sys
"""Modify description lines of fasta file.

example usage: python modify_database.py example_file.fasta
The above command will create a new file with the name 'example_file_mod.fasta"""

def modify(in_):
    """Modify the description line which does not contain ' OR '."""
    tokens_pipe = in_.split("|")
    name = tokens_pipe[0].split(",")
    mutation = tokens_pipe[1]
    id_ = name[0]
    
    new_line = id_.split("_")[0] + "_" + id_.split("_")[1]
    new_line += "|" + id_.split("_")[-1] + "_" + mutation + "|"
    
    if len(name) > 1:
        for k in name[1:]:
            new_line += k + ","
    
    for i in tokens_pipe[2:]:
        new_line = new_line + "," + i
    
    return new_line

def modify_or_line(line):
    """Modify description line which contains ' OR '."""
    tokens_pipe = line.split("|")
    name = tokens_pipe[1]
    mutation = tokens_pipe[2]
    new_line = ">" + tokens_pipe[0] + "|"

    if "," in name:
        name_tokens = name.split(",")
        new_line = new_line + name_tokens[0] + "_"
        new_line = new_line + mutation + "|"
        for k in name_tokens[1:]:
            new_line += k + ","

    for i in tokens_pipe[3:]:
        new_line += i + ","

    new_line = new_line[:-1]

    return new_line

filename = sys.argv[1]    # Input fasta file

outfile = filename.split(".")[0]
outfile += "_mod.fasta"    # Output fasta file

f_in = open(filename)
lines = f_in.readlines()
f_in.close()

f_out = open(outfile, 'w')

for i in range(0, len(lines), 2):
    # assumption in for-loop: every second line in fasta file starts with '>'
    if not lines[i].startswith(">"):
        raise AttributeError("line", i+1, "does not start with '>'")
    
    if " OR " not in lines[i]:
        new_line = modify(lines[i])
        
        f_out.write(new_line)
        f_out.write(lines[i+1])
    else:
        # if description line contains multiple descriptions seperated by " OR ",
        # save same sequence with multiple modified descriptions
        multi_desc = lines[i].split(" OR ")
        
        new_line = modify(multi_desc[0])
        new_line = new_line + "\r\n"
        f_out.write(new_line)
        f_out.write(lines[i+1])
        
        for j in range(1,len(multi_desc)):
            line = multi_desc[j].split("\r\n")[0]
            new_line = modify_or_line(line)
            new_line = new_line + "\r\n"
            f_out.write(new_line)
            f_out.write(lines[i+1])
        
    
f_out.close()
